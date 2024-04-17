// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#include "fix_nh_rom.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "force.h"
#include "group.h"
#include "irregular.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <fstream>
#include <sstream>

#include "eigen/Eigen/Eigen"
using namespace Eigen;

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double DELTAFLIP = 0.1;
static constexpr double TILTMAX = 1.5;
static constexpr double EPSILON = 1.0e-6;

enum{NOBIAS,BIAS};
enum{NONE,XYZ,XY,YZ,XZ};
enum{ISO,ANISO,TRICLINIC};

// --- TESTING PARALLELISM ----
#include <omp.h>

/* ---------------------------------------------------------------------- */

FixNHROM::FixNHROM(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg),
  phi(nullptr), y(nullptr), y_dot(nullptr), y_dot_dot(nullptr), x0(nullptr)
  // y_all(nullptr), y_dot_all(nullptr)
{
  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix nvt/nph/npt rom", error);

      dynamic_group_allow = 1;
      time_integrate = 1;

      modelorder = utils::inumeric(FLERR,arg[iarg + 1],false,lmp);
      
      // set up MPI

      // MPI_Comm_rank(world,&me);
      // MPI_Comm_size(world,&nprocs);

      // allocate arrays

      int nlocal = atom->nlocal;

      memory->create(phi, nlocal * 3, modelorder, "FixNHROM:phi");
      memory->create(y, modelorder, "FixNHROM:y");
      memory->create(y_dot, modelorder, "FixNHROM:y_dot");
      memory->create(y_dot_dot, modelorder, "FixNHROM:y_dot_dot");
      memory->create(x0, nlocal, 3, "FixNHROM:x0");

      // ----- TESTING MASS MATRIX ------
      memory->create(F, modelorder, "FixNVEROM:redm"); // reduced f
      memory->create(M, modelorder, modelorder, "FixNVEROM:redm"); // reduced m

      // ---- TESTING ZETA -----
      memory->create(lamda0, nlocal, 3, "FixNHROM:lamda0");

      // save initial atom positions

      double **x = atom->x;
      imageint *image = atom->image;
      int *tag = atom->tag;
      int iatom;
      
      // unmap atoms
      double **unwrap;
      memory->create(unwrap,atom->nlocal,3,"rob:unwrap");

      for (int i = 0; i < nlocal; i++) domain->unmap(x[i], image[i], unwrap[i]);

      for (int i = 0; i < nlocal; i++) {
        iatom = tag[i] - 1;
        x0[iatom][0] = unwrap[i][0];
        x0[iatom][1] = unwrap[i][1];
        x0[iatom][2] = unwrap[i][2];

        domain->x2lamda(x0[iatom], lamda0[iatom]);
      }

      memory->destroy(unwrap);
      
      // check if rob file is available and readable
      if (!platform::file_is_readable(arg[iarg + 2]))
        error->all(FLERR, fmt::format("Cannot open file {}: {}", arg[iarg + 2], utils::getsyserror()));

      read_rob(arg[iarg + 2], phi);

      // check if last entry is NaN
      if (phi[nlocal * 3 - 1][modelorder - 1] != phi[nlocal * 3 - 1][modelorder - 1]) {
        error->all(FLERR, "Reduced order basis is not filled");
      }
       // ---- ADD PROCEDURE FOR ASSEMBLING MASS MATRIX ---
      double *rmass = atom->rmass;
      double *mass = atom->mass;
      int *type = atom->type;
      double **MPhi;
      memory->create(MPhi, nlocal * 3, modelorder, "nve/rom:MPhi");

      for (int i=0; i < nlocal; i++) {
        iatom = tag[i] - 1;
        double m;
        if (rmass) m = rmass[i]; else m = mass[type[i]];
        for (int j=0; j < modelorder; j++) {
          MPhi[iatom][j] = m * phi[iatom][j];
          MPhi[iatom + nlocal][j] = m * phi[iatom + nlocal][j];
          MPhi[iatom + nlocal*2][j] = m * phi[iatom + nlocal*2][j];
        }
      }

      for (int i=0; i < modelorder; i++) {
        for (int j=0; j < modelorder; j++) {
          M[i][j] = 0;
          for (int k=0; k < nlocal*3; k++) {
            M[i][j] += phi[k][i] * MPhi[k][j];
          }
        }
      }

      A = convert_to_matrix(M, modelorder, modelorder);
      memory->destroy(M);

      iarg += 3;
    } else iarg++;
  }

  if (pstat_flag) {
    // thermo_virial = 1;
    // virial_global_flag = 1;
    // vflag_global = 1;
  }
}

int FixNHROM::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  mask |= PRE_FORCE_RESPA;
  if (pre_exchange_flag) mask |= PRE_EXCHANGE;
  return mask;
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts
------------------------------------------------------------------------- */

void FixNHROM::setup(int /*vflag*/)
{
  // tdof needed by compute_temp_target()

  t_current = temperature->compute_scalar();
  tdof = temperature->dof;

  // t_target is needed by NVT and NPT in compute_scalar()
  // If no thermostat or using fix nphug,
  // t_target must be defined by other means.

  if (tstat_flag && strstr(style,"nphug") == nullptr) {
    compute_temp_target();
  } else if (pstat_flag) {

    // t0 = reference temperature for masses
    // set equal to either ptemp or the current temperature
    // cannot be done in init() b/c temperature cannot be called there
    // is b/c Modify::init() inits computes after fixes due to dof dependence
    // error if T less than 1e-6
    // if it was read in from a restart file, leave it be

    if (t0 == 0.0) {
      if (p_temp_flag) {
        t0 = p_temp;
      } else {
        t0 = temperature->compute_scalar();
        if (t0 < EPSILON)
          error->all(FLERR,"Current temperature too close to zero, consider using ptemp keyword");
      }
    }
    t_target = t0;
  }

  if (pstat_flag) compute_press_target();

  if (pstat_flag) {
    zeta_virial();
    if (pstyle == ISO) pressure->compute_scalar();
    else pressure->compute_vector();
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  // masses and initial forces on thermostat variables

  if (tstat_flag) {
    eta_mass[0] = tdof * boltz * t_target / (t_freq*t_freq);
    for (int ich = 1; ich < mtchain; ich++)
      eta_mass[ich] = boltz * t_target / (t_freq*t_freq);
    for (int ich = 1; ich < mtchain; ich++) {
      eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1] -
                         boltz * t_target) / eta_mass[ich];
    }
  }

  // masses and initial forces on barostat variables

  if (pstat_flag) {
    double kt = boltz * t_target;
    double nkt = (atom->natoms + 1) * kt;

    for (int i = 0; i < 3; i++)
      if (p_flag[i])
        omega_mass[i] = nkt/(p_freq[i]*p_freq[i]);

    if (pstyle == TRICLINIC) {
      for (int i = 3; i < 6; i++)
        if (p_flag[i]) omega_mass[i] = nkt/(p_freq[i]*p_freq[i]);
    }

  // masses and initial forces on barostat thermostat variables

    if (mpchain) {
      etap_mass[0] = boltz * t_target / (p_freq_max*p_freq_max);
      for (int ich = 1; ich < mpchain; ich++)
        etap_mass[ich] = boltz * t_target / (p_freq_max*p_freq_max);
      for (int ich = 1; ich < mpchain; ich++)
        etap_dotdot[ich] =
          (etap_mass[ich-1]*etap_dot[ich-1]*etap_dot[ich-1] -
           boltz * t_target) / etap_mass[ich];
    }
  }
}

/* ----------------------------------------------------------------------
   1st half of Verlet update
------------------------------------------------------------------------- */

void FixNHROM::initial_integrate(int /*vflag*/)
{
  // update eta_press_dot

  if (pstat_flag && mpchain) nhc_press_integrate();

  // update eta_dot

  if (tstat_flag) {
    compute_temp_target();
    nhc_temp_integrate();
  }

  // need to recompute pressure to account for change in KE
  // t_current is up-to-date, but compute_temperature is not
  // compute appropriately coupled elements of mvv_current

  if (pstat_flag) {
    zeta_virial();
    if (pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  if (pstat_flag) {
    compute_press_target();
    nh_omega_dot();
    nh_v_press();
  }

  nve_v();

  // remap simulation box by 1/2 step

  if (pstat_flag) remap();

  nve_x();

  // remap simulation box by 1/2 step
  // redo KSpace coeffs since volume has changed

  if (pstat_flag) {
    remap();
    if (kspace_flag) force->kspace->setup();
  }
}

/* ----------------------------------------------------------------------
   2nd half of Verlet update
------------------------------------------------------------------------- */

void FixNHROM::final_integrate()
{
  nve_v();

  // re-compute temp before nh_v_press()
  // only needed for temperature computes with BIAS on reneighboring steps:
  //   b/c some biases store per-atom values (e.g. temp/profile)
  //   per-atom values are invalid if reneigh/comm occurred
  //     since temp->compute() in initial_integrate()

  if (which == BIAS && neighbor->ago == 0)
    t_current = temperature->compute_scalar();

  if (pstat_flag) nh_v_press();

  // compute new T,P after velocities rescaled by nh_v_press()
  // compute appropriately coupled elements of mvv_current

  t_current = temperature->compute_scalar();
  tdof = temperature->dof;

  // need to recompute pressure to account for change in KE
  // t_current is up-to-date, but compute_temperature is not
  // compute appropriately coupled elements of mvv_current

  if (pstat_flag) {
    zeta_virial();
    if (pstyle == ISO) pressure->compute_scalar();
    else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  if (pstat_flag) nh_omega_dot();

  // update eta_dot
  // update eta_press_dot

  if (tstat_flag) nhc_temp_integrate();
  if (pstat_flag && mpchain) nhc_press_integrate();
}

/* ---------------------------------------------------------------------- */

Eigen::MatrixXd FixNHROM::convert_to_matrix(double **data, int rows, int columns)
{
    Eigen::MatrixXd eMatrix(rows, columns);
    for (int i = 0; i < rows; ++i)
        eMatrix.row(i) = Eigen::VectorXd::Map(&data[i][0], columns);
    return eMatrix;
}

/*  --------------------------------------------------------------------
    post force step for calculating virial
--------------------------------------------------------------------  */

void FixNHROM::zeta_virial()
{
  int i, j, iatom;
  int nlocal = atom->nlocal;
  int *tag = atom->tag;
  double **f = atom->f;

  // calculate projected force
  for (j = 0; j < modelorder; j++) {
    F[j] = 0;
    for (i = 0; i < nlocal; i++) {
      iatom = tag[i] - 1;
      F[j] += phi[iatom][j]          * f[i][0];
      F[j] += phi[iatom+nlocal][j]   * f[i][1];
      F[j] += phi[iatom+nlocal*2][j] * f[i][2]; 
    }
  }

  // scale positions
  double **x0temp;
  memory->create(x0temp, nlocal, 3, "FixNHROM:x0temp");

  for (i = 0; i < nlocal; i++) {
    iatom = tag[i] - 1;
    domain->lamda2x(lamda0[iatom], x0temp[iatom]);
  }

  // calculate projected position
  double *y0;
  memory->create(y0, modelorder, "FixNHROM:y0");

  for (j = 0; j < modelorder; j++) {
    y0[j] = 0;
    for (i = 0; i < nlocal; i++) {
      iatom = tag[i] - 1;
      y0[j] += phi[iatom][j]          * x0[i][0];
      y0[j] += phi[iatom+nlocal][j]   * x0[i][1];
      y0[j] += phi[iatom+nlocal*2][j] * x0[i][2]; 
    }
  }

  // compute virial
  virial[0] = 0;
  for (j = 0; j < modelorder; j++) {
    virial[0] += y0[j] * F[j];
  }
}

/* ----------------------------------------------------------------------
   perform half-step update of velocities in ROM
-----------------------------------------------------------------------*/

void FixNHROM::nve_v()
{
  int xflag = 0;

  compute_reduced_variables(xflag);

  // ---- TEST MASS CALCULATING y_dot_dot ----
  VectorXd b(modelorder);
  for (int i = 0; i < modelorder; i++) {
    b(i) = F[i];
  }
  VectorXd soln = A.colPivHouseholderQr().solve(b);
  for (int i = 0; i < modelorder; i++) {
    y_dot_dot[i] = soln(i);
  }

  // This section is heavily edited. Most atom group tests have been moved to the convert functions

  for (int i = 0; i < modelorder; i++) {
    y_dot[i] += dtf * y_dot_dot[i]; // changed from dtfm
  }

  // sum across processors
  // MPI_Allreduce(y_dot, y_dot_all, modelorder, MPI_DOUBLE, MPI_SUM, world);
  // y_dot = y_dot_all;

  update_physical_variables(xflag);
}

/* ----------------------------------------------------------------------
   perform full-step update of position with ROM
-----------------------------------------------------------------------*/

void FixNHROM::nve_x()
{
  int xflag = 1;
  imageint *image = atom->image;
  double **x = atom->x;
  double nlocal = atom->nlocal;

  // unwrap atom positions from simulation box and reset image box
  if (xflag) 
    for (int i = 0; i < nlocal; i++) {
      domain->unmap(x[i], image[i]);
      image[i] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
    }

  compute_reduced_variables(xflag);

  // ---- TEST MASS CALCULATING y_dot_dot ----
  VectorXd b(modelorder);
  for (int i = 0; i < modelorder; i++) {
    b(i) = F[i];
  }
  VectorXd soln = A.colPivHouseholderQr().solve(b);
  for (int i = 0; i < modelorder; i++) {
    y_dot_dot[i] = soln(i);
  }

  // This section is heavily edited. Most atom group tests have been moved to the convert functions

  for (int i = 0; i < modelorder; i++) { // dimension is modelorder
    y[i] += dtv * y_dot[i];
  }

  // sum across processors
  // MPI_Allreduce(y, y_all, modelorder, MPI_DOUBLE, MPI_SUM, world);
  // y = y_all;

  update_physical_variables(xflag);

  // wrap atom positions into simulation box
  if (xflag)
    for (int i = 0; i < nlocal; i++) {
      domain->remap(x[i], image[i]);
    }
}

/* ---------------------------------------------------------------------- 
    Reads the reduced-order basis from a designated file.
   ---------------------------------------------------------------------- */

void FixNHROM::read_rob(std::string robfile, double **robarray)
{
  utils::logmesg(lmp, "Reading reduced-order basis file {}\n", robfile);
  const int nlocal = atom->nlocal;

  try {
    std::ifstream file(robfile);
    std::string line;
    if (file.is_open()) {
      for (int i = 0; i < nlocal * 3; i++) {
        std::getline(file, line);
        std::stringstream ss(line);
        for (int j = 0; j < modelorder; j++) {
          ss >> robarray[i][j];
        }
      }
    }
  } catch (std::exception &e) {
    error->one(FLERR, "Error reading reduced-order basis: {}", e.what());
  }
}

/* ---------------------------------------------------------------------- 
    Converts from physical to reduced-order space. Updates y, y_dot, and
    y_dot_dot during initial integration, and y_dot and y_dot_dot during
    final integration.
   ---------------------------------------------------------------------- */

void FixNHROM::compute_reduced_variables(int xflag)
{
  int i,j,iatom;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom-> type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // ---- ZETA SCALING ----
  double zetax0[3];

  #pragma omp parallel for
  for (j = 0; j < modelorder; j++){
    if (xflag) y[j] = 0;
    y_dot[j] = 0;
    y_dot_dot[j] = 0.0;
    F[j] = 0.0;

    for (i = 0; i < nlocal; i++) {
      iatom = tag[i] - 1;

      if (xflag) {
        domain->lamda2x(lamda0[iatom], zetax0);
        for (int d = 0; d < 3; d++) y[j] += phi[iatom + nlocal*d][j] * (x[i][d] - zetax0[d]);
      }  
      
      y_dot[j] += phi[iatom][j]              * v[i][0];
      y_dot[j] += phi[iatom + nlocal][j]     * v[i][1];
      y_dot[j] += phi[iatom + nlocal*2][j]   * v[i][2];

      F[j] += phi[iatom][j]          * f[i][0];
      F[j] += phi[iatom+nlocal][j]   * f[i][1];
      F[j] += phi[iatom+nlocal*2][j] * f[i][2];               
    }
  }  
}

/* ---------------------------------------------------------------------- 
    Converts from reduced-order to physical space. Updates x and v during
    initial integration and only v during final integration.
   ---------------------------------------------------------------------- */

void FixNHROM::update_physical_variables(int xflag)
{
  int i,j,iatom;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // ---- ZETA SCALING ----
  double zetax0[3];

  #pragma omp parallel for
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      iatom = tag[i] - 1;

      if (xflag) {
        domain->lamda2x(lamda0[iatom], zetax0);
        for (int d = 0; d < 3; d++) x[i][d] = zetax0[d];
      }

      v[i][0] = 0.0;
      v[i][1] = 0.0;
      v[i][2] = 0.0;

      for (j = 0; j < modelorder; j++){

        if (xflag) {
          x[i][0] += phi[iatom][j]            * y[j];
          x[i][1] += phi[iatom + nlocal][j]   * y[j];
          x[i][2] += phi[iatom + nlocal*2][j] * y[j];  
        }

        v[i][0] += phi[iatom][j]              * y_dot[j];
        v[i][1] += phi[iatom + nlocal][j]     * y_dot[j];
        v[i][2] += phi[iatom + nlocal*2][j]   * y_dot[j];
      }
    }
}
