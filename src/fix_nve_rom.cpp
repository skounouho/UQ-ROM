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

#include "fix_nve_rom.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "respa.h"
#include "update.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "eigen/Eigen/Eigen"
using namespace Eigen;

using namespace LAMMPS_NS;
using namespace FixConst;

// --- TESTING PARALLELISM ----
#include <omp.h>

/* ---------------------------------------------------------------------- */

FixNVEROM::FixNVEROM(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg),
  phi(nullptr), y(nullptr), y_dot(nullptr), y_dot_dot(nullptr), x0(nullptr)
  // y_all(nullptr), y_dot_all(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix nve/rom", error);

  dynamic_group_allow = 1;
  time_integrate = 1;
  vector_flag = 1;
  extvector = 0;

  modelorder = utils::inumeric(FLERR,arg[3],false,lmp);

  // create arrays

  int nlocal = atom->nlocal;

  memory->create(phi, nlocal * 3, modelorder, "FixNVEROM:phi");
  memory->create(y, modelorder, "FixNVEROM:y");
  memory->create(y_dot, modelorder, "FixNVEROM:y_dot");
  memory->create(y_dot_dot, modelorder, "FixNVEROM:y_dot_dot");
  memory->create(x0, nlocal, 3, "FixNVEROM:x0");

  // ----- TESTING MASS MATRIX ------
  memory->create(F, modelorder, "FixNVEROM:redm"); // reduced f
  memory->create(M, modelorder, modelorder, "FixNVEROM:redm"); // reduced m

  // set compute variable size
  size_vector = modelorder;

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
  }

  memory->destroy(unwrap);

  // check if rob file is available and readable
  if (!platform::file_is_readable(arg[4]))
    error->all(FLERR, fmt::format("Cannot open file {}: {}", arg[4], utils::getsyserror()));
  
  read_rob(arg[4], phi);

  // check if phi matrix is filled and not NaN
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
}

/* ---------------------------------------------------------------------- */

Eigen::MatrixXd FixNVEROM::convert_to_matrix(double **data, int rows, int columns)
{
    Eigen::MatrixXd eMatrix(rows, columns);
    for (int i = 0; i < rows; ++i)
        eMatrix.row(i) = Eigen::VectorXd::Map(&data[i][0], columns);
    return eMatrix;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEROM::initial_integrate(int /*vflag*/)
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

  for (int i = 0; i < modelorder; i++) {
    y_dot[i] += dtf * y_dot_dot[i];
    y[i] += dtv * y_dot[i];
  }

  // sum across processors
  // MPI_Allreduce(y, y_all, modelorder, MPI_DOUBLE, MPI_SUM, world);
  // MPI_Allreduce(y_dot, y_dot_all, modelorder, MPI_DOUBLE, MPI_SUM, world);
  // y = y_all;
  // y_dot = y_dot_all;

  update_physical_variables(xflag);

  // wrap atom positions into simulation box
  if (xflag)
    for (int i = 0; i < nlocal; i++) {
      domain->remap(x[i], image[i]);
    }

}

/* ---------------------------------------------------------------------- */

void FixNVEROM::final_integrate()
{
  int xflag = 0;

  compute_reduced_variables(xflag);

  // This section is heavily edited. Most atom group tests have been moved to the convert functions

  // ---- TEST MASS CALCULATING y_dot_dot ----
  VectorXd b(modelorder);
  for (int i = 0; i < modelorder; i++) {
    b(i) = F[i];
  }

  VectorXd soln = A.colPivHouseholderQr().solve(b);
  for (int i = 0; i < modelorder; i++) {
    y_dot_dot[i] = soln(i);
  }

  for (int i = 0; i < modelorder; i++) {
    y_dot[i] += dtf * y_dot_dot[i];
  }

  // sum across processors
  // MPI_Allreduce(y_dot, y_dot_all, modelorder, MPI_DOUBLE, MPI_SUM, world);
  // y_dot = y_dot_all;

  update_physical_variables(xflag);

}

/* ----------------------------------------------------------------------
    return a single element of y[modelorder]
   ---------------------------------------------------------------------- */

double FixNVEROM::compute_vector(int n)
{
  return y[n];
}

/* ---------------------------------------------------------------------- 
    Reads the reduced-order basis from a designated file.
   ---------------------------------------------------------------------- */

void FixNVEROM::read_rob(std::string robfile, double **robarray)
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

void FixNVEROM::compute_reduced_variables(int xflag)
{
  int i,j,iatom;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom-> type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  #pragma omp parallel for
  for (j = 0; j < modelorder; j++){
    if (xflag) y[j] = 0;
    y_dot[j] = 0;
    F[j] = 0.0;

    for (i = 0; i < nlocal; i++) {
      iatom = tag[i] - 1;

      if (xflag) {
        y[j] += phi[iatom][j]                * (x[i][0] - x0[iatom][0]);
        y[j] += phi[iatom + nlocal][j]       * (x[i][1] - x0[iatom][1]);
        y[j] += phi[iatom + nlocal*2][j]     * (x[i][2] - x0[iatom][2]);
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

void FixNVEROM::update_physical_variables(int xflag)
{
  int i,j,iatom;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  #pragma omp parallel for
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      iatom = tag[i] - 1;

      if (xflag) {
        x[i][0] = x0[iatom][0];
        x[i][1] = x0[iatom][1];
        x[i][2] = x0[iatom][2];
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
