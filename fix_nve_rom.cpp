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
#include "force.h"
#include "respa.h"
#include "update.h"

// ****************** ADDED *********************

#include "domain.h"
#include "memory.h"
#include "text_file_reader.h"
#include "error.h"
#include <iostream>
#include <fstream>
#include <sstream>

// *************** END ADDED *******************

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEROM::FixNVEROM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) utils::missing_cmd_args(FLERR, "fix nve/rom", error);

  dynamic_group_allow = 1;
  time_integrate = 1;

  // ****************** ADDED *********************

  modelorder = utils::inumeric(FLERR,arg[3],false,lmp);

  int nlocal = atom->nlocal;

  phi = nullptr;
  mean = nullptr;
  A = nullptr;
  V = nullptr;
  X = nullptr;

  memory->create(phi, nlocal * 3, modelorder, "FixNVEROM:phi");
  memory->create(mean, nlocal * 3, 1, "FixNVEROM:mean");
  memory->create(A, nlocal * 3, "FixNVEROM:A");
  memory->create(V, nlocal * 3, "FixNVEROM:V");
  memory->create(X, nlocal * 3, "FixNVEROM:X");
  
  read_rob(arg[4]);
  read_mean(arg[5]);

  // *************** END ADDED *******************
}

/* ---------------------------------------------------------------------- */

int FixNVEROM::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVEROM::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (utils::strmatch(update->integrate_style,"^respa"))
    step_respa = (dynamic_cast<Respa *>(update->integrate))->step;
  
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEROM::initial_integrate(int /*vflag*/)
{
  double dtfm;

  // update v and x of atoms in group
  
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  
  //******************* ADDED ******************

  double boxlo[3], boxhi[3];
  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];
  boxhi[0] = domain->boxhi[0];
  boxhi[1] = domain->boxhi[1];
  boxhi[2] = domain->boxhi[2];

  double length_x = boxhi[0]-boxlo[0];
  double length_y = boxhi[1]-boxlo[1];
  double length_z = boxhi[2]-boxlo[2];

  remove_periodic_boundary(length_x, length_y, length_z);

  convert_physical_to_reduced();

  //**************** END ADDED *****************

  // This section is heavily edited. Most atom tests have been moved to the convert functions

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < modelorder; i++) // dimension is modelorder
    if (mask[0] & groupbit) { // only works with mask[0]
      V[i] += dtf * A[i]; // changed from dtfm
      X[i] += dtv * V[i];
    }
  
  //******************* ADDED ******************

  convert_reduced_to_physical();
  reinstate_periodic_boundary(length_x, length_y, length_z);

  //**************** END ADDED *****************
}

/* ---------------------------------------------------------------------- */

void FixNVEROM::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;

  //******************* ADDED ******************

  convert_physical_to_reduced();

  //**************** END ADDED *****************

  // This section is heavily edited. Most atom tests have been moved to the convert functions

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < modelorder; i++) // dimension is modelorder
    if (mask[0] & groupbit) { // only works with mask[0]
      V[i] += dtf * A[i]; // changed from dtfm
    }

  //******************* ADDED ******************

  convert_reduced_to_physical();

  //**************** END ADDED *****************
}

/* ---------------------------------------------------------------------- */

void FixNVEROM::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEROM::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEROM::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- 
    Consider combining the file reading functions into one
   ---------------------------------------------------------------------- */

void FixNVEROM::read_rob(std::string robfile)
{
  utils::logmesg(lmp, "Reading reduced-order basis file {}\n", robfile);
  const int nlocal = atom->nlocal;

  // try {
  //   TextFileReader reader(file, "reduced-order basis");
  //   for (int i = 0; i < nlocal * 3; i++) {
  //     //reader.next_dvector(phi[i], modelorder);
  //     ValueTokenizer values = reader.next_values(modelorder);
  //     for (int j = 0; j < modelorder; j++) {
  //       phi[i][j] = values.next_double();
  //     }
  //   }
  // } catch (std::exception &e) {
  //   error->one(FLERR, "Error reading file with reduced-order basis: {}", e.what());
  // }

  try {
    std::ifstream file(robfile);
    std::string line;
    if (file.is_open()) {
      for (int i = 0; i < nlocal * 3; i++) {
        std::getline(file, line);
        std::stringstream ss(line);
        for (int j = 0; j < modelorder; j++) {
          ss >> phi[i][j];
        }
      }
    }
  } catch (std::exception &e) {
    error->one(FLERR, "Error reading reduced-order basis: {}", e.what());
  }
}

/* ---------------------------------------------------------------------- 
    Consider formating the mean matrix as a nlocal x 3 matrix
   ---------------------------------------------------------------------- */

void FixNVEROM::read_mean(std::string file)
{
  utils::logmesg(lmp, "Reading ROM means file {}\n", file);
  const int nlocal = atom->nlocal;

  try {
    TextFileReader reader(file, "ROM means");
    for (int i = 0; i < nlocal * 3; i++) {
      reader.next_dvector(mean[i], 1);
    }
  } catch (std::exception &e) {
    error->one(FLERR, "Error reading file with ROM means: {}", e.what());
  }
}

/* ---------------------------------------------------------------------- 
    Code from Hao Zhang. Consider:
      - using in-built LAMMPS code to toggle periodic boundaries
*/

void FixNVEROM::remove_periodic_boundary(double length_x, double length_y, double length_z)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  imageint *image = atom->image;
  imageint ximage, yimage, zimage;

  for (int i = 0; i < nlocal; i++){
    ximage = image[i] & IMGMASK;
    yimage = (image[i] >> IMGBITS) & IMGMASK;
    zimage = image[i] >> IMG2BITS;

    x[i][0] = x[i][0] + length_x * (ximage-512);
    x[i][1] = x[i][1] + length_y * (yimage-512);
    x[i][2] = x[i][2] + length_z * (zimage-512);
  }
}

/* ---------------------------------------------------------------------- 
    Code from Hao Zhang. Consider:
      - using in-built LAMMPS code to toggle periodic boundaries
      - try not to repeat this code from the previous function
*/

void FixNVEROM::reinstate_periodic_boundary(double length_x, double length_y, double length_z)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  imageint *image = atom->image;
  imageint ximage, yimage, zimage;
  
  for (int i=0; i<nlocal; i++){
      ximage = image[i] & IMGMASK;
      yimage = (image[i] >> IMGBITS) & IMGMASK;
      zimage = image[i] >> IMG2BITS;

      x[i][0] = x[i][0] - length_x * (ximage-512);
      x[i][1] = x[i][1] - length_y * (yimage-512);
      x[i][2] = x[i][2] - length_z * (zimage-512);
  }
}

/* ---------------------------------------------------------------------- 
    Code from Hao Zhang. Converts from physical to reduced-order space.
    Assumes that tensors for reduced order X, V, and F are defined.
*/

void FixNVEROM::convert_physical_to_reduced()
{
  
  int atomid;
  
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom-> type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (int i = 0; i < modelorder; i++){
    X[i] = 0.0;
    V[i] = 0.0;
    A[i] = 0.0;

    if (rmass) {
      for (int j = 0; j < nlocal; j++){
        atomid = tag[j] - 1;
        X[i] += phi[atomid][i]* (x[j][0] - mean[atomid][0]) + phi[atomid+nlocal][i]* (x[j][1] - mean[atomid+nlocal][0]) + phi[atomid+nlocal*2][i]* (x[j][2] - mean[atomid+nlocal*2][0]);
        V[i] += phi[atomid][i]*v[j][0] + phi[atomid+nlocal][i]*v[j][1] + phi[atomid+nlocal*2][i]*v[j][2];
        A[i] += (phi[atomid][i]*f[j][0] + phi[atomid+nlocal][i]*f[j][1] + phi[atomid+nlocal*2][i]*f[j][2]) / rmass[j];                 
      }
    }
    else {
      for (int j = 0; j < nlocal; j++){
        atomid = tag[j] - 1;
        X[i] += phi[atomid][i]* (x[j][0] - mean[atomid][0]) + phi[atomid+nlocal][i]* (x[j][1] - mean[atomid+nlocal][0]) + phi[atomid+nlocal*2][i]* (x[j][2] - mean[atomid+nlocal*2][0]);
        V[i] += phi[atomid][i]*v[j][0] + phi[atomid+nlocal][i]*v[j][1] + phi[atomid+nlocal*2][i]*v[j][2];
        A[i] += (phi[atomid][i]*f[j][0] + phi[atomid+nlocal][i]*f[j][1] + phi[atomid+nlocal*2][i]*f[j][2]) / mass[type[j]];                 
      }
    }

  }
}

/* ---------------------------------------------------------------------- 
    Code from Hao Zhang. Converts from physical to reduced-order space.
    Assumes that tensors for reduced order X, V, and F are defined.
*/

void FixNVEROM::convert_reduced_to_physical()
{
  double **x = atom->x;
  double **v = atom->v;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (int i=0; i<nlocal; i++){
    int atomid = tag[i] - 1;

    x[i][0] = mean[atomid][0];
    x[i][1] = mean[atomid+nlocal][0];
    x[i][2] = mean[atomid+2*nlocal][0];

    v[i][0] = 0.0;
    v[i][1] = 0.0;
    v[i][2] = 0.0;

    for (int j=0; j<modelorder; j++){
      x[i][0] += phi[atomid][j] * X[j];
      x[i][1] += phi[atomid+nlocal][j] * X[j];
      x[i][2] += phi[atomid+nlocal*2][j] * X[j];

      v[i][0] += phi[atomid][j] * V[j];
      v[i][1] += phi[atomid+nlocal][j] * V[j];
      v[i][2] += phi[atomid+nlocal*2][j] * V[j];
    }
  }
}