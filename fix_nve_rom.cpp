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
#include "error.h"
#include <iostream>
#include <fstream>
#include <sstream>

// *************** END ADDED *******************

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEROM::FixNVEROM(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix nve/rom", error);

  dynamic_group_allow = 1;
  time_integrate = 1;

  modelorder = utils::inumeric(FLERR,arg[3],false,lmp);

  // create arrays

  int nlocal = atom->nlocal;

  phi = nullptr;
  y = nullptr;
  y_dot = nullptr;
  y_dot_dot = nullptr;  
  x0 = nullptr;

  memory->create(phi, nlocal * 3, modelorder, "FixNVEROM:phi");
  memory->create(y_dot_dot, nlocal * 3, "FixNVEROM:y_dot_dot");
  memory->create(y_dot, nlocal * 3, "FixNVEROM:y_dot");
  memory->create(y, nlocal * 3, "FixNVEROM:y");
  memory->create(x0, nlocal, 3, "FixNVEROM:x0");

  // save initial atom positions

  double **xinit = atom->x;
  int *tag = atom->tag;
  int iatom;

  for (int i = 0; i < nlocal; i++) {
    iatom = tag[i] - 1;
    x0[iatom][0] = xinit[i][0];
    x0[iatom][1] = xinit[i][1];
    x0[iatom][2] = xinit[i][2];
  }
  
  read_rob(arg[4], phi);
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEROM::initial_integrate(int /*vflag*/)
{
  int xflag = 1;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  
  compute_reduced_variables(xflag);

  // This section is heavily edited. Most atom tests have been moved to the convert functions

  for (int i = 0; i < modelorder; i++) // dimension is modelorder
    if (mask[0] & groupbit) { // only works with mask[0]
      y_dot[i] += dtf * y_dot_dot[i]; // changed from dtfm
      y[i] += dtv * y_dot[i];
    }

  update_physical_variables(xflag);
}

/* ---------------------------------------------------------------------- */

void FixNVEROM::final_integrate()
{
  int xflag = 0;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  compute_reduced_variables(xflag);

  // This section is heavily edited. Most atom tests have been moved to the convert functions

  for (int i = 0; i < modelorder; i++) // dimension is modelorder
    if (mask[0] & groupbit) { // only works with mask[0]
      y_dot[i] += dtf * y_dot_dot[i]; // changed from dtfm
    }

  update_physical_variables(xflag);

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
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom-> type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (j = 0; j < modelorder; j++){
    if (xflag) y[j] = 0;
    y_dot[j] = 0;
    y_dot_dot[j] = 0.0;

    if (rmass) {
      for (i = 0; i < nlocal; i++){
        iatom = tag[i] - 1;

        if (xflag) {
          y[j] += phi[iatom][j]                * (x[i][0] - x0[iatom][0]);
          y[j] += phi[iatom + nlocal][j]       * (x[i][1] - x0[iatom][1]);
          y[j] += phi[iatom + nlocal*2][j]     * (x[i][2] - x0[iatom][2]);
        }

        y_dot[j] += phi[iatom][j]              * v[i][0];
        y_dot[j] += phi[iatom + nlocal][j]     * v[i][1];
        y_dot[j] += phi[iatom + nlocal*2][j]   * v[i][2];

        y_dot_dot[j] += phi[iatom][j]          * f[i][0] / rmass[i];
        y_dot_dot[j] += phi[iatom+nlocal][j]   * f[i][1] / rmass[i];
        y_dot_dot[j] += phi[iatom+nlocal*2][j] * f[i][2] / rmass[i];                 
      }
    }
    else {
      for (i = 0; i < nlocal; i++){
        iatom = tag[i] - 1;

        if (xflag) {
          y[j] += phi[iatom][j]                * (x[i][0] - x0[iatom][0]);
          y[j] += phi[iatom + nlocal][j]       * (x[i][1] - x0[iatom][1]);
          y[j] += phi[iatom + nlocal*2][j]     * (x[i][2] - x0[iatom][2]);
        }

        y_dot[j] += phi[iatom][j]              * v[i][0];
        y_dot[j] += phi[iatom + nlocal][j]     * v[i][1];
        y_dot[j] += phi[iatom + nlocal*2][j]   * v[i][2];

        y_dot_dot[j] += phi[iatom][j]          * f[i][0] / mass[type[i]];
        y_dot_dot[j] += phi[iatom+nlocal][j]   * f[i][1] / mass[type[i]];
        y_dot_dot[j] += phi[iatom+nlocal*2][j] * f[i][2] / mass[type[i]];                 
      }
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
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++){
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