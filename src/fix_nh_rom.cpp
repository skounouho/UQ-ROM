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
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "respa.h"
#include "update.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNHROM::FixNHROM(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix nvt/nph/npt rom", error);
      dynamic_group_allow = 1;
      time_integrate = 1;

      modelorder = utils::inumeric(FLERR,arg[iarg + 1],false,lmp);
      
      // set up MPI

      MPI_Comm_rank(world,&me);
      MPI_Comm_size(world,&nprocs);

      // create arrays

      int nlocal = atom->nlocal;

      phi = nullptr;
      y = nullptr;
      y_dot = nullptr;
      y_dot_dot = nullptr;  
      x0 = nullptr;

      memory->create(phi, nlocal * 3, modelorder, "FixNHROM:phi");
      memory->create(y, modelorder, "FixNHROM:y");
      memory->create(y_dot, modelorder, "FixNHROM:y_dot");
      memory->create(y_dot_dot, modelorder, "FixNHROM:y_dot_dot");
      memory->create(x0, nlocal, 3, "FixNHROM:x0");

      // MPI variables
      memory->create(y_all, modelorder, "FixNHROM:y_all");
      memory->create(y_dot_all, modelorder, "FixNHROM:y_dot_all");

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
      
      read_rob(arg[iarg + 2], phi);

      if (phi[nlocal * 3 - 1][modelorder - 1] != phi[nlocal * 3 - 1][modelorder - 1]) {
        error->all(FLERR, "Reduced order basis is not filled");
      }

      iarg += 3;
    } else iarg++;
  }
}

/* ----------------------------------------------------------------------
   perform half-step update of velocities in ROM
-----------------------------------------------------------------------*/

void FixNHROM::nve_v()
{
  int xflag = 1;

  compute_reduced_variables(xflag);

  // This section is heavily edited. Most atom group tests have been moved to the convert functions

  for (int i = 0; i < modelorder; i++) {
    y_dot[i] += dtf * y_dot_dot[i]; // changed from dtfm
  }

  MPI_Allreduce(y_dot, y_dot_all, modelorder, MPI_DOUBLE, MPI_SUM, world);
  y_dot = y_dot_all;

  update_physical_variables(xflag);
}

/* ----------------------------------------------------------------------
   perform full-step update of position with ROM
-----------------------------------------------------------------------*/

void FixNHROM::nve_x()
{
  int xflag = 1;

  compute_reduced_variables(xflag);

  // This section is heavily edited. Most atom group tests have been moved to the convert functions

  for (int i = 0; i < modelorder; i++) { // dimension is modelorder
    y[i] += dtv * y_dot[i];
  }

  MPI_Allreduce(y, y_all, modelorder, MPI_DOUBLE, MPI_SUM, world);
  y = y_all;

  update_physical_variables(xflag);
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

  for (j = 0; j < modelorder; j++){
    if (xflag) y[j] = 0;
    y_dot[j] = 0;
    y_dot_dot[j] = 0.0;

    if (rmass) {
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

          y_dot_dot[j] += phi[iatom][j]          * f[i][0] / rmass[i];
          y_dot_dot[j] += phi[iatom+nlocal][j]   * f[i][1] / rmass[i];
          y_dot_dot[j] += phi[iatom+nlocal*2][j] * f[i][2] / rmass[i];                 
        }
    }
    else {
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

void FixNHROM::update_physical_variables(int xflag)
{
  int i,j,iatom;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

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
