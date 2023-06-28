/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nve/rom,FixNVEROM);
// clang-format on
#else

#ifndef LMP_FIX_NVE_ROM_H
#define LMP_FIX_NVE_ROM_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEROM : public FixNVE {
 public:
  FixNVEROM(class LAMMPS *, int, char **);

  void initial_integrate(int) override;
  void final_integrate() override;
  double compute_vector(int) override;

 protected:
  int me, nprocs;
  int modelorder;      // order of reduced order model
  double **phi;        // reduced-order basis
  double **x0;         // initial position of particles
  double *y;           // reduced order position
  double *y_dot;       // reduced order velocity
  double *y_dot_dot;   // reduced order acceleration

  // implementing MPI run capability
  double *y_all;       // reduced order position across processors
  double *y_dot_all;   // reduced order velocity across processors

  // helper methods

  void read_rob(std::string, double**);
  void compute_reduced_variables(int);
  void update_physical_variables(int);
};

}    // namespace LAMMPS_NS

#endif
#endif