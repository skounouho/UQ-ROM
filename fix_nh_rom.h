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

#ifndef LMP_FIX_NH_ROM_H
#define LMP_FIX_NH_ROM_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNHROM : public FixNH {
 public:
  FixNHROM(class LAMMPS *, int, char **);

  int modelorder;
  double **phi;
  double **x0;         // initial position of particles
  double *y;           // reduced order position
  double *y_dot;       // reduced order velocity
  double *y_dot_dot;   // reduced order acceleration 

 protected:
  double inertia;

//   void nve_v() override;
  void nve_x() override;

  //******************* ADDED ******************

  // helper methods
  
  void read_rob(std::string, double**);
  void compute_reduced_variables(int);
  void update_physical_variables(int);
   
};

}    // namespace LAMMPS_NS

#endif
