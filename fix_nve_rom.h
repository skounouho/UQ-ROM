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

//   int setmask() override;
//   void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
//   void initial_integrate_respa(int, int, int) override;
//   void final_integrate_respa(int, int) override;
//   void reset_dt() override;


  //******************* ADDED ******************

  void read_rob(std::string);
  void read_mean(std::string);
  void remove_periodic_boundary(double, double, double);
  void reinstate_periodic_boundary(double, double, double);
  void convert_physical_to_reduced();
  void convert_reduced_to_physical();

  //**************** END ADDED *****************
 
 protected:
//   double dtv, dtf;
//   double *step_respa;
//   int mass_require;

   // *************************** ADDED ****************************

   int modelorder;
   double **phi;
   double **mean;
   double *A;
   double *V;
   double *X;

   //**************** END ADDED *****************

};

}    // namespace LAMMPS_NS

#endif
#endif