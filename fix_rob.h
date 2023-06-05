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
FixStyle(rob,FixROB);
// clang-format on
#else

#ifndef LMP_FIX_ROB_H
#define LMP_FIX_ROB_H

#include "fix.h"

#include "Eigen/Dense"

namespace LAMMPS_NS {

class FixROB : public Fix {
 public:
  FixROB(class LAMMPS *, int, char **);

  int setmask() override;
  void end_of_step() override;
  void post_run() override;
  Eigen::MatrixXd ConvertToEigenMatrix(double **, int, int);
 
 protected:
   int modelorder;
   int nsnapshots;
   double **snapshots;
   double **phi;
   double **initial;
   std::string robfilename;

};

}    // namespace LAMMPS_NS

#endif
#endif