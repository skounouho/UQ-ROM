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
#include <Eigen/Eigen>

namespace LAMMPS_NS {

class FixROB : public Fix {
 public:
  FixROB(class LAMMPS *, int, char **);

  int setmask() override;
  void end_of_step() override;
  void post_run() override;
  Eigen::MatrixXd convert_to_matrix(double **, int, int);
  Eigen::BDCSVD<Eigen::MatrixXd> compute_svd();
  void write_phi(Eigen::BDCSVD<Eigen::MatrixXd>);

 private:
   int modelorder;
   int nsnapshots;
   double **snapshots;
   double **x0;
   char *robfilename;
   int fullorderflag;
   int nmodels;
};

}    // namespace LAMMPS_NS

#endif
#endif