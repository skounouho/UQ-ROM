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
FixStyle(rob/sample,FixROBSample);
// clang-format on
#else

#ifndef LMP_FIX_ROB_SAMPLE_H
#define LMP_FIX_ROB_SAMPLE_H

#include "fix_rob.h"
#include <Eigen/Core>
#include <Eigen/SVD>

namespace LAMMPS_NS {

class FixROBSample : public FixROB {
 public:
  FixROBSample(class LAMMPS *, int, char **);

  virtual int setmask() override;
  virtual void end_of_step() override;
  virtual void post_run() override;
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
};

}    // namespace LAMMPS_NS

#endif
#endif