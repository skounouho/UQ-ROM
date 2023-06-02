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

 protected:
  double inertia;

  void nve_v() override;
  void nve_x() override;

  //******************* ADDED ******************

  void read_rob(std::string, double**);
  void read_mean(std::string, double**);
  void convert_physical_to_reduced(double *, double *, double *);
  void convert_reduced_to_physical(double *, double *);

 protected:

   int modelorder;
   double **phi;
   double **mean;
   double *A;
   double *V;
   double *X;
};

}    // namespace LAMMPS_NS

#endif
