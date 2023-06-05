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

#include "fix_rob.h"

#include "atom.h"
#include "force.h"
#include "respa.h"
#include "update.h"

// ****************** ADDED *********************

#include "memory.h"
#include "error.h"
#include <fstream>
#include <Eigen/Dense>

#include "iostream" // for debugging

using namespace Eigen;

// *************** END ADDED *******************

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixROB::FixROB(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix rob", error);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0)
    error->all(FLERR,"Illegal fix rob nevery value: {}", nevery);
  
  robfilename = utils::strdup(arg[5]);

  int nlocal = atom->nlocal;

  if (strcmp(arg[4], "full") == 0)  modelorder = nlocal * 3;
    else modelorder = utils::inumeric(FLERR,arg[4],false,lmp);

  phi = nullptr;
  snapshots = nullptr;

  memory->create(phi, nlocal * 3, modelorder, "FixROB:phi");
  memory->create(snapshots, 1, nlocal * 3, "FixROB:snapshots");
  memory->create(initial, nlocal, 3, "FixROB:initial");

  double **x = atom->x;
  int *tag = atom->tag;
  int atomid;

  for (int i = 0; i < nlocal; i++) {
    atomid = tag[i] - 1;
    initial[atomid][0] = x[i][0];
    initial[atomid][1] = x[i][1];
    initial[atomid][2] = x[i][2];
  }

  nsnapshots = 0;
}

/* ---------------------------------------------------------------------- */

int FixROB::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixROB::end_of_step()
{

  int atomid;

  double **x = atom->x;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  if (nsnapshots > 0) {
    memory->grow(snapshots, nsnapshots + 1, nlocal * 3, "FixROB:snapshots");
  }

  for (int i = 0; i < nlocal; i++) {
    atomid = tag[i] - 1;
    snapshots[nsnapshots][atomid] = x[i][0];
    snapshots[nsnapshots][atomid + nlocal] = x[i][1];
    snapshots[nsnapshots][atomid + nlocal*2] = x[i][2];
  }
  
  nsnapshots++;
}

/* ---------------------------------------------------------------------- */

void FixROB::post_run()
{
  utils::logmesg(lmp, "Collected {} snapshots\n", nsnapshots);

  // shift positions by initial

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < nsnapshots; j++) {
      snapshots[j][i] = snapshots[j][i] - initial[i][0];
      snapshots[j][i + nlocal] = snapshots[j][i + nlocal] - initial[i][1];
      snapshots[j][i + nlocal*2] = snapshots[j][i + nlocal*2] - initial[i][2];
    }
  }

  // convert data to Eigen matrix

  MatrixXd snapshotmatrix = ConvertToEigenMatrix(snapshots, nsnapshots, nlocal * 3);
  snapshotmatrix.transposeInPlace();
  Eigen::BDCSVD<Eigen::MatrixXd> svd(snapshotmatrix, Eigen::ComputeFullU);
  MatrixXd U = svd.matrixU();

  try {
    std::ofstream file(robfilename);
    if (file.is_open()) {
      for (int i = 0; i < nlocal * 3; i++) {
        file << U(i, 0); // write first values
        for (int j = 1; j < modelorder; j++) {
          file << " " << U(i, j);
        }
        file << '\n';
      }
    }
  } catch (std::exception &e) {
    error->one(FLERR, "Error writing reduced-order basis: {}", e.what());
  }

  // save singular values

  VectorXd S = svd.singularValues();
  try {
    std::ofstream file("lambda.txt");
    if (file.is_open()) {
      for (int i = 0; i < S.size(); i++) {
        file << S(i) << '\n';
      }
    }
  } catch (std::exception &e) {
    error->one(FLERR, "Error writing singular values: {}", e.what());
  }

}

/* ---------------------------------------------------------------------- */

Eigen::MatrixXd FixROB::ConvertToEigenMatrix(double **data, int rows, int columns)
{
    Eigen::MatrixXd eMatrix(rows, columns);
    for (int i = 0; i < rows; ++i)
        eMatrix.row(i) = Eigen::VectorXd::Map(&data[i][0], columns);
    return eMatrix;
}