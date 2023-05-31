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
  if (narg != 7) utils::missing_cmd_args(FLERR, "fix rob", error);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0)
    error->all(FLERR,"Illegal fix rob nevery value: {}", nevery);
  
  robfilename = utils::strdup(arg[5]);
  meanfilename = utils::strdup(arg[6]);

  int nlocal = atom->nlocal;

  if (strcmp(arg[4], "full") == 0)  modelorder = nlocal * 3;
    else modelorder = utils::inumeric(FLERR,arg[4],false,lmp);

  phi = nullptr;
  means = nullptr;
  snapshots = nullptr;

  memory->create(phi, nlocal * 3, modelorder, "FixROB:phi");
  memory->create(means, nlocal * 3, 1, "FixROB:means");
  memory->create(snapshots, 1, nlocal * 3, "FixROB:snapshots");

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

void FixROB::init()
{
  end_of_step();
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
  // std::cout << "I made " << nsnapshots << " snapshots\n";

  // std::cout << "I made it to the end! Here are the first three values I found:" << '\n';
  // std::cout << snapshots[0][268 + 272] << ", " << snapshots[1000][268 + 272] << ", " << snapshots[2000][268 + 272] << "\n\n";

  // calculate means

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal * 3; i++) {
    means[i][0] = 0;
    for (int j = 0; j < nsnapshots; j++) {
      means[i][0] += snapshots[j][i];
    }
    means[i][0] = means[i][0] / nsnapshots;
  }

  // output mean
  
  try {
    std::ofstream file(meanfilename);
    if (file.is_open()) {
      for (int i = 0; i < nlocal * 3; i++) {
        file << means[i][0] << '\n';
      }
    }
  }
  catch (std::exception &e) {
    error->one(FLERR, "Error writing means: {}", e.what());
  }
  

  // std::cout << "I made past the means! Here are the first three values I found:" << '\n';
  // std::cout << means[0][0] << ", " << means[1][0] << ", " << means[2][0] << "\n\n";

  // subtract means

  for (int i = 0; i < nlocal * 3; i++) {
    for (int j = 0; j < nsnapshots; j++) {
      snapshots[j][i] = snapshots[j][i] - means[i][0];
    }
  }

  // convert data to Eigen matrix

  MatrixXd snapshotmatrix = ConvertToEigenMatrix(snapshots, nsnapshots, nlocal * 3);
  snapshotmatrix.transposeInPlace();
  // snapshotmatrix /= sqrt(nsnapshots - 1); // normalizing
  Eigen::BDCSVD<Eigen::MatrixXd> svd(snapshotmatrix, Eigen::ComputeFullU);
  MatrixXd U = svd.matrixU();

  // Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> pod(snapshotmatrix);
  // MatrixXd U = pod.householderQ();

  // std::cout << svd.singularValues() << '\n';
  // std::cout << svd.rank() << '\n';;

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


}

/* ---------------------------------------------------------------------- */

Eigen::MatrixXd FixROB::ConvertToEigenMatrix(double **data, int rows, int columns)
{
    Eigen::MatrixXd eMatrix(rows, columns);
    for (int i = 0; i < rows; ++i)
        eMatrix.row(i) = Eigen::VectorXd::Map(&data[i][0], columns);
    return eMatrix;
}