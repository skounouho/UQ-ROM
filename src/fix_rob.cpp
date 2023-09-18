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
#include "error.h"
#include "force.h"
#include "modify.h"
#include "memory.h"
#include "respa.h"
#include "update.h"
#include <fstream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/SVD>

using namespace Eigen;
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

  int nlocal = atom->nlocal;
  fullorderflag = 0;

  if (strcmp(arg[4], "full") == 0) fullorderflag = 1;
    else modelorder = utils::inumeric(FLERR,arg[4],false,lmp);
  
  robfilename = utils::strdup(arg[5]);

  // create arrays

  snapshots = nullptr;
  x0 = nullptr;
  memory->create(snapshots, 1, nlocal * 3, "FixROB:snapshots");
  memory->create(x0, nlocal, 3, "FixNVEROM:x0");

  // initialize snapshot count

  nsnapshots = 0;

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

  int iatom;
  double **x = atom->x;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  if (update->ntimestep == 0) {

    // save initial atom positions and skip timestep 0
    // this is for rerun simulations

    for (int i = 0; i < nlocal; i++) {
      iatom = tag[i] - 1;
      x0[iatom][0] = x[i][0];
      x0[iatom][1] = x[i][1];
      x0[iatom][2] = x[i][2];
    }
    return;
  }

  if (nsnapshots > 0) {
    memory->grow(snapshots, nsnapshots + 1, nlocal * 3, "FixROB:snapshots");
  }

  // shift by initial position and save

  for (int i = 0; i < nlocal; i++) {
    iatom = tag[i] - 1;
    snapshots[nsnapshots][iatom] = x[i][0] - x0[iatom][0];
    snapshots[nsnapshots][iatom + nlocal] = x[i][1] - x0[iatom][1];
    snapshots[nsnapshots][iatom + nlocal*2] = x[i][2] - x0[iatom][2];
  }

  nsnapshots++;
}

/* ---------------------------------------------------------------------- */

void FixROB::post_run()
{
  if (nsnapshots <= 0) {
    error->warning(FLERR,"Did not collect any snapshots, unable to compute reduced order basis");
    return;
  }
  BDCSVD<MatrixXd, Eigen::ComputeThinU> svd = compute_svd();
  write_phi(svd);
}

/* ---------------------------------------------------------------------- */

Eigen::MatrixXd FixROB::convert_to_matrix(double **data, int rows, int columns)
{
    Eigen::MatrixXd eMatrix(rows, columns);
    for (int i = 0; i < rows; ++i)
        eMatrix.row(i) = Eigen::VectorXd::Map(&data[i][0], columns);
    return eMatrix;
}

/* ---------------------------------------------------------------------- */

Eigen::BDCSVD<MatrixXd, Eigen::ComputeThinU> FixROB::compute_svd()
{
  int nlocal = atom->nlocal;
  utils::logmesg(lmp, "\nCollected {} snapshots\n", nsnapshots);

  // convert data to Eigen matrix

  MatrixXd snapshotmatrix = convert_to_matrix(snapshots, nsnapshots, nlocal * 3);

  // tranpose and compute singular value decomposition

  snapshotmatrix.transposeInPlace();
  Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU> svd(snapshotmatrix);

  return svd;
}

/* ---------------------------------------------------------------------- */

void FixROB::write_phi(Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU> svd) {
  int i,j;
  int nlocal = atom->nlocal;

  // get U matrix of svd [X] = [U][S][V]'

  MatrixXd U = svd.matrixU();

  // get singular values

  VectorXd S = svd.singularValues();
  int nsingularvals = svd.nonzeroSingularValues();
  utils::logmesg(lmp, "The rank is approximately {}\n", nsingularvals);

  // cap model order to number of singular values (rank)
  
  if (fullorderflag) modelorder = nsingularvals;
  else if (modelorder > nsingularvals) modelorder = nsingularvals;

  // print result to file

  try {
    std::ofstream file(robfilename);
    file << std::fixed;
    file << std::setprecision(16);
    if (file.is_open()) {
      for (i = 0; i < nlocal * 3; i++) {
        file << U(i,0); // write first values
        for (j = 1; j < modelorder; j++) {
          file << " " << U(i,j);
        }
        file << '\n';
      }
      file << "\n# singular values\n";
      for (i = 0; i < S.size(); i++) {
        file << S(i) << " ";
      }
    }
  } catch (std::exception &e) {
    error->one(FLERR, "Error writing reduced-order basis: {}", e.what());
  }
}