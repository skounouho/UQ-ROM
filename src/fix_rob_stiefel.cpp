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

#include "fix_rob_stiefel.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "memory.h"
#include "respa.h"
#include "update.h"

#include <complex>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <random>
#include <sstream>
#include <string>
#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>

using namespace Eigen;
using namespace LAMMPS_NS;
using namespace FixConst;

// fix id group rob/stiefel nsamples modelorder file1 file2 file3 sampleformat
// fix 1 all rob/stiefel 20 25 phi_lj.txt phi_lj2.txt phi_lj3.txt phi_global.txt sample/phi_sample{}.txt

/* ---------------------------------------------------------------------- */

FixROBStiefel::FixROBStiefel(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix rob/stiefel", error);

  nlocal = atom->nlocal;

  nsamples = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nsamples <= 0) error->all(FLERR,"Illegal fix rob/stiefel nsamples value: {}", nsamples);
  modelorder = utils::inumeric(FLERR,arg[4],false,lmp);
  if (modelorder <= 0) error->all(FLERR,"Illegal fix rob/stiefel modelorder value: {}", modelorder);
  
  store_files(narg - 6, &arg[5]);
  if (nfile < 2) error->all(FLERR,"No global reduced order basis file specified for fix rob/stiefel");
  nmodels = nfile - 1;

  sampleformat = utils::strdup(arg[narg - 1]);

  // generate and print samples

  // create arrays

  MatrixXd *rob = new MatrixXd[nfile];
  MatrixXd *stiefel_samples = new MatrixXd[nsamples];

  // read rob files into Eigen matrices
  for (int i = 0; i < nfile; i++) {
    rob[i] = read_rob(files[i]);
  }

  // generate samples
  MatrixXd u0 = rob[nmodels];
  generate_samples(rob, u0, stiefel_samples);

  // save all samples to text files
  utils::logmesg(lmp, "Printing bases\n");
  for (int isample = 0; isample < nsamples; isample++) {
    std::string filename = fmt::format(sampleformat, isample + 1);
    write_rob(filename, stiefel_samples[isample]);
  }
}

/* ---------------------------------------------------------------------- */

int FixROBStiefel::setmask()
{
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixROBStiefel::store_files(int nstr, char **str)
{
  nfile = nstr;
  files = new char*[nfile];
  
  // loop through files, no wildcards

  for (int i = 0; i < nfile; i++) {
    files[i] = utils::strdup(str[i]);
  }
}

/* ---------------------------------------------------------------------- 
    Reads the reduced-order basis from a designated file.
   ---------------------------------------------------------------------- */

Eigen::MatrixXd FixROBStiefel::read_rob(std::string robfile)
{
  utils::logmesg(lmp, "Reading reduced-order basis file {}\n", robfile);
  Eigen::MatrixXd phi(nlocal * 3, modelorder);
  phi.setZero();

  try {
    std::ifstream file(robfile);
    std::string line;
    if (file.is_open()) {
      for (int i = 0; i < nlocal * 3; i++) {
        std::getline(file, line);
        std::stringstream ss(line);
        for (int j = 0; j < modelorder; j++) {
          ss >> phi(i,j);
        }
      }
    }
  } catch (std::exception &e) {
    error->one(FLERR, "Error reading reduced-order basis: {}", e.what());
  }

  return phi;
}

/* ---------------------------------------------------------------------- */

void FixROBStiefel::write_rob(std::string robfilename, Eigen::MatrixXd U)
{
  int i, j;

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
    }
  } catch (std::exception &e) {
    std::cout << "Error writing reduced-order basis";
  }
}

/* ------------------------------------------------------------------
    Input arguments:
      U0, U1 : points on St(N,n)
         tau : convergence threshold
      
    Output arguments:
       Delta : Log^{St}_U0(U1),
               i.e., tangent vector such that Exp^St_U0(Delta) = U1
   ------------------------------------------------------------------ */
inline Eigen::MatrixXd FixROBStiefel::stiefel_log(Eigen::MatrixXd u0, Eigen::MatrixXd u1, double tau)
{
  // get dimensions
  int p = u0.cols();

  // step 1
  MatrixXd M = u0.transpose()*u1;

  // step 2, thin qr of normal component of u1
  Eigen::HouseholderQR<MatrixXd> qr1(u1 - u0 * M); // Householder QR is not the most stable numerically
  MatrixXd Q = qr1.householderQ();
  Q.conservativeResize(NoChange, M.cols()); // make Q thin

  MatrixXd N = qr1.matrixQR();
  N.conservativeResize(M.rows(), NoChange); // make N thin
  N = N.triangularView<Upper>(); // make triangular

  // step 3, orthogonal completion
  MatrixXd MN(M.rows() + N.rows(), M.cols()); // concatenate M and N vertically
  MN << M, N;
  Eigen::HouseholderQR<MatrixXcd> qr2(MN);
  MatrixXcd V = qr2.householderQ();

  // "Procrustes preprocessing"
  BDCSVD<MatrixXcd> svd(V(seq(p, p * 2 - 1), seq(p, p * 2 - 1)), ComputeFullU | ComputeFullV);
  MatrixXcd D = svd.matrixU();
  MatrixXcd R = svd.matrixV();
  V(all, seq(p, p * 2 - 1)) = V(all, seq(p, p * 2 - 1)) * (R * D.transpose());
  V << MN, V(all, seq(p, p * 2 - 1));

  // step 4, for-loop

  MatrixXcd LV;

  for (int k=0; k < 10000; k++) {

    // step 5
    LV = V.log();
    MatrixXcd C = LV(seq(p, p * 2 - 1), seq(p, p * 2 - 1));

    // steps 6 - 8: convergence check

    BDCSVD<MatrixXcd> svdC(C);
    double normC = svdC.singularValues()(0);

    if (normC < tau) {
      break;
    }

    // step 9
    MatrixXcd phi = (-C).exp();
    V(all, seq(p, p * 2 - 1)) = V(all, seq(p, p * 2 - 1)) * phi;
  }

  // prepare output
  // upon convergence, we have logm(V) =
  //    A = LV(0:p-1, 0:p-1);   B = LV(p+1:2*p, 1:p)
  // delta = U0 * A + Q * B

  MatrixXcd delta = u0 * LV(seq(0,p - 1), seq(0,p - 1)) + Q * LV(seq(p, p * 2 - 1), seq(0, p - 1));

  return delta.real();
}

/* Tested the function above with the inputs
    u0 = [1 0
          0 1
          0 0]
     
    u1 = [0.8944 0.7071
          0      0.7071
          0.4472 0 ]
    
    tau = 1e-4
    
    and got the result
    Delta = [ 0.0515345  0.782797
              -0.118107  -0.273923
               0.456774  -0.197914 ]
    
    which matches the output from Hao's stiefel_log funciton. */

/* ------------------------------------------------------------------
    Input arguments:
         U0 : base point on St(N/n)
      delta : tangent vector on TuSt(N,n)
   ------------------------------------------------------------------ */
inline Eigen::MatrixXd FixROBStiefel::stiefel_exp(Eigen::MatrixXd u0, Eigen::MatrixXd delta)
{
  MatrixXd A = u0.transpose() * delta;
  int r = u0.cols();

  // thin qr decomposition
  Eigen::HouseholderQR<MatrixXd> qr(delta - u0 * A);
  MatrixXd Q = qr.householderQ();
  Q.conservativeResize(NoChange, A.cols()); // make Q thin
  MatrixXd R = qr.matrixQR();
  R.conservativeResize(A.rows(), NoChange); // make N thin
  R = R.triangularView<Upper>(); // make triangular

  // set up matrix
  MatrixXd AR(A.rows() + R.cols(), A.cols() + R.cols());
  AR << A, -R.transpose(), R, MatrixXd::Zero(r,r);

  // eigenvalue decomposition
  Eigen::EigenSolver<MatrixXd> eig(AR,true);
  MatrixXcd V = eig.eigenvectors();
  MatrixXcd D = eig.eigenvalues().asDiagonal();

  // MATLAB is not consistent when it creates eigenvectors,
  // so the eigenvectors computed by Eigen may have differen
  // values. However, the result is the same.

  MatrixXcd MN = V * D.exp() * V.adjoint() * MatrixXd::Identity(r * 2,r);
  MatrixXd M = MN(seq(0,r - 1), all).real();
  MatrixXd N = MN(seq(r, last), all).real();

  MatrixXd u = u0 * M + Q * N;

  return u;
}

/* Tested the function above with the inputs
    u0 = [1 0
          0 1
          0 0]
     
    delta = [ 0.0515345  0.782797
              -0.118107  -0.273923
               0.456774  -0.197914 ]
    
    and got the result
    u =  [ 0.796972  -0.0639456
          -0.249794   0.364992
           0.425455   0.175382 ]
    
    which matches the output from Hao's stiefel_exp funciton. */

inline void FixROBStiefel::generate_samples(Eigen::MatrixXd *rob, Eigen::MatrixXd u0, MatrixXd *stiefel_samples)
{
  int i, isample;

  // project onto the tangential space
  // assume that the Log is well-defined on all the points
  MatrixXd Vs[nmodels + 1];
  double tau = 1e-4;
  for (i = 0; i < nmodels; i++) {
    Vs[i] = stiefel_log(u0, rob[i], tau);
  }

  Vs[nmodels] = MatrixXd::Zero(nlocal * 3, modelorder); // project base point to zero

  // solve quadratic programming problem

  MatrixXd X(nlocal * 3 * modelorder, nmodels);
  for (i = 0; i < nmodels; i++) {
    X.col(i) = Vs[i].reshaped().transpose();
  }

  X.transposeInPlace();
  
  // set up matrices and constraints
  MatrixXd H = X * X.transpose();
  VectorXd f = VectorXd::Zero(nmodels);
  MatrixXd Aeq = MatrixXd::Ones(nmodels, 1);
  VectorXd beq = VectorXd::Ones(1);

  // set lower bound to 0
  MatrixXd CI = MatrixXd::Identity(nmodels, nmodels);
  VectorXd ci0 = VectorXd::Zero(nmodels);

  VectorXd beta(nmodels);

  utils::logmesg(lmp, "Solving quadratic programming problem\n");
  solve_quadprog(H, f, Aeq, -beq, CI, ci0, beta);

  // print gamma values
  utils::logmesg(lmp, "gamma = ");
  for (i = 0; i < beta.size(); i++) {
    utils::logmesg(lmp, "  {}", beta(i));
  }
  
  // Dirichlet sampling
  MatrixXd w(nsamples, nmodels);

  // start generator
  std::mt19937 gen(41);

  utils::logmesg(lmp, "\nSampling gamma distribution\n");
  // loop through shape values
  for (i = 0; i < nmodels; i++) {
    double a = beta(i);
    std::gamma_distribution<> d(a, 1);
    for (isample = 0; isample < nsamples; isample++) {
      w(isample, i) = d(gen);
    }
  }

  // ensure rows are a valid probability distribution
  for (isample = 0; isample < nsamples; isample++) {
    w.row(isample) = w.row(isample) / w.row(isample).sum();
  }

  // generate samples
  utils::logmesg(lmp, "Generating samples\n");
  MatrixXd tangential_samples = w * X;

  // project back to Stiefel manifold
  MatrixXd delta[nsamples];

  utils::logmesg(lmp, "Retracting samples from Stiefel manifold\n");
  for (isample = 0; isample < nsamples; isample++) {
    delta[isample] = tangential_samples.row(isample).reshaped(nlocal * 3, modelorder);
    stiefel_samples[isample] = stiefel_exp(u0, delta[isample]);
  }
}


