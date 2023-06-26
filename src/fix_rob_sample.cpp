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

#include "fix_rob_sample.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "memory.h"
#include "respa.h"
#include "update.h"
#include <fstream>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Cholesky>

using namespace Eigen;
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixROBSample::FixROBSample(LAMMPS *lmp, int narg, char **arg) :
  FixROB(lmp, narg, arg)
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
    // save initial atom positions
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
  BDCSVD<MatrixXd> svd = compute_svd();
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

Eigen::BDCSVD<MatrixXd> FixROB::compute_svd()
{
  int nlocal = atom->nlocal;
  utils::logmesg(lmp, "\nCollected {} snapshots\n", nsnapshots);

  // convert data to Eigen matrix

  MatrixXd snapshotmatrix = convert_to_matrix(snapshots, nsnapshots, nlocal * 3);

  // tranpose and compute singular value decomposition

  snapshotmatrix.transposeInPlace();
  Eigen::BDCSVD<Eigen::MatrixXd> svd(snapshotmatrix, Eigen::ComputeThinU);

  return svd;
}

/* ---------------------------------------------------------------------- */

void FixROB::write_phi(Eigen::BDCSVD<Eigen::MatrixXd> svd) {
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

void write_rob(std::string robfilename, Eigen::MatrixXd U)
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "eiquadprog.hpp"
#include <random>
#include <string>
 
using namespace Eigen;

// CONSTANTS

const int nlocal = 8141;
const int modelorder = 25;
const int nmodels = 3;
const int nsamples = 20;

inline Eigen::MatrixXd ConvertPhiToEigenMatrix(double data[nlocal * 3][modelorder])
{
  int rows = nlocal * 3;
  int columns = modelorder;
  Eigen::MatrixXd eMatrix(rows, columns);
  for (int i = 0; i < rows; ++i)
    eMatrix.row(i) = Eigen::VectorXd::Map(&data[i][0], columns);
  return eMatrix;
}

/* ------------------------------------------------------------------
    Input arguments:
      U0, U1 : points on St(N,n)
         tau : convergence threshold
      
    Output arguments:
       Delta : Log^{St}_U0(U1),
               i.e., tangent vector such that Exp^St_U0(Delta) = U1
   ------------------------------------------------------------------ */
inline Eigen::MatrixXd stiefel_log(Eigen::MatrixXd u0, Eigen::MatrixXd u1, double tau)
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
inline Eigen::MatrixXd stiefel_exp(Eigen::MatrixXd u0, Eigen::MatrixXd delta)
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

inline void generate_samples(Eigen::MatrixXd rob[nmodels + 1], Eigen::MatrixXd u0, MatrixXd stiefel_samples[nsamples])
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

  std::cout << "Solving quadratic programming problem\n";
  solve_quadprog(H, f, Aeq, -beq, CI, ci0, beta);

  // issue with matching solution
  std::cout << "beta = " << beta.transpose() << "\n";
  
  // Dirichlet sampling
  MatrixXd w(nsamples, nmodels);

  // start generator
  std::mt19937 gen(41);

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
  MatrixXd tangential_samples = w * X;

  // project back to Stiefel manifold
  MatrixXd delta[nsamples];

  for (isample = 0; isample < nsamples; isample++) {
    delta[isample] = tangential_samples.row(isample).reshaped(nlocal * 3, modelorder);
    stiefel_samples[isample] = stiefel_exp(u0, delta[isample]);
  }
}

/* ----------------------------------------------------------------------------
    Quadractic Programming solver
   ---------------------------------------------------------------------------- */

/*
 FILE eiquadprog.hh
 
 NOTE: this is a modified of uQuadProg++ package, working with Eigen data structures. 
       uQuadProg++ is itself a port made by Angelo Furfaro of QuadProg++ originally developed by 
       Luca Di Gaspero, working with ublas data structures. 

 The quadprog_solve() function implements the algorithm of Goldfarb and Idnani 
 for the solution of a (convex) Quadratic Programming problem
by means of a dual method.
	 
The problem is in the form:

min 0.5 * x G x + g0 x
s.t.
    CE^T x + ce0 = 0
    CI^T x + ci0 >= 0
	 
 The matrix and vectors dimensions are as follows:
     G: n * n
		g0: n
				
		CE: n * p
	 ce0: p
				
	  CI: n * m
   ci0: m

     x: n
 
 The function will return the cost of the solution written in the x vector or
 std::numeric_limits::infinity() if the problem is infeasible. In the latter case
 the value of the x vector is not correct.
 
 References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
             strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

 Notes:
  1. pay attention in setting up the vectors ce0 and ci0. 
	   If the constraints of your problem are specified in the form 
	   A^T x = b and C^T x >= d, then you should set ce0 = -b and ci0 = -d.
  2. The matrix G is modified within the function since it is used to compute
     the G = L^T L cholesky factorization for further computations inside the function. 
     If you need the original matrix G you should make a copy of it and pass the copy
     to the function.
    
 
 The author will be grateful if the researchers using this software will
 acknowledge the contribution of this modified function and of Di Gaspero's
 original version in their research papers.


LICENSE

Copyright (2011) Benjamin Stephens
Copyright (2010) Gael Guennebaud
Copyright (2008) Angelo Furfaro
Copyright (2006) Luca Di Gaspero


This file is a porting of QuadProg++ routine, originally developed
by Luca Di Gaspero, exploiting uBlas data structures for vectors and
matrices instead of native C++ array.

uquadprog is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

uquadprog is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with uquadprog; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

namespace Eigen {

// namespace internal {

template<typename Scalar>
inline Scalar distance(Scalar a, Scalar b)
{
	Scalar a1, b1, t;
	a1 = std::abs(a);
	b1 = std::abs(b);
	if (a1 > b1) 
	{
		t = (b1 / a1);
		return a1 * std::sqrt(1.0 + t * t);
	}
	else
		if (b1 > a1)
		{
			t = (a1 / b1);
			return b1 * std::sqrt(1.0 + t * t);
		}
	return a1 * std::sqrt(2.0);
}

// }

inline void compute_d(VectorXd &d, const MatrixXd& J, const VectorXd& np)
{
  d = J.adjoint() * np;
}

inline void update_z(VectorXd& z, const MatrixXd& J, const VectorXd& d,  int iq)
{
  z = J.rightCols(z.size()-iq) * d.tail(d.size()-iq);
}

inline void update_r(const MatrixXd& R, VectorXd& r, const VectorXd& d, int iq) 
{
  r.head(iq)= R.topLeftCorner(iq,iq).triangularView<Upper>().solve(d.head(iq));
}

bool add_constraint(MatrixXd& R, MatrixXd& J, VectorXd& d, int& iq, double& R_norm);
void delete_constraint(MatrixXd& R, MatrixXd& J, VectorXi& A, VectorXd& u,  int p, int& iq, int l);

/* solve_quadprog2 is used when the Cholesky decomposition of the G matrix is precomputed */
double solve_quadprog2(LLT<MatrixXd,Lower> &chol,  double c1, VectorXd & g0,  
                      const MatrixXd & CE, const VectorXd & ce0,  
                      const MatrixXd & CI, const VectorXd & ci0, 
                      VectorXd& x);

/* solve_quadprog is used for on-demand QP solving */
inline double solve_quadprog(MatrixXd & G,  VectorXd & g0,  
                      const MatrixXd & CE, const VectorXd & ce0,  
                      const MatrixXd & CI, const VectorXd & ci0, 
                      VectorXd& x){
						  
  LLT<MatrixXd,Lower> chol(G.cols());
  double c1;

  /* compute the trace of the original matrix G */
  c1 = G.trace();

  /* decompose the matrix G in the form LL^T */
  chol.compute(G);

  return solve_quadprog2(chol, c1, g0, CE, ce0, CI, ci0, x);

}

/* solve_quadprog2 is used for when the Cholesky decomposition of G is pre-computed */
inline double solve_quadprog2(LLT<MatrixXd,Lower> &chol,  double c1, VectorXd & g0,  
                      const MatrixXd & CE, const VectorXd & ce0,  
                      const MatrixXd & CI, const VectorXd & ci0, 
                      VectorXd& x)
{
  int i, j, k, l; /* indices */
  int ip, me, mi;
  int n=g0.size();   
  int p=CE.cols(); 
  int m=CI.cols();
  MatrixXd R(g0.size(),g0.size()), J(g0.size(),g0.size());
  
 
  VectorXd s(m+p), z(n), r(m + p), d(n),  np(n), u(m + p);
  VectorXd x_old(n), u_old(m + p);
  double f_value, psi, c2, sum, ss, R_norm;
  const double inf = std::numeric_limits<double>::infinity();
  double t, t1, t2; /* t is the step length, which is the minimum of the partial step length t1 
    * and the full step length t2 */
  VectorXi A(m + p), A_old(m + p), iai(m + p), iaexcl(m+p);
  int q;
  int iq, iter = 0;
 	
  me = p; /* number of equality constraints */
  mi = m; /* number of inequality constraints */
  q = 0;  /* size of the active set A (containing the indices of the active constraints) */
  
  /*
   * Preprocessing phase
   */
	
	
 
  /* initialize the matrix R */
  d.setZero();
  R.setZero();
	R_norm = 1.0; /* this variable will hold the norm of the matrix R */
  
	/* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
  // J = L^-T
  J.setIdentity();
  J = chol.matrixU().solve(J);
	c2 = J.trace();
#ifdef TRACE_SOLVER
 print_matrix("J", J, n);
#endif
  
	/* c1 * c2 is an estimate for cond(G) */
  
	/* 
   * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x 
   * this is a feasible point in the dual space
	 * x = G^-1 * g0
   */
  x = chol.solve(g0);
  x = -x;
	/* and compute the current solution value */ 
	f_value = 0.5 * g0.dot(x);
#ifdef TRACE_SOLVER
  std::cerr << "Unconstrained solution: " << f_value << std::endl;
  print_vector("x", x, n);
#endif
  
	/* Add equality constraints to the working set A */
  iq = 0;
	for (i = 0; i < me; i++)
	{
    np = CE.col(i);
    compute_d(d, J, np);
		update_z(z, J, d,  iq);
		update_r(R, r, d,  iq);
#ifdef TRACE_SOLVER
		print_matrix("R", R, iq);
		print_vector("z", z, n);
		print_vector("r", r, iq);
		print_vector("d", d, n);
#endif
    
    /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint 
      becomes feasible */
    t2 = 0.0;
	if (std::abs(z.dot(z)) > std::numeric_limits<double>::epsilon()) // i.e. z != 0
      t2 = (-np.dot(x) - ce0(i)) / z.dot(np);
    
    x += t2 * z;

    /* set u = u+ */
    u(iq) = t2;
    u.head(iq) -= t2 * r.head(iq);
    
    /* compute the new solution value */
    f_value += 0.5 * (t2 * t2) * z.dot(np);
    A(i) = -i - 1;
    
    if (!add_constraint(R, J, d, iq, R_norm))
    {
      // FIXME: it should raise an error
      // Equality constraints are linearly dependent
      return f_value;
    }
  }
  
	/* set iai = K \ A */
	for (i = 0; i < mi; i++)
		iai(i) = i;
  
l1:	iter++;
#ifdef TRACE_SOLVER
  print_vector("x", x, n);
#endif
  /* step 1: choose a violated constraint */
	for (i = me; i < iq; i++)
	{
	  ip = A(i);
		iai(ip) = -1;
	}
	
	/* compute s(x) = ci^T * x + ci0 for all elements of K \ A */
	ss = 0.0;
	psi = 0.0; /* this value will contain the sum of all infeasibilities */
	ip = 0; /* ip will be the index of the chosen violated constraint */
	for (i = 0; i < mi; i++)
	{
		iaexcl(i) = 1;
		sum = CI.col(i).dot(x) + ci0(i);
		s(i) = sum;
		psi += std::min(0.0, sum);
	}
#ifdef TRACE_SOLVER
  print_vector("s", s, mi);
#endif

    
	if (std::abs(psi) <= mi * std::numeric_limits<double>::epsilon() * c1 * c2* 100.0)
	{
    /* numerically there are not infeasibilities anymore */
    q = iq;
		return f_value;
  }
    
  /* save old values for u, x and A */
   u_old.head(iq) = u.head(iq);
   A_old.head(iq) = A.head(iq);
   x_old = x;
    
l2: /* Step 2: check for feasibility and determine a new S-pair */
	for (i = 0; i < mi; i++)
	{
		if (s(i) < ss && iai(i) != -1 && iaexcl(i))
		{
			ss = s(i);
			ip = i;
		}
	}
  if (ss >= 0.0)
  {
    q = iq;
    return f_value;
  }
    
  /* set np = n(ip) */
  np = CI.col(ip);
  /* set u = (u 0)^T */
  u(iq) = 0.0;
  /* add ip to the active set A */
  A(iq) = ip;

#ifdef TRACE_SOLVER
	std::cerr << "Trying with constraint " << ip << std::endl;
	print_vector("np", np, n);
#endif
    
l2a:/* Step 2a: determine step direction */
  /* compute z = H np: the step direction in the primal space (through J, see the paper) */
  compute_d(d, J, np);
  update_z(z, J, d, iq);
  /* compute N* np (if q > 0): the negative of the step direction in the dual space */
  update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
  std::cerr << "Step direction z" << std::endl;
		print_vector("z", z, n);
		print_vector("r", r, iq + 1);
    print_vector("u", u, iq + 1);
    print_vector("d", d, n);
    print_ivector("A", A, iq + 1);
#endif
    
  /* Step 2b: compute step length */
  l = 0;
  /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
  t1 = inf; /* +inf */
  /* find the index l s.t. it reaches the minimum of u+(x) / r */
  for (k = me; k < iq; k++)
  {
    double tmp;
    if (r(k) > 0.0 && ((tmp = u(k) / r(k)) < t1) )
    {
      t1 = tmp;
      l = A(k);
    }
  }
  /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
  if (std::abs(z.dot(z))  > std::numeric_limits<double>::epsilon()) // i.e. z != 0
    t2 = -s(ip) / z.dot(np);
  else
    t2 = inf; /* +inf */

  /* the step is chosen as the minimum of t1 and t2 */
  t = std::min(t1, t2);
#ifdef TRACE_SOLVER
  std::cerr << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") ";
#endif
  
  /* Step 2c: determine new S-pair and take step: */
  
  /* case (i): no step in primal or dual space */
  if (t >= inf)
  {
    /* QPP is infeasible */
    // FIXME: unbounded to raise
    q = iq;
    return inf;
  }
  /* case (ii): step in dual space */
  if (t2 >= inf)
  {
    /* set u = u +  t * [-r 1) and drop constraint l from the active set A */
    u.head(iq) -= t * r.head(iq);
    u(iq) += t;
    iai(l) = l;
    delete_constraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
    std::cerr << " in dual space: " 
      << f_value << std::endl;
    print_vector("x", x, n);
    print_vector("z", z, n);
		print_ivector("A", A, iq + 1);
#endif
    goto l2a;
  }
  
  /* case (iii): step in primal and dual space */
  
  x += t * z;
  /* update the solution value */
  f_value += t * z.dot(np) * (0.5 * t + u(iq));
  
  u.head(iq) -= t * r.head(iq);
  u(iq) += t;
#ifdef TRACE_SOLVER
  std::cerr << " in both spaces: " 
    << f_value << std::endl;
	print_vector("x", x, n);
	print_vector("u", u, iq + 1);
	print_vector("r", r, iq + 1);
	print_ivector("A", A, iq + 1);
#endif
  
  if (t == t2)
  {
#ifdef TRACE_SOLVER
    std::cerr << "Full step has taken " << t << std::endl;
    print_vector("x", x, n);
#endif
    /* full step has taken */
    /* add constraint ip to the active set*/
		if (!add_constraint(R, J, d, iq, R_norm))
		{
			iaexcl(ip) = 0;
			delete_constraint(R, J, A, u, p, iq, ip);
#ifdef TRACE_SOLVER
      print_matrix("R", R, n);
      print_ivector("A", A, iq);
#endif
			for (i = 0; i < m; i++)
				iai(i) = i;
			for (i = 0; i < iq; i++)
			{
				A(i) = A_old(i);
				iai(A(i)) = -1;
				u(i) = u_old(i);
			}
			x = x_old;
      goto l2; /* go to step 2 */
		}    
    else
      iai(ip) = -1;
#ifdef TRACE_SOLVER
    print_matrix("R", R, n);
    print_ivector("A", A, iq);
#endif
    goto l1;
  }
  
  /* a patial step has taken */
#ifdef TRACE_SOLVER
  std::cerr << "Partial step has taken " << t << std::endl;
  print_vector("x", x, n);
#endif
  /* drop constraint l */
	iai(l) = l;
	delete_constraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
  print_matrix("R", R, n);
  print_ivector("A", A, iq);
#endif
  
  s(ip) = CI.col(ip).dot(x) + ci0(ip);

#ifdef TRACE_SOLVER
  print_vector("s", s, mi);
#endif
  goto l2a;
}


inline bool add_constraint(MatrixXd& R, MatrixXd& J, VectorXd& d, int& iq, double& R_norm)
{
 int n=J.rows();
#ifdef TRACE_SOLVER
  std::cerr << "Add constraint " << iq << '/';
#endif
	int i, j, k;
	double cc, ss, h, t1, t2, xny;
	
  /* we have to find the Givens rotation which will reduce the element
		d(j) to zero.
		if it is already zero we don't have to do anything, except of
		decreasing j */  
	for (j = n - 1; j >= iq + 1; j--)
	{
    /* The Givens rotation is done with the matrix (cc cs, cs -cc).
			 If cc is one, then element (j) of d is zero compared with element
			 (j - 1). Hence we don't have to do anything. 
			 If cc is zero, then we just have to switch column (j) and column (j - 1) 
			 of J. Since we only switch columns in J, we have to be careful how we
			 update d depending on the sign of gs.
			 Otherwise we have to apply the Givens rotation to these columns.
			 The i - 1 element of d has to be updated to h. */
		cc = d(j - 1);
		ss = d(j);
		h = distance(cc, ss);
		if (h == 0.0)
			continue;
		d(j) = 0.0;
		ss = ss / h;
		cc = cc / h;
		if (cc < 0.0)
		{
			cc = -cc;
			ss = -ss;
			d(j - 1) = -h;
		}
		else
			d(j - 1) = h;
		xny = ss / (1.0 + cc);
		for (k = 0; k < n; k++)
		{
			t1 = J(k,j - 1);
			t2 = J(k,j);
			J(k,j - 1) = t1 * cc + t2 * ss;
			J(k,j) = xny * (t1 + J(k,j - 1)) - t2;
		}
	}
  /* update the number of constraints added*/
	iq++;
  /* To update R we have to put the iq components of the d vector
    into column iq - 1 of R
    */
  R.col(iq-1).head(iq) = d.head(iq);
#ifdef TRACE_SOLVER
  std::cerr << iq << std::endl;
#endif
  
	if (std::abs(d(iq - 1)) <= std::numeric_limits<double>::epsilon() * R_norm)
		// problem degenerate
		return false;
	R_norm = std::max<double>(R_norm, std::abs(d(iq - 1)));
	return true;
}


inline void delete_constraint(MatrixXd& R, MatrixXd& J, VectorXi& A, VectorXd& u,  int p, int& iq, int l)
{

  int n = R.rows();
#ifdef TRACE_SOLVER
  std::cerr << "Delete constraint " << l << ' ' << iq;
#endif
	int i, j, k, qq;
	double cc, ss, h, xny, t1, t2;
  
	/* Find the index qq for active constraint l to be removed */
  for (i = p; i < iq; i++)
  if (A(i) == l)
  {
    qq = i;
    break;
  }
      
  /* remove the constraint from the active set and the duals */
  for (i = qq; i < iq - 1; i++)
  {
    A(i) = A(i + 1);
    u(i) = u(i + 1);
    R.col(i) = R.col(i+1);
  }
      
  A(iq - 1) = A(iq);
  u(iq - 1) = u(iq);
  A(iq) = 0; 
  u(iq) = 0.0;
  for (j = 0; j < iq; j++)
    R(j,iq - 1) = 0.0;
  /* constraint has been fully removed */
  iq--;
#ifdef TRACE_SOLVER
  std::cerr << '/' << iq << std::endl;
#endif 
  
  if (iq == 0)
    return;
  
  for (j = qq; j < iq; j++)
  {
    cc = R(j,j);
    ss = R(j + 1,j);
    h = distance(cc, ss);
    if (h == 0.0)
      continue;
    cc = cc / h;
    ss = ss / h;
    R(j + 1,j) = 0.0;
    if (cc < 0.0)
    {
      R(j,j) = -h;
      cc = -cc;
      ss = -ss;
    }
    else
      R(j,j) = h;
    
    xny = ss / (1.0 + cc);
    for (k = j + 1; k < iq; k++)
    {
      t1 = R(j,k);
      t2 = R(j + 1,k);
      R(j,k) = t1 * cc + t2 * ss;
      R(j + 1,k) = xny * (t1 + R(j,k)) - t2;
    }
    for (k = 0; k < n; k++)
    {
      t1 = J(k,j);
      t2 = J(k,j + 1);
      J(k,j) = t1 * cc + t2 * ss;
      J(k,j + 1) = xny * (J(k,j) + t1) - t2;
    }
  }
}

}

#endif
