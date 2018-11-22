/*
################################################################################
EM algorithm for mixture of multidimensional gaussians
developed by Luca Galbusera on September/October 2015.
Based on the tutorial "EM Demystified: An Expectation-Maximization Tutorial"
by Yihua Chen and Maya R. Gupta. UWEE Technical Report Number UWEETR-2010-0002
################################################################################
*/

#ifndef __EM__
#define __EM__

#include <iostream>
#include <cmath>
#include <RcppEigen.h>

using namespace Eigen;
using namespace std;

class C_EM{
public:
  //General variables
  int n, c, d; //Number of data components and dimensions
  int max_iter;
  double tol;
  Matrix<double, Dynamic, Dynamic> y; //Observed data
  double log_lik; //log likelihood used to check convergence

  //Model variables
  Matrix<double, Dynamic, Dynamic> g_ij; //gamma_ij
  Matrix<double, Dynamic, 1> n_j; //sum_i gamma_ij

  //Model parameters
  Matrix<double, Dynamic, Dynamic> mu;     //vector of vector of means
  Matrix<double, Dynamic, 1> w ;           //vector of vector of weights
  Matrix<double, Dynamic, Dynamic> *sigma; //vector of covariance matrices
  double Delta; //Range of uniform distribution

  //Constructor
  C_EM(const int _n, const int _c, const int _d,
       Matrix<double, Dynamic, Dynamic> _y, Matrix<double, Dynamic, Dynamic> _mu,
       Matrix<double, Dynamic, Dynamic> _w, Matrix<double, Dynamic, Dynamic> _sigma[],
       double _tol=0.001, int _max_iter=200, double _Delta=0);

  //Print the current values of the parameters
  void print_values(void);
  //Computes the log-likelihood
  double LL(void);
  
  //Gaussian distribution
  double phi(VectorXd y, VectorXd mu, Matrix<double, Dynamic, Dynamic> sigma);
  
  //Start the em algorithm
  void start_em(void);

  //Return the posterior distribution
  Matrix<double, Dynamic, Dynamic> Posterior();

private:
  void E_step(void);
  void M_step(void);
};

#endif
