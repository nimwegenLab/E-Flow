/*######################################################################
 Written by Luca Galbusera on November 2015 at Nimwegen's lab, Basel.
 This program performs the EM update for a mixture of c gaussians and
one uniform distribution. The width of the uniform if not specified is
computed from the data, while the starting values.

The parameter Delta (the uniform width) is optional.
The means must be stored in a matrix, one column for each gaussian component
so they will form a matrix dimXcomp

The main functions of the class EM are
- EM.start(): start the iterative EM
- EM.print_values(): print the current parameters
- EM.LL(): print the log-likelihood
- EM.Posterior(): print the posterior probability for each data to belong to a component
######################################################################*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include "em_class.hpp"

using namespace Rcpp;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const int max_iter = 500; //max number of iterations before stopping the algorithm
const double tol = 1e-10; //tollerance (difference in likelihood to stop the em)


//' Fit of a mixture of gaussians and one uniform using EM
//'
//' This program performs the EM update for a mixture of different gaussians and
//' one uniform distribution. The width of the uniform if not specified is
//' computed from the data. The number of gaussians in the mixture is determined
//' by the length of the parameter \code{mu_}.
//'
//' @param y_ The data to be fit with the mixture. It must be a matrix where each
//'   column is a variable and each row is a dimension.
//' @param mu_ Vector giving the initial values for the means of the gaussians.
//'   Its length determines the number of gaussians to be included in the
//'   mixture.
//' @param w_ Vector giving the initial wights of the gaussians and the uniform.
//'   It must sum to 1.
//' @param list_sigma_ List of matrices each giving the starting covariance
//'   matrix for each gaussians.
//' @param delta_ Dimensions of the uniform distribution. If it's 0 then the
//'   dimensions are the difference between the maximum and the minum of the data
//'   in each direction.
//' @param verbose_ Should some information be printed?
//'
//' @return A list containing the fitted statistics of the mixture and a matrix
//'   containing the posterior probability for each point to come from each
//'   component of the mixture (one column is one component, one row is one data
//'   point).
//'
//' @export
// [[Rcpp::export]]
SEXP EM_mixture(SEXP y_, SEXP mu_, SEXP w_, SEXP list_sigma_, SEXP delta_, SEXP verbose_){
  //Map R ojects to Eigen
  Map<MatrixXd> y  (as<Map<MatrixXd> >(y_));
  Map<MatrixXd> mu (as<Map<MatrixXd> >(mu_));  //means of each component by column
  Map<VectorXd> w  (as<Map<VectorXd> >(w_));
  List list_sigma(list_sigma_);
  double delta = as<double>(delta_);
  double verbose = as<bool>(verbose_);

  //Infer parameters from the shape of the matrices
  int n = y.rows();  //number of data
  int d = y.cols();  //number of dimensions
  int c = mu.cols(); //number of gaussians in the mixture (without the uniform)

  //Create the sigma matrices
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *sigma;
  sigma = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>[c];
  for(int i=0; i<list_sigma.size(); i++){
    sigma[i] = list_sigma[i];
    if(sigma[i].rows()!=d && sigma[i].cols()!=d){
      Rcerr << "Dimensions of the variance matrices is not correct" << endl;
      return wrap(0);
    }
  }

  //Check consistency of the parameters
  if(list_sigma.size()!=c){
    Rcerr << "Number of sigma matrices different from number of components" << endl;
    return wrap(0);
  }
  if(d!=mu.rows()){
    Rcerr << "Mean and data dimensions not compatible" << endl;
    return wrap(0);}
  else if(c!=w.size()-1){ //w has one more component for the uniform
    Rcerr << "Number of components in the mean and in the weights not compatible" << endl;
    return wrap(0);
  }

  //Cout some information
  if(verbose){
    Rcerr << "Number of data: " << n << endl;
    Rcerr << "Number of dimensions: " << d << endl;
    Rcerr << "Number of gaussian components: " << c << endl;
  }
  C_EM em(n, c, d, y, mu, w, sigma, tol, max_iter, delta);

  if(verbose) em.print_values();
  em.start_em();

  //Create the elements to be returned
  List cov;
  for(int i=0; i<c; i++)
    cov.push_back(em.sigma[i]);

  return Rcpp::List::create(Rcpp::Named("post") = em.Posterior(),
                            Rcpp::Named("mean") = em.mu,
                            Rcpp::Named("weights") = em.w,
                            Rcpp::Named("covariances") = cov,
                            Rcpp::Named("LL") = em.LL());
}
