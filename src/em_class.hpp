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
       double _tol=0.001, int _max_iter=200, double _Delta=0) :
    n(_n), c(_c), d(_d),
    y(_y), mu(_mu), w(_w), sigma(_sigma),
    tol(_tol), max_iter(_max_iter), Delta(_Delta)
  {
    //Create matrix g_ij and n_j with c+1 components (c gaussians and 1 uniform)
    g_ij.resize(n, c+1);
    n_j.resize(c+1);

    //If Delta is zero, then estimate it from the data
    if(Delta==0){
      double Delta_min, Delta_max;
      Delta = 1;
      for(int i=0; i<d; i++)
	Delta *= y.col(i).maxCoeff()-y.col(i).minCoeff();
    }


    //Check that the variance matrix is invertible
    for(int i=0; i<c; i++){
      if(abs(sigma[i].determinant())<1e-13)
	cerr << "The variance of the " << i << "-th component is not invertible!" << endl;
    }


  }


  void print_values(void){
    cerr << "Weigths w (the last is the uniform):" << endl;
    cerr << w << endl << endl;

    cerr << "Delta: " << Delta << endl;

    cerr << "Means mu (by columns):" << endl;
    cerr << mu << endl << endl;

    cerr << "Covariance matrices sigma: " << endl;
    for(int i=0; i<c; i++)
    cerr << sigma[i] << endl << endl;
  }

  double LL(void){ //Computes the log-likelihood
    double inner_sum=0;
    double outer_sum=0;
    for(int i=0; i<n; i++){
      inner_sum=0;
      for(int j=0; j<c; j++)
	inner_sum += w(j)*phi(y.row(i), mu.col(j), sigma[j]);
      outer_sum += log(inner_sum+w(c)/Delta);
    }

    return outer_sum/n;
  }


  //Gaussian distribution
  double phi(VectorXd y, VectorXd mu,
	     Matrix<double, Dynamic, Dynamic> sigma){
    double norm = pow(2.*M_PI, d/2.)*sqrt(abs(sigma.determinant()));
    double exponent = (0.5*(y-mu).transpose()*sigma.inverse()*(y-mu))(0);
    return exp(-exponent)/norm;
  }

  void start_em(void){
    double ll_old = 1e10;
    double ll_new = LL();
    int num_iter = 0;
    while(abs(ll_new-ll_old)>tol && num_iter++ < max_iter){
      ll_old=ll_new;
      E_step();
      M_step();
      ll_new = LL();
    }
    //print_values();
    return;
  }

  Matrix<double, Dynamic, Dynamic> Posterior(){
    return g_ij;
  }

private:
  void E_step(void){
    //Estimate g_ij
    double norm=0;
    for(int i=0; i<n; i++){
      norm=w(c)/Delta;   //Uniform
      g_ij(i,c) = w(c)/Delta; //Uniform
      for(int j=0; j<c; j++){ //Gaussians
	g_ij(i,j) = w(j)*phi(y.row(i), mu.col(j), sigma[j]);
	norm += w(j)*phi(y.row(i), mu.col(j), sigma[j]);
      }
      g_ij.row(i)/=(norm);
    }
    //Compute n_j
    n_j = g_ij.colwise().sum();

    return;
  }

  void M_step(void){
    //Update weigths
    w = n_j/n;

    //Update means
    Matrix<double, Dynamic, 1> sum1;
    sum1.resize(d);
    for(int j=0; j<c; j++){
      sum1=MatrixXd::Zero(d,1);
      for(int i=0; i<n; i++){
	sum1 += g_ij(i,j)*y.row(i);
      }
      mu.col(j) = sum1/n_j(j);
    }

    //Update covariances
    double sum2=0;
    for(int j=0; j<c; j++){
      for(int a=0; a<d; a++){
	for(int b=0; b<d; b++){
	  sum2=0;
	  for(int i=0; i<n; i++)
	    sum2 += g_ij(i,j)*(y(i,a)-mu(a,j))*(y(i,b)-mu(b,j));
	  sigma[j](a,b) = sum2/n_j(j);
	}
      }
    }

    return;
  }
};




#endif
