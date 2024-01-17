#include <iostream>
#include "math.h"
#include <Rcpp.h>
#include "Algos.h"

using namespace Rcpp;
using namespace std;
/* OP: method = "OP" 
 * PELT: method = "PELT" 
 * GeomFPOP(R-type:random/random)  method = "GeomFPOP"
 * chFPOP: method = "chFPOP"
 */

//converting parameter ("method") to a numeric value.
unsigned int typeAlgo(std::string method) {
  unsigned int type_algo = INFINITY;
  if (method == "OP") { type_algo = 0; }
  if (method == "PELT") { type_algo = 1; }
  if (method == "GeomFPOP") { type_algo = 2; }
  if (method == "chFPOP") { type_algo = 3; }
  return type_algo;
}

//' @title get_changepoints
//'
//' @description Multiple Changepoint Detection using methods: OP, PELT, GeomFPOP.
//' @param data is a matrix of data (p-rows x n-columns, where p - dimension, n - length).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param method is the algorithm: 'OP', 'PELT', 'GeomFPOP (rectangle:random/random)' or 'chFPOP'.
//' @param showNbCands is the logical parameter (if "true", than to show the number of candidates at each iteration).
//' @return a list of  elements  = (changes, means, UnpenalizedCost, NumberOfCandidates).
//'
//' \describe{
//' \item{\code{changes}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
//' \item{\code{means}}{is the list of successive means for the p-variate time series.}
//' \item{\code{UnpenalizedCost}}{is a number equal to the global cost.}
//' \item{\code{NumberOfCandidates}}{is a number of candidates at each iteration (vector).}
//' }
//'
//' @examples
//' N <- 100000
//' Chpt <-5000
//' Means <-  matrix(c(0,1,1,10), nrow = 2) 
//' Noise <- 1
//' Dim <- 2
//' Penalty <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
//' Dim <- 3
//' Penalty <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = 3, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 3), noise = 1)
//' get_changepoints(data = time_series, penalty = Penalty, method = 'OP',showNbCands = FALSE)
//' get_changepoints(data = time_series, penalty = Penalty, method = 'PELT', showNbCands = FALSE)
//' get_changepoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', showNbCands = FALSE)
//' get_changepoints(data = time_series, penalty = Penalty, method = 'chFPOP', showNbCands = FALSE)

// [[Rcpp::export]]
List get_changepoints(Rcpp::NumericMatrix data, double penalty, std::string method = "GeomFPOP", bool showNbCands = false) {
  unsigned int type_algo = typeAlgo(method);
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {throw std::range_error("Penalty should be a non-negative number!");}
  if(type_algo == INFINITY){throw std::range_error("This combination of parameters is not available.");}
//  if( ((unsigned int)data.nrow() < 1) ||  ((unsigned int)data.nrow() > 20)) {throw std::range_error("The dimension of time series can not exceed 20.");}
  //----------------------------------------------------------------------------
  unsigned int p = (unsigned int)data.nrow();
  
  if (p == 1){
    Algos<1> X = Algos<1>(data, penalty);
    return X.algosOP(type_algo, showNbCands);
  } else 
  if (p == 2){
    Algos<2> X = Algos<2>(data, penalty);
    return X.algosOP(type_algo, showNbCands);
  } else if (p == 3){
    Algos<3> X = Algos<3>(data, penalty);
    return X.algosOP(type_algo, showNbCands);
  } else if (p == 4){
    Algos<4> X = Algos<4>(data, penalty);
    return X.algosOP(type_algo, showNbCands);
  } 
  return NULL;
}
