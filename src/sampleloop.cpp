#include <Rcpp.h>
#include <random>  // Include the random library

std::random_device rd;   // Initialize a random device
std::mt19937 g(rd());    // Use Mersenne Twister random number generator

// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix samploop(Rcpp::NumericMatrix a, Rcpp::NumericVector b) {

  // clone vector and matrix to leave alone
  int ar = a.nrow();
  int ac = a.ncol();

  for (int j = 0; j<ac; j++){
    std::shuffle(b.begin(), b.end(), g);
    for (int i = 0; i<ar; i++){
      a(i, j) = b[i];
    }
  }
  return a;
}
