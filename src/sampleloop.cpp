#include <Rcpp.h>

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix samploop(Rcpp::NumericMatrix a, Rcpp::NumericVector b, int c) {

    // clone vector and matrix to leave alone
	int ar = a.nrow();
	int ac = a.ncol();

		for (int j = 0; j<ac; j++){
		  std::random_shuffle(b.begin(), b.end(), randWrapper);
			for (int i = 0; i<ar; i++){
				a(i, j) = b[i];
			}
		}
		return a;
	}
