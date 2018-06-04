# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' An RcppArmadillo function to compute the
//' determinant of a Toeplitz matrix
//'
//' @param c A numeric vector, which is the first row of the Toeplitz matrix
//' @return A list containing the log determinant and two vectors z and v
//' to be passed to trenchInv to calculate the inverse.
// [[Rcpp::export]]
List trenchDetcpp(arma::vec c) {
  
  // Create storage
  int N = c.size();
  arma::vec xi(N-1);
  arma::vec z(N-1);
  arma::vec v(N);
  double beta = 1;
  double alpha;
  double l = 0;
  double logdet;
  
  // preliminary computations
  xi = c.subvec(1, N-1)/c(0);
  z(0) =  -xi(0);
  alpha = -xi(0);
  int K = xi.size(); 
  
  // Durbins Algorithm
  for (int i = 0; i < (K - 1); i++) {
      beta = (1 - pow(alpha, 2)) * beta ;
      l = l + log(beta) ;
      if (i == 0) {
        alpha = - (xi(i + 1) + xi(0) * z(0)) / beta ;
        z(0) = z(0) + alpha * z(0) ;
      } else {
        alpha = - (xi(i + 1) + arma::dot(flipud(xi.subvec(0, i)), z.subvec(0, i)) ) / beta ;
        z.subvec(0, i) = z.subvec(0, i) + alpha * flipud(z.subvec(0, i)) ;
      }
      z(i+1) = alpha  ;
  }
  
  beta = (1 - pow(alpha, 2)) * beta ;
  l = l + log(beta) ;
  
  logdet = l + N * log(c(0));
  
  v(N-1) = 1 / ((1 + arma::dot(xi, z)) * c(0)) ;
  v.subvec(0,N-2) = v(N-1) * flipud(z.subvec(0, N - 2));
  
  return List::create(Rcpp::Named("logdet") = logdet,
                      Rcpp::Named("z") = z,
                      Rcpp::Named("v") = v);
}
//' An RcppArmadillo function to compute the
//' inverse of a Toeplitz matrix
//'
//' @param v A numeric vector, which is calculated from trenchDetcpp
//' @return The inverse of the toeplitz matrix
// [[Rcpp::export]]
arma::mat trenchInvcpp(arma::vec v) {
  
  //storage
  int N = v.size();
  int i;
  int j;
  arma::mat C(N, N);
  arma::mat trenchInv;
  
  //Durbin Algorithm
  C.row(0) = flipud(v).t();
  C.col(0) = flipud(v);
  C.row(N - 1) = v.t();
  C.col(N - 1) = v;
  for(i = 1; i < floor( (N - 1) / 2 ) + 1; i++){
    for(j = 1; j < N - i; j++){
      C(i, j) = C(i - 1, j - 1) + (v(N - 1 - j) * v(N - 1 - i) - v(i - 1) * v(j - 1)) / v(N - 1) ;
      C(j, i) = C(i, j);
      C(N - 1 - i , N - 1 - j ) = C(i, j);
      C(N - 1 - j , N - 1 - i ) = C(i, j);
    }
  } 
  
  trenchInv = C;
  return(trenchInv) ; 
  
}
