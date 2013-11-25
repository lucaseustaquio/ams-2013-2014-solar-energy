#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
double fn_opt_mae(NumericMatrix x, NumericVector y, NumericVector params) {
    
    int rows = y.size();
    int cols = params.size();
        
    double mae_err = 0.0;
    for (int r = 0; r < rows; ++r) {
      double y_pred = 0.0;
      for (int c = 0; c < cols; ++c) {
        y_pred += params[c]*x(r,c);
      }
      double err = y[r]-y_pred;
      if (err < 0) err = -err;
      mae_err += err;
    }
     
    return mae_err/(double)rows;
}

