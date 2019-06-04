#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector find_range(DataFrame data, DataFrame flightlines) 
{
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  NumericVector T = data["gpstime"];
  
  NumericVector x = flightlines["X"];
  NumericVector y = flightlines["Y"];
  NumericVector z = flightlines["Z"];
  NumericVector t = flightlines["gpstime"];
  
  NumericVector R(X.size());
  NumericVector::iterator it;
  double dx, dy, dz, r;
  int j;
  
  for (int i = 0 ; i < X.size() ; i++)
  {
    it = std::lower_bound(t.begin(), t.end(), T[i]);
    j  = it - t.begin();
    
    r  = 1 - (t[j+1]-T[i])/(t[j+1]-t[j]);
    
    dx = X[i] - (x[j] + (x[j+1] - x[j])*r);
    dy = Y[i] - (y[j] + (y[j+1] - y[j])*r);
    dz = Z[i] - (z[j] + (z[j+1] - z[j])*r);
    
    R[i] = std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  return R;
}