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
    
    if (it == t.end())
    {
      dx = X[i] - x[j];
      dy = Y[i] - y[j];
      dz = Z[i] - z[j];
    }
    else
    {
      j  = it - t.begin();
      r  = 1 - (t[j+1]-T[i])/(t[j+1]-t[j]);
      
      dx = X[i] - (x[j] + (x[j+1] - x[j])*r);
      dy = Y[i] - (y[j] + (y[j+1] - y[j])*r);
      dz = Z[i] - (z[j] + (z[j+1] - z[j])*r);
    }
    
    /*Rcout << std::setprecision(10) << 
      "\nj = " << j << 
      "\nt = " << T[i+1] << 
      "\nx = [" << x[j] << " " << x[j+1] << "]" <<
      "\ny = [" << y[j] << " " << y[j+1] << "]" <<
      "\nz = [" << z[j] << " " << z[j+1] << "]" <<
      "\nt = [" << t[j] << " " << t[j+1] << "]" <<
      "\ndx = " << dx << 
      "\ndy = " << dy << 
      "\ndz = " << dz << 
      "\nR  = " << std::sqrt(dx*dx + dy*dy + dz*dz)<<  std::endl;*/
    
    R[i] = std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  return R;
}

// [[Rcpp::export]]
DataFrame get_sensor_position_at(NumericVector at, DataFrame sensor_pos) 
{
  NumericVector x = sensor_pos["X"];
  NumericVector y = sensor_pos["Y"];
  NumericVector z = sensor_pos["Z"];
  NumericVector t = sensor_pos["gpstime"];
  
  NumericVector X(at.size());
  NumericVector Y(at.size());
  NumericVector Z(at.size());
  NumericVector dt(at.size());
  
  NumericVector::iterator it;
  double dx, dy, dz, r;
  int j;
  
  for (int i = 0 ; i < at.size() ; i++)
  {
    it = std::lower_bound(t.begin(), t.end(), at[i]);
    
    if (it == t.end())
    {
      X[i] = NA_REAL;
      Y[i] = NA_REAL;
      Z[i] = NA_REAL;
      dt[i] = NA_REAL;
    }
    else
    {
      j  = it - t.begin();
      r  = 1 - (t[j+1]-at[i])/(t[j+1]-t[j]);
      
      X[i] = (x[j] + (x[j+1] - x[j])*r);
      Y[i] = (y[j] + (y[j+1] - y[j])*r);
      Z[i] = (z[j] + (z[j+1] - z[j])*r);
      dt[i] = t[j] - at[i];
    }
  }
  
  return  DataFrame::create(_["gpstime"] = at, _["dt"]= dt, _["X"]= X, _["Y"]= Y, _["Z"] = Z);
}