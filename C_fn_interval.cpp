#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List C_fn_interval(int nbellapsed, double bin, int step, int nbpairs, NumericVector somme_cumulative)
{
  IntegerVector seq = seq_len(nbellapsed)-1;
  NumericVector dseq = as<NumericVector>(seq);
  
  NumericVector intervalles_debut = (bin + (double)step) * dseq;
  NumericVector intervalles_fin   = intervalles_debut + bin;
  IntegerVector indexes_debut(nbellapsed);
  IntegerVector indexes_fin(nbellapsed);

  int n = somme_cumulative.length();
  int j = 0;
  int k = 0;
  
  for (int i = 0 ; i < n ; i++)
  {
    if (j < nbellapsed && somme_cumulative[i] > intervalles_debut[j])
    {
      indexes_debut[j] = i+1;
      j++;
    }
    
    if (k < nbellapsed && somme_cumulative[i] > intervalles_fin[k])
    {
      indexes_fin[k] = i+1;
      k++;
    }
  }
  
  return List::create(Rcpp::Named("debut") = indexes_debut, Named("fin") = indexes_fin);
}

// [[Rcpp::export]]
NumericVector fast_cumsum_diff(NumericVector x)
{
  int n = x.length();
  NumericVector y(n-1);
  y[0] = x[1]-x[0];
  
  for(int i = 1 ; i < n-1 ; i++)
  {
    y[i] = y[i-1] + x[i+1] - x[i];
  }
  
  return y;
}

// [[Rcpp::export]]
IntegerVector fast_unique(IntegerVector x)
{
  return unique(x);
}
