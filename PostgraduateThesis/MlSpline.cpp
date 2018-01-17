#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector ml_tridisolve(NumericVector a, NumericVector b, NumericVector cc,
                            NumericVector d)
{
  NumericVector x = d; double mu = 0.0;
  int n = x.size();
  
  for(int j = 0; j < (n-1); ++j)
  {
    mu = a[j]/b[j];
    b[j+1] = b[j+1] - mu * cc[j];
    x[j+1] = x[j+1] - mu * x[j];
  }
  
  x[n-1] = x[n-1] / b[n-1];
  for(int i = (n-2); i > -1; --i)
  {
    x[i] = (x[i] - cc[i] * x[i+1]) / b[i];
  }
  
  return x;
}

// [[Rcpp::export]]
NumericVector ml_splineslopes(NumericVector h, NumericVector delta)
{
  int n = h.size() + 1;
  NumericVector a(n-1), b(n), cc(n-1), r(n), d(n);
  for(int i = 0; i < (n-2); ++i) //a[1:(n-2)] <- h[2:(n-1)]
  {
    a[i] = h[i+1];
  }
  a[n-2] = h[n-3] + h[n-2];
  
  b[0] = h[1];
  for(int j = 1; j < (n-1); ++j) //b[2:(n-1)] <- 2*(h[2:(n-1)] + h[1:(n-2)])
  {
    b[j] = 2 * (h[j] + h[j-1]);
  }
  b[n-1] = h[n-3];
  
  cc[0] = h[0] + h[1];
  for(int k = 1; k < (n-1); ++k) // cc[2:(n-1)] <- h[1:(n-2)]
  {
    cc[k] = h[k-1];
  }
  
  r[0] = ((h[0] + 2 * cc[0]) * h[1] * delta[0] + pow(h[0], 2) * delta[1]) / cc[0];
  // r[2:(n-1)] <- 3 * (h[2:(n-1)] * delta[1:(n-2)] + h[1:(n-2)] * delta[2:(n-1)])
  for(int l = 1; l < (n-1); ++l)
  {
    r[l] = 3 * (h[l] * delta[l-1] + h[l-1] * delta[l]);
  }
  r[n-1] = (pow(h[n-2], 2) * delta[n-3] + (2*a[n-2] + h[n-2]) * h[n-3] * delta[n-2])/a[n-2];
  
  d = ml_tridisolve(a, b, cc, r);
  return d;
}

// [[Rcpp::export]]
double ml_spline(NumericVector x, NumericVector y, double xi)
{
  if((xi > max(x)) | (xi < min(x)) | (R_IsNA(xi))) return NumericVector::get_na();
  else
  {
    int n = x.size();
    NumericVector h(n-1), delta(n-1), d(n), cc(n-1), b(n-1);
    h  = diff(x);
    delta = diff(y) / h;
    d = ml_splineslopes(h, delta);
    
    for(int i = 0; i < (n-1); ++i)
    {
      //cc <- (3*delta - 2*d[1:(n-1)] - d[2:n])/h
      cc[i] = (3 * delta[i] - 2* d[i] - d[i+1]) / h[i]; 
      //b  <- (d[1:(n-1)] - 2*delta + d[2:n])/h^2
      b[i] = (d[i] - 2 * delta[i] + d[i+1]) / pow(h[i], 2);
    }
    
    int Count = 0;
    for(int j = 1; j < (n-1); ++j)
    {
      if(x[j] <= xi) Count = j;
      else break;
    }
    
    double s = xi - x[Count];
    double v = y[Count] + s * (d[Count] + s * (cc[Count] + s * b[Count]));
    
    return v;
  }
}
