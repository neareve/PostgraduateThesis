#include <Rcpp.h>
#include <iostream>
#include <cmath>
using namespace Rcpp;
using namespace R;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
const double UMAX = 1 - 1e-10;
const double UMIN = 1e-10;
const double XEPS = 1e-4;

// [[Rcpp::export]]
NumericVector Hfunc(int family, NumericVector u, NumericVector v, 
                    NumericVector theta_all, double nu)
{
  int j = 0;
  int n = u.size();
  NumericVector out(n);
  NumericVector h(n);
  double x = 0.0;
  for(int i = 0; i < n; ++i)
  {
    if(u[i] < UMIN) u[i] = UMIN;
    else if(u[i] > UMAX) u[i] = UMAX;
    if(v[i] < UMIN) v[i] = UMIN;
    else if(v[i] > UMAX) v[i] = UMAX;
  }
  
  for(j = 0; j < n; ++j)
  {
    if((v[j] == 0) | (u[j] == 0)) h[j] = 0;
    else if(v[j] == 1) h[j] = u[j];
    else
    {
      double theta = theta_all[j]; 
      if(family == 0) //independent
      {
        h[j] = u[j];
      }
      else if(family == 1) //gaussian
      {
        x = (qnorm(u[j], 0.0, 1.0, 1, 0) - theta * qnorm(v[j], 0.0, 1.0, 1, 0))/
          sqrt(1.0 - pow(theta, 2.0));
        if(isfinite(x)) 
          h[j] = pnorm(x, 0.0, 1.0, 1, 0);
        else if((qnorm(u[j], 0.0, 1.0, 1, 0) - theta * qnorm(v[j], 0.0, 1.0, 1, 0)) < 0)
          h[j] = 0;
        else 
          h[j] = 1;
      }
      else if(family == 2) //student
      {
        double t1, t2, mu, sigma2;
        t1 = qt(u[j], nu, 1, 0); 
        t2 = qt(v[j], nu, 1, 0); 
        mu = theta * t2; 
        sigma2 = ((nu + t2 * t2) * (1.0 - theta * (theta))) / (nu + 1.0);
        h[j] = pt((t1 - mu) / sqrt(sigma2), nu + 1.0, 1, 0);
      }
      else if(family == 3) //clayton
      { 
        if(theta == 0) h[j] = u[j];
        if(theta < XEPS) h[j] = u[j];
        else
        { 
          x = pow(u[j], -theta) + pow(v[j], -theta) - 1.0;
          h[j] =   pow(v[j], -theta - 1.0) * pow(x, -1.0 - 1.0 / (theta));
          if(theta < 0)
          {
            if(x < 0) h[j] = 0; 
          }
        }
      }
      else if(family == 4) //gumbel
      {
        if(theta == 1) h[j] = u[j]; 
        else
        {
          h[j] = -(exp(-pow(pow(-log(v[j]), theta) +
            pow(-log(u[j]), theta), 1.0 / (theta))) * pow(pow(-log(v[j]), theta) +
            pow(-log(u[j]), theta), 1.0 / (theta) - 1.0) *
            pow(-log(v[j]), theta)) / (v[j] * log(v[j]));
        }
      }
      else if(family == 5) //frank
      {
        if(theta == 0) h[j] = u[j];
        else
        {
          h[j] = -(exp(theta) * (exp(theta * u[j]) - 1.0)) /
            (exp(theta * v[j] + theta * u[j]) - exp(theta * v[j] + theta) -
              exp(theta * u[j] + theta) + exp(theta));
        }
      }
      else if(family == 6) //joe
      {
        if(theta == 1) h[j]=u[j];
        else
        {
          h[j] = pow(pow(1.0 - u[j], theta) + pow(1.0 - v[j], theta) - 
            pow(1.0 - u[j], theta) * pow(1.0 - v[j], theta), 1.0/(theta) - 1) * 
            pow(1.0 - v[j], theta - 1.0) * (1 - pow(1 - u[j], theta));
        }
      }
      else if(family == 13) //rotated clayton (180бу)
      {
        if(theta == 0) h[j] = u[j];
        if(theta < XEPS) h[j] = u[j];
        else
        {
          u[j] = 1 - u[j]; 
          v[j] = 1 - v[j];
          x = pow(u[j], -theta) + pow(v[j], -theta) - 1.0;
          h[j] = pow(v[j], -theta - 1.0) * pow(x, -1.0 - 1.0 / (theta)); 
          // pow(v[j],-*theta-1.0)*pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1.0,-1.0-1.0/(*theta));
          h[j] = 1 - h[j];
          u[j] = 1 - u[j];
          v[j] = 1 - v[j];
        }
      }
      else if(family == 14) //rotated gumbel (180бу)
      {
        v[j] = 1 - v[j];
        u[j] = 1 - u[j];
        h[j] = -(exp(-pow(pow(-log(v[j]), theta) + pow(-log(u[j]), theta), 1.0 /
          (theta))) * pow(pow(-log(v[j]), theta) + pow(-log(u[j]), theta), 1.0 /
            (theta) - 1.0) * pow(-log(v[j]), theta)) / (v[j] * log(v[j]));
        h[j] = 1 - h[j];
        u[j] = 1 - u[j];
        v[j] = 1 - v[j];
      }
      else if(family == 16)
      {
        v[j] = 1 - v[j];
        u[j] = 1 - u[j];
        h[j] = pow(pow(1.0 - u[j], theta) + pow(1.0 - v[j], theta) - 
          pow(1.0 - u[j], theta) * pow(1.0 - v[j], theta), 1.0 / (theta) - 1) * 
          pow(1.0 - v[j], theta - 1.0) * (1 - pow(1 - u[j], theta));
        h[j] = 1 - h[j];
        u[j] = 1 - u[j];
        v[j] = 1 - v[j];
      }
    }
    out[j] = max(min(h[j], UMAX), UMIN);
  }
  return out;
}

// [[Rcpp::export]]
double HNumInv(int family, NumericVector u, NumericVector v, NumericVector theta, double nu)
{
  int br = 0;
  double tol = 1e-6, out = 0.0, fl = 0.0, fh = 0.0, val = 0.0;
  NumericVector ans = NumericVector::create(0.0);
  NumericVector x0 = NumericVector::create(UMIN);
  NumericVector x1 = NumericVector::create(UMAX);
  
  fl = Hfunc(family, x0, v, theta, nu)[0]; fl -= u[0];
  fh = Hfunc(family, x1, v, theta, nu)[0]; fh -= u[0];
  if(fabs(fl) <= tol) {ans[0] = x0[0]; br = 1;}
  if(fabs(fh) <= tol) {ans[0] = x1[0]; br = 1;}

  while(!br)
  {
    ans[0] = (x0[0] + x1[0]) / 2.0;
    val = Hfunc(family, ans, v, theta, nu)[0];
    val -= u[0];
    if(fabs(val) <= tol) br = 1;
    if(fabs(x0[0] - x1[0]) <= 1e-10) br = 1;
    if(val > 0.0) {x1[0] = ans[0]; fh = val;}
    else{x0[0] = ans[0]; fl = val;}
  }
  return out = ans[0];
}

// [[Rcpp::export]]
NumericVector Hinv(int family, NumericVector u, NumericVector v, 
                   NumericVector theta_all, double nu)
{
  int n = u.size(), j = 0, i = 0;
  NumericVector hinv(n), out(n);
  for(i = 0; i < n; ++i)
  {
    if(u[i] < UMIN) u[i] = UMIN;
    else if(u[i] > UMAX) u[i] = UMAX;
    if(v[i] < UMIN) v[i] = UMIN;
    else if(v[i] > UMAX) v[i] = UMAX;
  }
  
  for(j = 0; j < n; ++j)
  {
    double theta = theta_all[j]; 
    if(family == 0)
    {
      hinv[j] = u[j];
    }
    else if(family == 1) //gaussian
    {
      hinv[j] = pnorm(qnorm(u[j], 0.0, 1.0, 1,0) *
        sqrt(1.0 - pow(theta, 2.0)) + theta * qnorm(v[j], 0.0, 1.0, 1,0), 0.0, 1.0, 1,0);
    }
    else if(family == 2) //student
    {
      double temp1 = 0.0, temp2 = 0.0, mu = 0.0, var = 0.0;
      temp1 = qt(u[j],nu + 1.0, 1,0); 
      temp2 = qt(v[j],nu, 1, 0); 
      mu = theta * temp2; 
      var = ((nu + (temp2 * temp2)) * (1.0 - (theta * (theta)))) / (nu+1.0);
      hinv[j] = pt((sqrt(var) * temp1) + mu, nu, 1, 0);
    }
    else if(family == 3) //clayton
    {
      if(theta < XEPS) hinv[j] = u[j];
      else
        hinv[j] = pow(pow(u[j] * pow(v[j], theta + 1.0), -theta /
          (theta + 1.0)) + 1.0 - pow(v[j], -theta), -1.0 / (theta));
    }
    else if(family == 4) //gumbel - must turn to numerical inversion
    {
      NumericVector UTmp = NumericVector::create(u[j]);
      NumericVector Vtmp = NumericVector::create(v[j]);
      NumericVector ThetaTmp = NumericVector::create(theta);
      hinv[j] = HNumInv(family, UTmp, Vtmp, ThetaTmp, nu);
    }
    else if(family==5) //frank
    {
      hinv[j] = -1 / (theta) * log(1 - (1 - exp(-theta)) / 
        ((1 / u[j] - 1) * exp(-theta * v[j]) + 1));
    }
    else if(family == 6) //joe - numerical inversion
    {
      NumericVector UTmp = NumericVector::create(u[j]);
      NumericVector Vtmp = NumericVector::create(v[j]);
      NumericVector ThetaTmp = NumericVector::create(theta);
      hinv[j] = HNumInv(family, UTmp, Vtmp, ThetaTmp, nu);
    }
    else if(family == 13)
    {
      u[j] = 1 - u[j];
      v[j] = 1 - v[j];
      hinv[j] = pow(pow(u[j] * pow(v[j], theta + 1.0), -theta /
        (theta + 1.0)) + 1.0 - pow(v[j], -theta), -1.0 / (theta));
      hinv[j] = 1 - hinv[j];
      u[j] = 1 - u[j];
      v[j] = 1 - v[j];
    }
    else if(family == 14) //rotated gumbel (180бу) - must turn to numerical inversion
    {
      int jj = 4;
      u[j] = 1 - u[j];
      v[j] = 1 - v[j];
      NumericVector UTmp = NumericVector::create(u[j]);
      NumericVector Vtmp = NumericVector::create(v[j]);
      NumericVector ThetaTmp = NumericVector::create(theta);
      hinv[j] = HNumInv(jj, UTmp, Vtmp, ThetaTmp, nu);
      hinv[j] = 1 - hinv[j];
      u[j] = 1 - u[j];
      v[j] = 1 - v[j];
    }
    else if(family == 16) //rotated joe (180бу) - must turn to numerical inversion
    {
      int jj = 6;
      u[j] = 1 - u[j];
      v[j] = 1 - v[j];			
      NumericVector UTmp = NumericVector::create(u[j]);
      NumericVector Vtmp = NumericVector::create(v[j]);
      NumericVector ThetaTmp = NumericVector::create(theta);
      hinv[j] = HNumInv(jj, UTmp, Vtmp, ThetaTmp, nu);			
      hinv[j] = 1 - hinv[j];
      u[j] = 1 - u[j];
      v[j] = 1 - v[j];
    }
    out[j] = max(min(hinv[j], UMAX), UMIN); 
  }
  return out;
}