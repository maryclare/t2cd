
# Function that takes a value of the fractional differencing parameter d and a time series x
# and returns the log likelihood
fima.ll.auto = function(d, x) {
  # x: time series, allowed to be multivariate p-by-N
  # d: differencing parameter
  
  if (d > 1.5 | d < -2.5) {
    cat("Values of d not in supported range \n")
    break;
  } else {
    
    p = nrow(x)
    mu = rowMeans(x, na.rm = TRUE)
    
    ll = 0
    if (d < 0.5 & d > -0.5) {
      for (i in 1:p){
        ll_i = fima.ll(x[i,] - mu[i], dfrac = d, cond.first = 1)
        ll = ll + ll_i
      }
    } else if (d <= 1.5 & d >= 0.5) {
      for (i in 1:p){
        ll_i = fima.ll(x[i,-1] - x[i,-ncol(x)], 
                        dfrac = d - 1, cond.first = 0)
        ll = ll + ll_i
      }
    } else if (d <= -0.5 & d >= -1.5) {
      for (i in 1:p){
        ll_i = fima.ll(x[i,] - mu[i], dfrac = d + 1, theta = -1,
                        cond.first = 1)        
        ll = ll + ll_i
      }
    } 
    
    return(ll)
    
  }
  
}

### helper functions

# Function that gives conditional means and variances for likelihood computation,
# uses the Durbin-Levinson algorithm
solve.dl = function(cov.fun, z) {
  n = length(cov.fun) 
  C = matrix(0, nrow = n, ncol = n)
  m = v = rep(NA, n)
  gamma.x.0 = cov.fun[1]
  m[1] = 0
  C[1 + 1, 1] = cov.fun[2]/gamma.x.0
  m[2] = C[2, 1]*z[1]
  v[1] = gamma.x.0
  v[2] = v[1]*(1 - C[1 + 1, 1]^2)
  for (i in 2:(n - 1)) {
    C[1 + i, i] = cov.fun[i + 1]
    for (j in 1:(i - 1)) {
      C[1 + i, i] = C[1 + i, i] - C[1 + i-1, j]*cov.fun[abs(i - j) + 1]
    }
    C[1 + i, i] = C[1 + i, i]/v[i]
    for (j in (i-1):1) {
      C[1 + i, j] = C[1 + i - 1, j] - C[1 + i, i]*C[1 + i - 1, i - j]
    }
    m[1 + i] = sum(C[1 + i, 1:i]*z[i:1])
    if (i < n) {
      v[i + 1] = v[i]*(1 - C[1 + i, i]^2)
    }
  }
  
  return(list("m"=m, 
              "v" = v, 
              "C" = C))
}

# Function that returns the autocovariance function of a fractional differencing process
# for specified arguments ks
fi.cv = function(ks, d) {
  
  cvs = numeric(length(ks))
  if (d != 0) {
    for (k in unique(abs(ks))) {
      if (k == 0) {
        cvs[which(abs(ks) == k)] = gamma(1 - 2*d)/(gamma(1 - d)^2)
      } else {
        if (k + d <= 171 & k + 1 - d <= 171) {
          aa = gamma(k + d)
          bb = gamma(k + 1 - d)
          cc = exp(log(aa)-log(bb))
        } else {
          # Used Stirling's formula to make things behave when k is large
          cc = exp(log((k + d - 1)/(k - d))/2 + (k + d - 1)*log(k + d - 1) - (k - d)*log(k - d) - 2*d + 1)
        }
        cvs[which(abs(ks) == k)] = cc*gamma(1 - 2*d)/((gamma(1 - d)*gamma(d)))
      }
    }
  } else {
    cvs[ks == 0] = 1  
  }
  return(cvs)
}

# Function that returns the autocovariance function of a FIMA process, i.e. an ARFIMA(0, d, q) process
fima.cv = function(lag.max, d, theta) {
  cvs = numeric(lag.max + 1)
  for (k in 0:lag.max) {
    nz = which(theta != 0)
    q = length(theta)
    maacf = ARMAacf(ar = 0, ma = theta, lag.max = q)*(1 + sum(theta^2))
    which.nz = c(-nz, 0, nz)
    maacf.nz = maacf[abs(which.nz) + 1]
    cvs[k + 1] = sum(maacf.nz*fi.cv(k - which.nz, d = d))
  }
  return(cvs)
}

# Function that computes the log likelihood of an ARFIMA(0, p, q) process, conditioning
# on the first 'cond.first' observations
fima.ll = function (s, theta = 0, dfrac = 0, cond.first = 0) {
  
  z = s[!is.na(s)]
  n = length(z)
  r = fima.cv(lag.max = n - 1, d = dfrac, theta = theta)
  coef.vars = solve.dl(r, z)
  
  ee = (z - coef.vars$m)[(1 + cond.first):length(coef.vars$m)]
  vv = (coef.vars$v)[(1 + cond.first):length(coef.vars$m)]
  S = sum(ee^2/vv)
  logl = -(length(ee))/2 * log(S/(length(ee))) - sum(log(vv))/2
  return(logl)
}

