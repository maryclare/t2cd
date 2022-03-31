#' @export
# extracts series resistance according to selection of expt, frequency, gel, inf, nor, null, wou
extractdata = function(dat.com, dat.info, expt, freq, gel, inf, nor, null=0, wou=0){
  # dat.com, dat.info: data
  # expt: experiment number
  # freq: frequency
  # gel: gel type
  # inf: whether infected
  # nor: whether normal
  # null: whether well is empty
  # wou: whether wounded
  ind = dat.info$gel==gel & dat.info$inf==inf & dat.info$nor==nor & dat.info$null==null & dat.info$wou==wou
  res = dat.com[expt,'res',ind,freq,]
  tim = dat.com[expt,'tim',ind,freq,]
  # remove rows with all NA
  res = res[!apply(tim, 1, function(x){all(is.na(x))}),]
  tim = tim[!apply(tim, 1, function(x){all(is.na(x))}),]
  return(list(res=res,tim=tim))
}

# simulate series with first segment being a polynomial or gaussian process, and second segment being an ARFIMA model
# noise allowed to have different variance for the two segments, regime 1 noise allowed to be heteroscedastic
# if multivariate, assumes all p sequences change at same location tau+1 if only one tau.percent provided
# if multivariate, observations are aligned
#' @export
sim_simple = function(Tend=70, N=420, tau.percent=0.2, d=0.1, phip=0, piq=0,
                      p=1, pcoeff=NA, deg=5, regime1 = c('poly', 'gp'),
                      sig1=2, sig2=0.5, hetero1=F, share_d=T, seed=0){
  # Tend: ending t
  # N: number of observations
  # tau.percent: percentage of observations in regime 1
  # d, phip, phiq: ARFIMA parameters
  # p: dimensionality of time series
  # pcoeff, deg: polynomial parameters
  # regime1: choice of polynomial and gaussian process for regime 1
  # sig1, sig2: standard deviation of Gaussian noise for regime 1 and 2
  # hetero1: if True, standard deviation of noise in regime 1 varies between 0.1 and sig1
  # share_d: if True, all dimensions have the same d
  set.seed(seed)
  seq_reg2 = NULL
  M = matrix(NA, nrow = p, ncol = N)
  pcoeff.mat = matrix(NA, nrow = p, ncol = deg+1)
  tseq = seq(0, Tend, length.out = N)
  tim = matrix(rep(tseq, p), p, byrow = T)
  if (p>1 & length(tau.percent)==1){
    tau.percent = rep(tau.percent, p)
  }
  tau.idx1 = floor(tau.percent*N)
  tau = tseq[tau.idx1]
  for (i in 1:p){
    tau.idx1_i = tau.idx1[i]
    tau_i = tau[i]
    if (regime1 == 'poly'){
      if (any(is.na(pcoeff))){
        pcoeff_i = rnorm(deg+1, sd = 0.1)
      }else{
        pcoeff_i = pcoeff[i,]
      }
      pcoeff.mat[i,] = pcoeff_i
      for (k in 3:deg){
        # coefficients of degree 3 and up are shrunk
        pcoeff_i[k+1] = pcoeff_i[k+1]*(0.1^k)
      }
      seq_reg1 = rep(pcoeff_i[1], times = tau.idx1_i)
      for (j in 1:deg){
        seq_reg1 = seq_reg1 + pcoeff_i[j+1]*(tseq[1:tau.idx1_i]**j)
      }
    }else if (regime1 == 'gp'){
      seq_reg1 = 10*gaussprocess(m = tau.idx1_i)
    }
    if (phip>0 | piq>0){
      sim = sim.arf(N-tau.idx1_i, d, phip, piq, sig2)
      seq_arfima = sim$s
      armacoeff = list(ar = sim$phicoeff, ma = sim$picoeff)
    }else{
      if (share_d){
        sim = sim.fi(N-tau.idx1_i, d, sig2)
        seq_reg2 = sim$trend + seq_reg1[tau.idx1_i]
      }else{
        sim = sim.fi(N-tau.idx1_i, d+i*0.1, sig2)
        seq_reg2 = sim$trend + seq_reg1[tau.idx1_i]
      }
      seq_arfima = sim$s
      armacoeff = NULL
    }
    if (hetero1){
      # noise sd corresponds to sequence value, mapped to range 0.1 to sig1
      m = (sig1 - 0.1)/(max(seq_reg1) - min(seq_reg1))
      sig1.hetero = m*seq_reg1 + 0.1 - m*min(seq_reg1)
      noise1 = rnorm(tau.idx1_i, 0, sig1.hetero)
    }else{
      noise1 = rnorm(tau.idx1_i, 0, sig1)
    }
    seq_i = c(seq_reg1 + noise1, seq_reg1[tau.idx1_i] + seq_arfima)
    M[i,] = seq_i
  }
  return(list(res=M, tim=tim, pcoeff.mat=pcoeff.mat, tau=tau, d=d,
              armacoeff = armacoeff,
              gt1 = seq_reg1, gt2 = seq_reg2))
}

#' @export
# diffseries but keep mean
diffseries_keepmean = function(x, d){
  # x: input sequence
  # d: differencing parameter
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = length(x), ncol = 1)
  } else if (nrow(x) == 1) {
    x <- t(x)
  }
  N  = nrow(x) # assumes x is N-by-p
  p = ncol(x)
  x[is.na(x)] = 0
  dp_coeff = sapply(0:(N-1), function(k){return(choose(d, k)*(-1)**k)})
  s = matrix(nrow = 0, ncol = p)
  for (i in 1:N){
    s = rbind(s, colSums(matrix(x[i:1,]*dp_coeff[1:i], ncol = p)))
  }
  return(s)
}

#' @export
tr_fi <- function(s, eff.d, mu, dn_coeff = NULL) {
  if (is.null(dn_coeff)) {
    dn_coeff = sapply(1:length(s), function(k){return(choose(eff.d, k)*(-1)**(k+1))})
  }
  mu + sum((s - mu)*dn_coeff[1:length(s)])
}

# simulates univariate FI series
#' @export
sim_fi = function(N, eff.d, sig=1, mu = 0){

  noise = rnorm(N, sd = sig) # ordered from 1 to N
  dn_coeff = sapply(1:N, function(k){return(choose(eff.d, k)*(-1)**(k+1))})
  trend = c(mu)
  s = c(trend[length(trend)] + noise[1])
  for (i in 1:(N-1)){
    trend = c(trend,
              tr_fi(s = s[i:1], eff.d = eff.d, mu = mu, dn_coeff = dn_coeff))
    # trend = c(trend, mu + sum((s[i:1] - mu)*dn_coeff[1:i]))
    s = c(s, trend[length(trend)] + noise[i+1])

  }
  # Note - diffseries_keepmean(s - mu, eff.d) recovers noise
  return(list(s = s, trend = trend))
}

#' @export
trend_fi = function(s, eff.d, mu = 0){

  N <- length(s)
  dn_coeff = sapply(1:N, function(k){return(choose(eff.d, k)*(-1)**(k+1))})
  trend = c(mu)
  for (i in 1:(N-1)){
    trend = c(trend,
              tr_fi(s = s[i:1], eff.d = eff.d,
                    mu = mu, dn_coeff = dn_coeff))

  }
  return(trend)
}

# Don't use unless I have to - need to check it
#' # simulates univariate ARFIMA series; first difference is ARFIMA if d>0.5
#' sim.arf = function(N, eff.d, phip, piq, sig=1){
#'   phicoeff = runif(phip)
#'   picoeff = runif(piq)
#'   if (eff.d>=0.5){
#'     dfrac = eff.d-1
#'     dint = 1
#'   }else{
#'     dfrac = eff.d
#'     dint = 0
#'   }
#'   s = arfima::arfima.sim(N, model = list(phi = phicoeff, dfrac = dfrac, dint = dint,
#'                                          theta = picoeff), sigma2 = sig^2)
#'   return(list(s = s, phicoeff = phicoeff, picoeff = picoeff))
#' }

# generates a gp
# modified from https://www.r-bloggers.com/r-function-for-simulating-gaussian-processes/
#' @export
gaussprocess = function(from = 0, to = 5, K = function(s, t) {10*exp(-0.5*(s-t)**2)},
                        m = 100) {

  tim = seq(from = from, to = to, length.out = m)
  Sigma = sapply(tim, function(s1) {
    sapply(tim, function(s2) {
      K(s1, s2)
    })
  })

  path = MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)

  return(path)
}


