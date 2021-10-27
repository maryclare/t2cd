# multivariate implementation for T2CD-sigmoid
# optimize the likelihood for d and tau
# option to initialize tau at multiple indices
#' @export
t2cd_sigmoid.mv = function(dat, t.max = 72, tau.range = c(10, 50),
                           init.tau = c(15, 30, 45), deg = 3, C = 1000,
                           seqby = 1, resd.seqby = 5, use_scale = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # init.tau: candidate taus to initialize learning
  # deg: degree for B-spline
  # C: regularization coefficient
  # seqby, resd.seqby: interval between knots
  # use_scale: if true, scale time series

  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }

  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)

  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  if (use_scale){
    res_mean = t(scale(t(res), center = F)) # scaling
  }else{
    res_mean = res
  }
  N = ncol(res)
  p = nrow(res)

  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))
  tim_cp = matrix(tim[,tau.idx], nrow = nrow(tim))

  # optimize for spline component
  ll.1.mat = matrix(ncol = N, nrow = 0)
  for (k in 1:p){
    fit1 = refitWLS(tim[k,], res[k,], deg = deg)
    resd1 = res[k,] - fit1$fit.vals
    var.resd1 = fit1$var.resd
    ll.1.vec = dnorm(resd1, log = TRUE, sd = sqrt(var.resd1))
    ll.1.mat = rbind(ll.1.mat, ll.1.vec)
  }

  # determine range of d to search and find initializing parameters
  init.param = matrix(NA, p, 4)
  init.d = rep(NA, p)
  for (k in 1:p){
    res_k = t2cd_sigmoid(list(res=matrix(res[k,], 1), tim=matrix(tim[k,], 1)),
                         t.max, tau.range, init.tau, deg, C,
                         seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale)
    init.param[k,] = res_k$par
    init.d[k] = res_k$d
  }

  dflag = 'original'
  x.2 = res_mean
  init.d.2 = init.d

  # loglikelihood, penalty to enforce tau within tau.range
  negloglik_partial_pen = function(param){
    m = param[1:p]
    dfrac = param[p+1]

    # weights
    wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
    wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
               wt_cp,
               matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))

    x.2m = x.2-m
    diff_p = t(diffseries_keepmean(t(wt*(x.2m)), dfrac))

    neglogL = -sum((1-wt)*ll.1.mat) + 0.5*log(2*pi)*sum(wt) + 0.5*sum(wt) +
      sum(0.5*rowSums(wt)*log(rowSums(wt*diff_p^2)/rowSums(wt))) -
      p*C*sum(wt_cp[,ncol(wt_cp)] - wt_cp[,1])

    return(neglogL)
  }

  # loglikelihood
  loglik = function(param){
    m = param[1:p]
    dfrac = param[p+1]

    # weights
    wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
    wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
               wt_cp,
               matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))

    x.2m = x.2-m
    diff_p = t(diffseries_keepmean(t(wt*(x.2m)), dfrac))

    logL = sum((1-wt)*ll.1.mat) - 0.5*log(2*pi)*sum(wt) - 0.5*sum(wt) -
      sum(0.5*rowSums(wt)*log(rowSums(wt*diff_p^2)/rowSums(wt)))

    return(logL)
  }

  alpha0 = init.param[,1]
  alpha1 = init.param[,2]
  m = init.param[,3]
  optim_params = optim(par = c(m, mean(init.d.2)),
                       fn = negloglik_partial_pen, method = "BFGS")

  opt_param = c(alpha0, alpha1, optim_params$par)
  opt_logL = loglik(optim_params$par)
  opt_d = opt_param[3*p+1]
  univ_d = init.param[,4]

  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))
  opt_tau.idx = apply(wt, 1, function(x){return(which(x>=0.5, arr.ind = TRUE)[1]-1)})
  opt_tau = c()
  for (k in 1:p){
    opt_tau = c(opt_tau, tim[k,opt_tau.idx[k]])
  }

  return(list(res = res, tim = tim, tau.idx = tau.idx, idx = opt_tau.idx,
              d = opt_d, univ_d = univ_d, tau = opt_tau, param = opt_param, logL = opt_logL,
              dflag = dflag))
}

# plot sequences and fitted lines
#' @export
plot.t2cd_sigmoid.mv = function(results, tau.range = c(10, 50), deg = 3,
                             seqby = 1, resd.seqby = 5, return_plot = TRUE){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
  dflag = results$dflag
  res_mean = t(scale(t(res), center = F)) # scaling
  N = ncol(res_mean)
  p = nrow(res_mean)

  fit1 = var.resd1 = matrix(nrow = 0, ncol = N)
  for (k in 1:p){
    fitwls = refitWLS(tim[k,], res_mean[k,], deg = deg,
                      seqby = seqby, resd.seqby = resd.seqby)
    fit1 = rbind(fit1, fitwls$fit.vals)
    var.resd1 = rbind(var.resd1, fitwls$var.resd)
  }

  # points in change range
  tim_cp = matrix(tim[,tau.idx], nrow = nrow(tim))

  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_param = results$param

  # fitted values
  alpha0 = opt_param[1:p]
  alpha1 = opt_param[p + 1:p]
  m = opt_param[2*p + 1:p]

  # weights
  wt_cp <- matrix(nrow = nrow(tim_cp), ncol = ncol(tim_cp))
  var.resd <- fit.vals <- wt <- matrix(nrow = nrow(tim), ncol = ncol(tim))
  for (i in 1:p) {
    wt_cp[i, ] = sigmoid(alpha0[i]+alpha1[i]*tim_cp[i, ])
    wt[i, ] <- c(rep(wt_cp[i,1], tau.idx[1]-1),
                 wt_cp[i, ],
                 rep(wt_cp[i,ncol(wt_cp)], N-tau.idx[length(tau.idx)]))


  # update variables if using original or first difference
  if (dflag == 'original'){
    x.2 = res_mean[i, ]
    d = opt_d
    diff_p = t(diffseries_keepmean(t(wt[i, ]*(x.2-m[i])), d))
  }else{
    x.2 = c(0, t(diff(t(res_mean[i, ]), 1))-m[i])
    d = opt_d - 1
    diff_p = t(diffseries_keepmean(t(wt[i, ]*x.2), d))
  }

  if (dflag == 'original'){
    mu.2 = wt[i, ]*(x.2-m[i]) - diff_p
    fit.vals[i, ] = (mu.2 + wt[i, ]*m[i])*attributes(res_mean)$'scaled:scale'[i] +
      (1-wt[i, ])*fit1[i, ]*attributes(res_mean)$'scaled:scale'[i]
  }else{
    diff.2 = wt[i, ]*x.2 - diff_p
    mu.2 = matrix(nrow = 1, ncol = 0)
    for (j in 1:N){
      mu.2 = c(mu.2, res_mean[i,j] + diff.2[j])
    }
    fit.vals[i, ] = (wt[i, ]*mu.2)*attributes(res_mean)$'scaled:scale'[i] +
      (1-wt[i, ])*fit1[i, ]*attributes(res_mean)$'scaled:scale'[i]
  }
  var.resd[i, ] <- var.resd1[i, ]*attributes(res_mean)$'scaled:scale'[i]^2
}

  return(list(fit.vals1 = fit.vals,
              var.resd1 = var.resd,
              wt = wt, N = N))
}
