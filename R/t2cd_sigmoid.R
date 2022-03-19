# loglikelihood, penalty to enforce tau within tau.range
#' @export
negloglik_pen_res_sigmoid = function(param, tim_cp, tau.idx, N, p, x.2, dflag,
                                     ll.1.mat, C){
  alpha0 = param[1]
  alpha1 = param[2]
  m = param[3]
  dfrac = param[4]

  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))

  if (dflag == 'original'){
    x.2m = x.2-m
  }else{
    x.2m = cbind(rep(0, p), x.2-m)
  }
  diff_p = t(diffseries_keepmean(t(wt*(x.2m)), dfrac))

  neglogL = -sum((1-wt)*ll.1.mat) + 0.5*log(2*pi)*sum(wt) + 0.5*sum(wt) +
    sum(0.5*rowSums(wt)*log(rowSums(wt*diff_p^2)/rowSums(wt))) -
    p*C*sum(wt_cp[,ncol(wt_cp)] - wt_cp[,1])

  return(neglogL)
}

# loglikelihood
#' @export
loglik_res_sigmoid = function(param, tim_cp, tau.idx, wt_cp, N, p, x.2, dflag, ll.1.mat){
  alpha0 = param[1]
  alpha1 = param[2]
  m = param[3]
  dfrac = param[4]

  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))

  if (dflag == 'original'){
    x.2m = x.2-m
  }else{
    x.2m = cbind(rep(0, p), x.2-m)
  }
  diff_p = t(diffseries_keepmean(t(wt*(x.2m)), dfrac))

  logL = sum((1-wt)*ll.1.mat) - 0.5*log(2*pi)*sum(wt) - 0.5*sum(wt) -
    sum(0.5*rowSums(wt)*log(rowSums(wt*diff_p^2)/rowSums(wt)))

  return(logL)
}

# T2CD-sigmoid method
# wrapper around search_dtau_sigmoid to fit with both ranges of d
# and pick the fit with highest likelihood
#' @export
t2cd_sigmoid = function(dat, t.max = 72, tau.range = c(10, 50),
                        init.tau = c(15, 30, 45), deg = 3, C = 1000,
                        seqby = 1, resd.seqby = 5, use_scale = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # init.tau: candidate taus to initialize learning
  # deg: degree for B-spline
  # C: regularization coefficient
  # segby, resd.seqby: interval between knots
  # use_scale: if true, scale time series
  res1 = search_dtau_sigmoid(dat, t.max, tau.range, init.tau, deg, C, dflag = 'original',
                             seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale)
  # increase C if change point not in candidate range
  multiplier = 2
  while (is.na(res1$tau)){
    C_new = multiplier*C
    res1 = search_dtau_sigmoid(dat, t.max, tau.range, init.tau, deg, C_new, dflag = 'original',
                               seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale)
    multiplier = multiplier + 1
  }
  return(res1)
}

### helper and plotting functions

# sigmoid
#' @export
sigmoid = function(a){return(1/(1+exp(-a)))}
#' @export
softmax = function(a,b){
  return(exp(a)/(exp(a) + exp(b)))
}

# optimize the likelihood for d and tau
# option to initialize tau at multiple indices
#' @export
search_dtau_sigmoid = function(dat, t.max = 72, tau.range = c(10, 50),
                               init.tau = c(15, 30, 45), deg = 3, C = 1000, dflag = 'fdiff',
                               seqby = 1, resd.seqby = 5, use_scale = T){
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
  N = ncol(res_mean)
  p = nrow(res_mean)

  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))
  tim_cp = matrix(tim[,tau.idx], nrow = nrow(tim))

  tol = 0.01
  # determine range of d to search
  if (is.null(dflag)){
    tau_last = tau.idx[length(tau.idx)]
    dfrac_p = rep(0, p)
    for (k in 1:p){
      x.last = res_mean[k,(tau_last+1):N]
      arf = arfima::arfima(x.last)
      dfrac_p[k] = arf$modes[[1]]$dfrac
    }
    if ((0.5-mean(dfrac_p)) < tol){
      dflag = 'fdiff' # likelihood evaluated on first differences
    }else{
      dflag = 'original' # likelihood evaluated on original series
    }
  }

  # optimize for spline component
  ll.1.mat = matrix(ncol = N, nrow = 0)
  for (k in 1:p){
    fit1 = refitWLS(tim[k,], res_mean[k,], deg = deg,
                    seqby = seqby, resd.seqby = resd.seqby)
    resd1 = res_mean[k,] - fit1$fit.vals
    var.resd1 = fit1$var.resd
    ll.1.vec = dnorm(resd1, log = TRUE, sd = sqrt(var.resd1))
    ll.1.mat = rbind(ll.1.mat, ll.1.vec)
  }

  if (dflag == 'original'){
    x.2 = res_mean
    lastN = N
  }else{
    x.2 = t(diff(t(res_mean), 1))
    lastN = N-1
  }

  opt_logL = -Inf
  for (tau_i in init.tau){
    tau_i.idx = which.min(apply(tim <= tau_i, 2, all) == T)
    optim_i = optim(par = c(-tau_i, 1, mean(x.2[tau_i.idx:lastN]), 0),
                    fn = negloglik_pen_res_sigmoid, method = "BFGS",
                    tim_cp = tim_cp, tau.idx = tau.idx, N = N, p = p, x.2 = x.2,
                    dflag = dflag, ll.1.mat = ll.1.mat, C = C)
    logL_i = loglik_res_sigmoid(optim_i$par,
                                tim_cp = tim_cp, tau.idx = tau.idx, N = N,
                                p = p, x.2 = x.2, dflag = dflag, ll.1.mat = ll.1.mat)
    if (logL_i > opt_logL){
      opt_logL = logL_i
      optim_params = optim_i
    }
  }

  opt_param = optim_params$par
  alpha0 = opt_param[1]
  alpha1 = opt_param[2]
  if (dflag == 'original'){
    opt_d = opt_param[4]
  }else{
    opt_d = opt_param[4] + 1
  }

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

  # range of change location
  opt_taurange1.idx = apply(wt, 1, function(x){return(which(x>=0.1, arr.ind = TRUE)[1])})
  opt_taurange2.idx = apply(wt, 1, function(x){return(which(x>=0.9, arr.ind = TRUE)[1]-1)})
  opt_taurange1 = c()
  opt_taurange2 = c()
  for (k in 1:p){
    opt_taurange1 = c(opt_taurange1, tim[k,opt_taurange1.idx[k]])
    opt_taurange2 = c(opt_taurange2, tim[k,opt_taurange2.idx[k]])
  }
  opt_taurange1[is.na(opt_taurange1)] = tau.range[1]
  opt_taurange2[is.na(opt_taurange2)] = tau.range[2]

  return(list(res = res, tim = tim, tau.idx = tau.idx,
              d = opt_d, tau = opt_tau, m = opt_param[3], idx = opt_tau.idx,
              tau.range1 = opt_taurange1, tau.range2 = opt_taurange2,
              param = opt_param, logL = opt_logL, dflag = dflag))
}

# plot sequences and fitted lines
#' @export
plot_t2cd_sigmoid = function(results, tau.range = c(10, 50), deg = 3,
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
  alpha0 = opt_param[1]
  alpha1 = opt_param[2]
  m = opt_param[3]

  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(ifelse(wt_cp[,ncol(wt_cp)] != 0, wt_cp[,ncol(wt_cp)], 1),
                        N-tau.idx[length(tau.idx)]), p, byrow = F))

  # update variables if using original or first difference
  if (dflag == 'original'){
    x.2 = res_mean
    d = opt_d
    diff_p = t(diffseries_keepmean(t(wt*(x.2-m)), d))
  }else{
    x.2 = cbind(0, t(diff(t(res_mean), 1))-m)
    d = opt_d - 1
    diff_p = t(diffseries_keepmean(t(wt*x.2), d))
  }

  if (dflag == 'original'){
    mu.2 = wt*(x.2-m) - diff_p
    fit.vals = (mu.2 + wt*m)*attributes(res_mean)$'scaled:scale' +
      (1-wt)*fit1*attributes(res_mean)$'scaled:scale'
  }else{
    diff.2 = wt*x.2 - diff_p
    mu.2 = matrix(nrow = p, ncol = 0)
    for (i in 1:N){
      mu.2 = cbind(mu.2, res_mean[,i] + diff.2[,i])
    }
    fit.vals = (wt*mu.2)*attributes(res_mean)$'scaled:scale' +
      (1-wt)*fit1*attributes(res_mean)$'scaled:scale'
  }

  # plotting
  if (return_plot){
    plot(tim[1,], res[1,], ylim = c(min(c(res, fit.vals)), max(c(res, fit.vals))), type = 'l',
         main = paste('Values fitted with d: ', round(opt_d,3), ' tau: ', round(mean(opt_tau),3)),
         xlab = 'Time (hour)', ylab = 'Resistance (ohm)')
    if (p > 1){
      for (k in 2:p){
        lines(tim[k,], res[k,])
      }
    }
    if (is.na(opt_tau[1])){
      lines(tim[1,], fit.vals[1,], col = "blue", lwd = 1)
    }else{
      opt_tau.idx = which(tim[1,] == opt_tau[1])
      lines(tim[1,1:opt_tau.idx], fit.vals[1,1:opt_tau.idx], col = "blue", lwd = 1)
      lines(tim[1,(opt_tau.idx):ncol(fit.vals)], fit.vals[1,(opt_tau.idx):ncol(fit.vals)], col = "green", lwd = 1)
      abline(v = opt_tau[1], lty = 2, col = "red")
    }
    if (p > 1){
      for (k in 2:p){
        if (is.na(opt_tau[k])){
          lines(tim[k,], fit.vals[k,], col = "blue", lwd = 1)
        }else{
          opt_tau.idx = which(tim[k,] == opt_tau[k])
          lines(tim[k,1:opt_tau.idx], fit.vals[k,1:opt_tau.idx], col = "blue", lwd = 1)
          lines(tim[k,(opt_tau.idx):ncol(fit.vals)], fit.vals[k,(opt_tau.idx):ncol(fit.vals)], col = "green", lwd = 1)
          abline(v = opt_tau[k], lty = 2, col = "red")
        }
      }
    }
    abline(v = tau.range, lty = 1, col = "red")
  }

  if (p ==1){
    opt_tau.idx = results$idx
    return(list(fit.vals = fit.vals,
                fit.vals1 = fit1, fit.vals2 = mu.2,
                var.resd1 = var.resd1,
                wt = wt, scaling = attributes(res_mean)$'scaled:scale'))
  }else{
    return(list(fit.vals = fit.vals,
                var.resd1 = var.resd1*attributes(res_mean)$'scaled:scale'^2,
                wt = wt, scaling = attributes(res_mean)$'scaled:scale'))
  }
}

# parametric bootstrap using outputs from t2cd_step and plot.t2cd_step
#' @export
bootstrap_sample_sigmoid = function(results, plot_results, seed = 0){

  set.seed(seed)
  res = results$res
  tim = results$tim
  N = length(res)

  # regime 1
  fit.vals1 = plot_results$fit.vals1*plot_results$scaling
  var.resd1 = plot_results$var.resd1*plot_results$scaling^2
  noise1 = rnorm(N, 0, sqrt(var.resd1))

  # regime 2
  opt_d = results$d
  m = results$m
  wt = plot_results$wt
  x.2 = t(scale(t(res), center = F)) # scaling
  diff_p = t(diffseries_keepmean(t(wt*(x.2-m)), opt_d))
  sd.resd2 = sqrt(rowSums(wt*diff_p^2)/rowSums(wt))
  sim = sim.fi(N, opt_d, sd.resd2, mu = m)
  seq_fi = sim$s

  samp = c((1-wt)*(fit.vals1 + noise1) + wt*(seq_fi)*plot_results$scaling)

  return(list(res=matrix(samp, nrow=1), tim=tim))
}
