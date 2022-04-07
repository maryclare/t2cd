# loglikelihood, penalty to enforce tau within tau.range
# loglikelihood
#' @export
loglik_res_sigmoid = function(param, tim_cp, tau.idx, x.2, ll.1,
                              pen = FALSE, C = NULL, use_t = FALSE){
  alpha0 = param[1]
  alpha1 = param[2]

  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  # MCG Change: Updated this to make sure that the last weights aren't equal to 0
  wt = c(rep(wt_cp[1], tau.idx[1]-1),
         wt_cp,
         rep(ifelse(wt_cp[length(wt_cp)] != 0, wt_cp[length(wt_cp)], 1),
             length(x.2)-tau.idx[length(tau.idx)]))

  dfrac = param[3]

  if (use_t) {
    df = param[length(param)]
  }

  if ((length(param) == 3 & !use_t)) {
    m <- get.m(x.2 = x.2, dfrac = dfrac, wt = wt)
    sd <- get.sd(x.2 = x.2, dfrac = dfrac, mean = m, wt = wt)
  } else if (length(param) == 4 & !use_t) {
    m <- param[4]
    sd <- get.sd(x.2 = x.2, dfrac = dfrac, mean = m, wt = wt)
  } else if ((length(param) == 5) & !use_t | use_t) {
    m <- param[4]
    sd <- param[5]
  }


  logL = sum((1-wt)*ll.1) +
    ifelse(!use_t,
           loglik_res_step(par = c(dfrac, m, sd), x.2 = x.2, wt = wt),
           loglik_t_res_step(par = c(dfrac, m, sd, df), x.2 = x.2, wt = wt)) +
    ifelse(pen, C*sum(wt_cp[length(wt_cp)] - wt_cp[1]), 0)

  return(logL)
}

# T2CD-sigmoid method
# wrapper around search_dtau_sigmoid to fit with both ranges of d
# and pick the fit with highest likelihood
#' @export
t2cd_sigmoid = function(dat, t.max = 72, tau.range = c(10, 50),
                        init.tau = c(15, 30, 45), deg = 3, C = 1000,
                        seqby = 1, resd.seqby = 5, use_scale = TRUE,
                        use_t = FALSE, prof = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # init.tau: candidate taus to initialize learning
  # deg: degree for B-spline
  # C: regularization coefficient
  # segby, resd.seqby: interval between knots
  # use_scale: if true, scale time series
  res1 = search_dtau_sigmoid(dat, t.max, tau.range, init.tau, deg, C,
                             seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale,
                             use_t = use_t, prof = prof)
  # increase C if change point not in candidate range
  multiplier = 2
  while (is.na(res1$tau)){
    C_new = multiplier*C
    res1 = search_dtau_sigmoid(dat, t.max, tau.range, init.tau, deg, C_new,
                               seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale,
                               use_t = use_t, prof = prof)
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
                               init.tau = c(15, 30, 45), deg = 3, C = 1000,
                               seqby = 1, resd.seqby = 5, use_scale = T,
                               use_t = FALSE, prof = TRUE){
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

  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))
  tim_cp = matrix(tim[,tau.idx], nrow = nrow(tim))

  # optimize for spline component

  fit1 = refitWLS(tim, res_mean, deg = deg,
                  seqby = seqby, resd.seqby = resd.seqby)
  resd1 = res_mean - fit1$fit.vals
  var.resd1 = fit1$var.resd
  ll.1 = dnorm(resd1, log = TRUE, sd = sqrt(var.resd1))

  x.2 = res_mean
  lastN = N

  opt_logL = -Inf
  for (tau_i in init.tau){
    tau_i.idx = which.min(apply(tim <= tau_i, 2, all) == T)
    start <- c(-tau_i, 1, 0)
    lower <- rep(-Inf, length(start))
    upper <- rep(Inf, length(start))
    if (!prof & !use_t) {
      start <- c(start, mean(x.2[tau_i.idx:lastN]), sd(x.2[tau_i.idx:lastN]))
      lower <- c(lower, rep(-Inf, 2))
      upper <- c(upper, rep(Inf, 2))
    }
    if (use_t) {
      start <- c(start, mean(x.2[tau_i.idx:lastN]), sd(x.2[tau_i.idx:lastN]), 3)
      lower <- c(lower, -Inf, -Inf, 2 + 10^(-14))
      upper <- c(upper, Inf, Inf, 300)
    }
    optim_i = optim(par = start,
                    fn = loglik_res_sigmoid, method = "L-BFGS-B",
                    tim_cp = tim_cp, tau.idx = tau.idx, x.2 = x.2,
                    ll.1 = ll.1, pen = TRUE, C = C,
                    control = list("fnscale" = -1), use_t = use_t,
                    lower = lower, upper = upper)
    logL_i = loglik_res_sigmoid(optim_i$par,
                                tim_cp = tim_cp, tau.idx = tau.idx,
                                x.2 = x.2, ll.1 = ll.1, use_t = use_t)
    if (logL_i > opt_logL){
      opt_logL = logL_i
      optim_params = optim_i
    }
  }

  opt_param = optim_params$par
  alpha0 = opt_param[1]
  alpha1 = opt_param[2]
  opt_d = opt_param[3]


  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  # MCG Change: Updated this to make sure that the last weights aren't equal to 0
  wt = c(rep(wt_cp[1], tau.idx[1]-1),
         wt_cp,
         rep(ifelse(wt_cp[length(wt_cp)] != 0, wt_cp[length(wt_cp)], 1),
             length(x.2)-tau.idx[length(tau.idx)]))
  opt_tau.idx = which(wt>=0.5, arr.ind = TRUE)[1]-1
  opt_tau = tim[opt_tau.idx[1]]

  # range of change location
  opt_taurange1.idx = which(wt>=0.1, arr.ind = TRUE)[1]-1
  opt_taurange2.idx = which(wt>=0.9, arr.ind = TRUE)[1]-1
  opt_taurange1 = tim[opt_taurange1.idx[1]]
  opt_taurange2 = tim[opt_taurange2.idx[1]]
  opt_taurange1[is.na(opt_taurange1)] = tau.range[1]
  opt_taurange2[is.na(opt_taurange2)] = tau.range[2]

  if (use_t) {
    opt_m <- opt_param[4]
    opt_sd <- abs(opt_param[5])
    opt_df <- opt_param[6]
  } else {
    if (prof) {
      opt_m <- get.m(x.2 = x.2, dfrac = opt_d, wt = wt)
      opt_sd <- get.sd(x.2 = x.2, dfrac = opt_d, mean = opt_m, wt = wt)
    } else if (!prof) {
      opt_m <- opt_param[4]
      opt_sd <- abs(opt_param[5])
    }
    opt_df = NA
  }

  return(list(res = res, tim = tim, tau.idx = tau.idx, m = opt_m,
              sd = opt_sd, df = opt_df,
              d = opt_d, tau = opt_tau, idx = opt_tau.idx,
              tau.range1 = opt_taurange1, tau.range2 = opt_taurange2,
              param = opt_param, logL = opt_logL))
}

# plot sequences and fitted lines
#' @export
plot_t2cd_sigmoid = function(results, tau.range = c(10, 50), deg = 3,
                             seqby = 1, resd.seqby = 5, return_plot = TRUE){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
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
  # MCG Change: Updated this to make sure that the last weights aren't equal to 0
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(ifelse(wt_cp[,ncol(wt_cp)] != 0, wt_cp[,ncol(wt_cp)], 1),
                        N-tau.idx[length(tau.idx)]), p, byrow = F))

  # update variables if using original or first difference
  x.2 = res_mean
  d = opt_d
  diff_p = t(diffseries_keepmean(t(wt*(x.2-m)), d))

  mu.2 = wt*(x.2-m) - diff_p
  fit.vals = (mu.2 + wt*m)*attributes(res_mean)$'scaled:scale' +
  (1-wt)*fit1*attributes(res_mean)$'scaled:scale'

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
