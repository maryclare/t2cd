# grid search for tau given FI parameters
# return the likelihood for each combination
#' @export
search_tau_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                           dflag = 'original', seqby = 1, resd.seqby = 5,
                           use_arf = FALSE, use_scale = TRUE,
                           fix.d = NULL, fix.m = NULL){
  stopifnot(use_arf==FALSE)

  # select data below t.max
  if (is.na(t.max)){
    t.max = max(dat$tim, na.rm = T)
  }

  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(tim.ind == T)

  res = dat$res[t.maxidx]
  if (use_scale){
    res_mean = scale(res[t.maxidx], center = F) # scaling
  }else{
    res_mean = matrix(res[t.maxidx], ncol = 1)
  }
  tim = dat$tim[t.maxidx]
  N = length(res_mean)

  # initialize result vectors
  tau.idx = which(tim >= tau.range[1] & tim <= tau.range[2])
  M = c()
  d = c()
  m = c()

  # iterate through each tau, return log-likelihood
  for (j in 1:length(tau.idx)){
    tau_j = tau.idx[j]

    # optimize for polynomial component
    x.1 = res_mean[1:tau_j]
    t.1 = tim[1:tau_j]
    n.1 = length(x.1)

    fit1 = refitWLS(t.1, x.1, deg = deg, seqby = seqby, resd.seqby = resd.seqby)
    resd1.1 = x.1 - fit1$fit.vals
    var.resd1.1 = fit1$var.resd
    ll.1 = sum(dnorm(resd1.1, log = TRUE, sd = sqrt(var.resd1.1)))

    # optimize for ARFIMA
    if (dflag == 'original'){
      x.2 = res_mean[(tau_j+1):N]
    }else{
      x.2 = diff(res_mean[(tau_j):N], 1)
    }
    # optimizing (why are we optimizing again here??)
    optim.2 = optim(par = c(fix.m, fix.d),
                    fn = negloglik_res_step, method = "BFGS", x.2 = x.2)
    if (dflag == 'original'){
      d = c(d, optim.2$par[2])
    }else{
      d = c(d, optim.2$par[2] + 1)
    }
    m = c(m, optim.2$par[1])
    ll.2 = loglik_res_step(optim.2$par, x.2 = x.2)

    M = c(M, ll.1 + ll.2)
  }

  # tau and d at maximum log-likelihood
  M_df = data.frame(tau = tim[tau.idx], M = M, d = d, m = m)
  max.idx = which.max(M)
  max.tau = tim[tau.idx[max.idx]]
  max.d = d[max.idx]
  max.m = m[max.idx]

  return(list(M_df = M_df, res = res, tim = tim, tau.idx = tau.idx,
              tau = max.tau, d = max.d, m = max.m,
              idx = tau.idx[max.idx], logL = M[max.idx],
              dflag = dflag))

}

#' @export
get.m <- function(x.2, dfrac) {
  optim(x.2[1], function(m, x.2, dfrac) {
    diff_p = t(diffseries_keepmean(t(x.2-m), dfrac))
    sum(diff_p^2)
  }, x.2 = x.2, dfrac = dfrac, method = "Brent",
  lower = min(x.2), upper = max(x.2))$par
}

# helper functions for loglikelihood
#' @export
negloglik_res_step_mv = function(param, x.2){
  dfrac = param

  n.2 = rowSums(!is.na(x.2))
  m <- rep(NA, nrow(x.2))
  for (k in 1:length(m)) {
    m[k] <- get.m(x.2 = na.omit(x.2[k, ]), dfrac = dfrac)
  }
  diff_p = t(diffseries_keepmean(t(x.2-m), dfrac))

  diff_p[diff_p==0] = NA
  neglogL = sum(n.2*log(rowMeans(diff_p^2, na.rm = T)))

  return(neglogL)
}

# multivariate implementation for T2CD-step
# grid search for d and tau
#' @export
t2cd_step_mv = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                        seqby = 1, resd.seqby = 5,
                        use_arf = TRUE, use_scale = TRUE,
                        maxiter=10, tol=1e-6){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # deg: degree for B-spline
  # seqby, resd.seqby: interval between knots
  # use_arf: if true, use arfima estimates from arfima package
  # use_scale: if true, scale time series
  # maxiter: maximum number of iterations solving for tau and FI parameters
  # tol: stops iterating if difference between the most recent objective values is below tol

  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }

  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)

  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  N = ncol(res)
  p = nrow(res)

  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))

  # alternate between optimizing for tau and FI parameters
  iter_k = 1
  iter_flag = TRUE
  neglogL_prev = Inf
  init.d = init.m = fix.d = fix.m = init.tau = init.idx = tau = idx = rep(NA, p)
  hist.d = hist.neglogL = c()

  while (iter_flag){
    # retain all data in regime 2
    reg2 = matrix(NA, nrow = p, ncol = N)
    for (k in 1:p){
      cat("k=", k, "\n")
      if (iter_k==1){
        res_k = t2cd_step(list(res=res[k,], tim=tim[k,]), t.max, tau.range, deg,
                          seqby = seqby, resd.seqby = resd.seqby,
                          use_arf = use_arf, use_scale = use_scale)
        init.d[k] = res_k$d
        init.m[k] = res_k$m
        init.tau[k] = res_k$tau
        init.idx[k] = res_k$idx
      }else{
        res_k = search_tau_step(list(res=res[k,], tim=tim[k,]), t.max, tau.range, deg,
                                seqby = seqby, resd.seqby = resd.seqby,
                                use_arf = use_arf, use_scale = use_scale,
                                fix.d = d_current, fix.m = m_current[k])
        fix.d[k] = res_k$d
        fix.m[k] = res_k$m
        tau[k] = res_k$tau
        idx[k] = res_k$idx
      }
      reg2[k, (res_k$idx+1):N] = res[k, (res_k$idx+1):N]
    }

    # preprocessing
    dflag = 'original'
    if (use_scale){
      reg2 = t(scale(t(reg2), center = F))
    }
    x.2 = reg2
    if (iter_k==1){
      init.d.2 = init.d
      init.m.2 = init.m
      d_current = mean(init.d)
    }

    # optimizing
    optim_params = optim(par = d_current,
                         fn = negloglik, method = "Brent", lower = -10,
                         upper = 10)
    d_current = optim_params$par
    m_current = rep(NA, nrow(res))
    for (k in 1:length(m_current)) {
      m_current[k] <-  get.m(x.2 = na.omit(x.2[k, ]), dfrac = d_current)
    }
    neglogL_current = negloglik(optim_params$par)

    # update iter_flag
    if (iter_k==maxiter | abs(neglogL_prev-neglogL_current)<tol){
      iter_flag = FALSE
    }else{
      iter_k = iter_k + 1
    }
    # update
    hist.d = c(hist.d, d_current)
    hist.neglogL = c(hist.neglogL, neglogL_current)
    neglogL_prev = neglogL_current

  }

  return(list(res = res, tim = tim, tau.idx = tau.idx,
              tau = tau, idx = idx, d = d_current,
              m = m_current,
              univ_tau = init.tau, univ_idx = init.idx, univ_d = init.d,
              hist.d = hist.d, hist.neglogL = hist.neglogL,
              iter_k = iter_k, dflag = dflag))
}

# plot sequences and fitted lines
#' @export
plot_t2cd_step_mv = function(results, tau.range = c(10, 50), deg = 3,
                             use_arf = TRUE, use_scale = TRUE, return_plot = TRUE){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
  dflag = results$dflag
  if (use_scale){
    res_mean = t(scale(t(res), center = F)) # scaling
  }else{
    res_mean = matrix(res, ncol = 1)
  }
  N = ncol(res)

  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_idx = results$idx

  fit.vals <- var.resd1 <- var.resd1.1 <- mu <- fit1 <- matrix(nrow = nrow(res), ncol = ncol(res))

  for (k in 1:nrow(res)) {
    ### fitted values for first regime
    f1 <- refitWLS(tim[k, 1:opt_idx[k]], res_mean[k, 1:opt_idx[k]], deg = deg)
    fit1[k, 1:opt_idx[k]] = f1$fit.vals
    var.resd1.1[k, 1:opt_idx[k]] = f1$var.resd
    mu[k, ] = c(f1$fit.vals, rep(NA, N-opt_idx[k]))

    # fitted values for second regime
    if (use_arf){
      # # Has not been updated for MV setting
      # if (dflag == 'original'){
      #   x.2 = res_mean[(opt_idx+1):N]
      #   arf = arfima::arfima(x.2)
      #   fit = arf$modes[[1]]$muHat + arf$modes[[1]]$fitted
      #   mu[(opt_idx+1):N] = fit
      # }else{
      #   x.2 = res_mean[(opt_idx):N]
      #   n.2 = length(x.2) - 1
      #   arf = arfima::arfima(diff(x.2, 1))
      #   diff.2 = arf$modes[[1]]$muHat + arf$modes[[1]]$fitted
      #   mu.2 = c()
      #   for (i in 1:n.2){
      #     mu.2 = c(mu.2, x.2[i] + diff.2[i])
      #   }
      #   mu[(opt_idx+1):N] = mu.2
      # }
    }else{
      m = results$m[k]
      if (dflag == 'original'){
        x.2 = res_mean[k, (opt_idx[k]+1):N]
        diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), opt_d))
      }else{
        # # Has not been updated for MV setting
        # x.2 = c(0, diff(res_mean[k, (opt_idx+1):N]-m, 1))
        # diff_p = c(diffseries_keepmean(matrix(x.2, ncol = 1), opt_d - 1))
      }

      if (dflag == 'original'){
        mu.2 = (x.2 - m) - diff_p
        mu[k, (opt_idx[k]+1):N] = mu.2 + m
      }else{
        # # Has not been updated for MV setting
        # diff.2 = x.2 - diff_p
        # mu.2 = c()
        # for (i in (opt_idx+1):N){
        #   mu.2 = c(mu.2, res_mean[i] + diff.2[i-opt_idx])
        # }
        # mu[k, (opt_idx+1):N] = mu.2
      }
    }

    if (use_scale){
      fit.vals[k, ] = mu[k, ]*(attributes(res_mean)$'scaled:scale'[k])
      var.resd1[k, 1:length(f1$var.resd)] =
        f1$var.resd*(attributes(res_mean)$'scaled:scale'[k])^2
    }else{
      fit.vals[k, ] = mu[k, ]
      var.resd1[k, 1:length(f1$var.resd)] = f1$var.resd
    }
  }

  return(list(fit.vals = fit.vals,
              opt_idx = opt_idx, N = N,
              scaling = attributes(res_mean)$'scaled:scale',
              var.resd1 = var.resd1))
}

# parametric bootstrap using outputs from t2cd_step and plot.t2cd_step
#' @export
bootstrap_sample_step_mv = function(results, plot_results, seed = 0){

  set.seed(seed)
  res = results$res
  tim = results$tim
  N = ncol(res)
  opt_idx = results$idx

  samp <- matrix(nrow = nrow(res), ncol = ncol(res))

  for (k in 1:nrow(res)) {
    # regime 1
    fit.vals1 = plot_results$fit.vals[k, 1:opt_idx[k]]
    var.resd1 = plot_results$var.resd1[k, 1:opt_idx[k]]
    noise1 = rnorm(opt_idx[k], 0, sqrt(var.resd1))

    # regime 2
    opt_d = results$d
    m = results$m[k]
    res_mean = scale(res[k, ], center = F) # scaling
    diff_p = c(diffseries_keepmean(matrix(res_mean[(opt_idx[k]+1):N]-m, ncol = 1), opt_d))
    sd.resd2 = sqrt(mean(diff_p^2))
    sim = sim.fi(N-opt_idx[k], opt_d, sd.resd2, mu = m)
    seq_fi = sim$s

    samp[k, ] = c(fit.vals1 + noise1, (seq_fi)*plot_results$scaling[k])
  }

  return(list(res=samp, tim=tim))
}
#' @export
boot_k_mv = function(k, results, plot_results, methodname){

  if (methodname=='step'){
    samp = bootstrap_sample_step.mv(results, plot_results, seed = k)
    res = t2cd_step.mv(samp, use_arf = F)
  }else if (methodname=='sigmoid'){
    samp = bootstrap_sample_sigmoid.mv(results, plot_results, seed = k)
    res = t2cd_sigmoid.mv(samp)
  }
  tau_step = res$tau
  d_step = res$d

  return(list(tau=tau_step, d=d_step))
}
