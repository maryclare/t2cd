
# multivariate implementation for T2CD-step
# grid search for d and tau
t2cd_step.mv = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
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
    }else{
      init.d.2 = fix.d
      init.m.2 = fix.m
    }

    # helper functions for loglikelihood
    negloglik = function(param){
      m = param[1:p]
      dfrac = param[p+1]

      n.2 = rowSums(!is.na(x.2))
      diff_p = t(diffseries_keepmean(t(x.2-m), dfrac))

      diff_p[diff_p==0] = NA
      neglogL = sum(n.2*log(rowMeans(diff_p^2, na.rm = T)))

      return(neglogL)
    }

    # optimizing
    optim_params = optim(par = c(init.m.2, mean(init.d.2)),
                         fn = negloglik, method = "BFGS")
    d_current = optim_params$par[p + 1]
    m_current = optim_params$par[1:p]
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


# grid search for tau given FI parameters
# return the likelihood for each combination
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
    # helper functions for loglikelihood
    negloglik = function(param){
      m = param[1]
      dfrac = param[2]
      if (dfrac<=-0.5 | dfrac>=1.5){
        return(1e+10)
      }

      diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), dfrac))

      neglogL = sum(diff_p^2)

      return(neglogL)
    }
    loglik = function(param){
      m = param[1]
      dfrac = param[2]

      diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), dfrac))

      logL = sum(dnorm(diff_p, log = TRUE, sd = sd(diff_p)))

      return(logL)
    }

    # optimizing
    optim.2 = optim(par = c(fix.m, fix.d),
                    fn = negloglik, method = "BFGS")
    if (dflag == 'original'){
      d = c(d, optim.2$par[2])
    }else{
      d = c(d, optim.2$par[2] + 1)
    }
    m = c(m, optim.2$par[1])
    ll.2 = loglik(optim.2$par)

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

# plot sequences and fitted lines
plot.t2cd_step.mv = function(results, tau.range = c(10, 50), deg = 3,
                          use_arf = TRUE, use_scale = TRUE, return_plot = TRUE){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
  dflag = results$dflag
  if (use_scale){
    res_sd = c(apply(res, 1, function(r) {sqrt(sum(r^2)/(length(r) - 1))}))
    res_mean = t(apply(res, 1, function(r) {scale(r, center = F)})) # scaling
  }else{
    res_mean = matrix(c(res), nrow = nrow(res), ncol = ncol(res))
  }
  N = ncol(res)

  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_idx = results$idx

  ### fitted values for first regime
  fit.vals <- mu <- var.resd1 <- var.resd1.1 <- matrix(NA, nrow = nrow(res), ncol = ncol(res))
  for (i in 1:nrow(res)) {
    fit <- refitWLS(tim[i, 1:opt_idx[i]],
             res_mean[i, 1:opt_idx[i]], deg = deg)
    mu[i, 1:opt_idx[i]] <- fit$fit.vals
    var.resd1.1[i, 1:opt_idx[i]] <- fit$var.resd

    # fitted values for second regime
    if (use_arf){
      if (dflag == 'original'){
        x.2 = res_mean[i, (opt_idx[i]+1):N]
        arf = arfima::arfima(x.2)
        fit = arf$modes[[1]]$muHat + arf$modes[[1]]$fitted
        mu[i, (opt_idx[i]+1):N] = fit
      }else{
        x.2 = res_mean[i, (opt_idx[i]):N]
        n.2 = length(x.2) - 1
        arf = arfima::arfima(diff(x.2, 1))
        diff.2 = arf$modes[[1]]$muHat + arf$modes[[1]]$fitted
        mu.2 = c()
        for (i in 1:n.2){
          mu.2 = c(mu.2, x.2[i] + diff.2[i])
        }
        mu[i, (opt_idx[i]+1):N] = mu.2
      }
    } else{
      m = results$m[i]
      if (dflag == 'original'){
        x.2 = res_mean[i, (opt_idx[i]+1):N]
        diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), opt_d))
      } else{
        x.2 = c(0, diff(res_mean[i, (opt_idx[i]+1):N]-m, 1))
        diff_p = c(diffseries_keepmean(matrix(x.2, ncol = 1), opt_d - 1))
      }

      if (dflag == 'original'){
        mu.2 = (x.2 - m) - diff_p
        mu[i, (opt_idx[i]+1):N] = mu.2 + m
      }else{
        diff.2 = x.2 - diff_p
        mu.2 = c()
        for (i in (opt_idx+1):N){
          mu.2 = c(mu.2, res_mean[i] + diff.2[i-opt_idx])
        }
        mu[i, (opt_idx[i]+1):N] = mu.2
      }
    }

    if (use_scale){
      fit.vals[i, ] = mu[i, ]*res_sd[i]
      var.resd1[i, ] = var.resd1.1[i, ]*res_sd[i]^2
    }else{
      fit.vals[i, ] = mu[i, ]
      var.resd1[i, ] = var.resd1.1[i, ]
    }

  }

  return(list(fit.vals1 = fit.vals, fit.vals2 = fit.vals,
              opt_idx = opt_idx, N = N,
              var.resd1 = var.resd1))
}

# parametric bootstrap using outputs from t2cd_step and plot.t2cd_step
bootstrap_sample.mv = function(results, plot_results, seed = 0){

  set.seed(seed)
  res = results$res
  tim = results$tim
  N = ncol(res)
  opt_idx = results$idx

  samps <- matrix(nrow = nrow(res), ncol = ncol(res))

  for (i in 1:nrow(res)) {
   # regime 1
   fit.vals1 = plot_results$fit.vals1[i, 1:opt_idx[i]]
   var.resd1 = plot_results$var.resd1[i, 1:opt_idx[i]]
   noise1 = rnorm(opt_idx[i], 0, sqrt(var.resd1))

   # regime 2
   opt_d = results$d
   fit.vals2 = plot_results$fit.vals2[i, (opt_idx[i]+1):N]
   sd.resd2 = sd(res[(opt_idx[i]+1):N]-fit.vals2)
   sim = sim.fi(N-opt_idx[i], opt_d, sd.resd2)
   seq_fi = sim$s

   samps[i, ] = c(fit.vals1 + noise1, seq_fi + fit.vals1[opt_idx[i]])
  }

  return(list(res=samps, tim=tim))
}
