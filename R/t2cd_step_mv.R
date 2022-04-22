# helper functions for loglikelihood
# MCG Change: Compute profile likelihood, profiling out means
#' @export
loglik_res_step_mv = function(param, x, idx, ll.1,
                              use_t = FALSE){
  dfrac = param[1]
  if (use_t) {
    p <- nrow(x)
    mean <- param[1 + 1:p]
    sd <- param[1 + p + 1:p]
    df <- param[1 + 2*p + 1:p]
  }

  logL <- 0
  for (k in 1:nrow(x)) {
    if (!use_t) {
      par <- dfrac
    } else {
      par <- c(dfrac, mean[k], sd[k], df[k])
    }
    ll <- loglik_res(par = par, x = x[k, ],
                     sigmoid = FALSE, ll.1 = ll.1[k, ],
                     use_t = use_t, idx = idx[k])

    logL <- logL + ll
  }

  return(logL)
}

# multivariate implementation for T2CD-step
# grid search for d and tau
#' @export
t2cd_step_mv = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                        seqby = 1, resd.seqby = 5,
                        use_scale = TRUE,
                        maxiter=10, tol=1e-6, use_t = FALSE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # deg: degree for B-spline
  # seqby, resd.seqby: interval between knots
  # use_scale: if true, scale time series
  # maxiter: maximum number of iterations solving for tau and FI parameters
  # tol: stops iterating if difference between the most recent objective values is below tol

  if (!is.matrix(dat$res)) {
    dat$res <- matrix(dat$res, nrow = 1, ncol = length(dat$res))
  }
  if (!is.matrix(dat$tim)) {
    dat$tim <- matrix(dat$tim, nrow = 1, ncol = length(dat$tim))
  }

  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }

  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)

  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  if (use_scale){
    res_mean = t(scale(t(res), center = F)) # scaling
  }else{
    if (!is.matrix(res)) {
      res_mean = matrix(res, ncol = 1)
    } else {
      res_mean = res
    }
  }
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  N = ncol(res_mean)
  p = nrow(res_mean)

  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))

  # alternate between optimizing for tau and FI parameters
  iter_k = 1
  iter_flag = TRUE
  neglogL_prev = Inf
  init.d = init.m = init.sd = init.df =
    fix.d = fix.m = fix.sd = fix.df =
    init.tau = init.idx = tau = idx = rep(NA, p)
  hist.d = hist.neglogL = c()

  while (iter_flag){
    # retain all data in regime 2
    reg2 = matrix(NA, nrow = p, ncol = N)
    fit.vals <- var.resd <- ll.1 <- matrix(0, nrow = p, ncol = N)
    for (k in 1:p){
      cat("k=", k, "\n")
      if (iter_k==1){
        res_k = t2cd_step(list(res=res_mean[k,], tim=tim[k,]),
                          t.max, tau.range, deg,
                          seqby = seqby, resd.seqby = resd.seqby,
                          use_scale = FALSE, use_t = use_t)
        init.d[k] = res_k$d
        init.m[k] = res_k$m
        init.sd[k] = res_k$sd
        init.df[k] = res_k$df
        init.tau[k] = res_k$tau
        init.idx[k] = res_k$idx
        ll.1[k, 1:length(res_k$ll.1)] <- res_k$ll.1
        fit.vals[k, 1:length(res_k$fit.vals)] <- res_k$fit.vals
        var.resd[k, 1:length(res_k$var.resd)] <- res_k$var.resd
      }else{
        res_k = search_dtau_step(list(res=res_mean[k,], tim=tim[k,]), t.max,
                                tau.range, deg,
                                seqby = seqby, resd.seqby = resd.seqby,
                                use_scale = FALSE,
                                fix.d = d_current, fix.m = m_current[k],
                                fix.sd = sd_current[k], fix.df = df_current[k],
                                use_t = use_t, cp.only = TRUE)
        fix.d[k] = res_k$d
        fix.m[k] = res_k$m
        fix.sd[k] = res_k$sd
        fix.df[k] = res_k$df
        tau[k] = res_k$tau
        idx[k] = res_k$idx
        ll.1[k, 1:length(res_k$ll.1)] <- res_k$ll.1
        fit.vals[k, 1:length(res_k$fit.vals)] <- res_k$fit.vals
        var.resd[k, 1:length(res_k$var.resd)] <- res_k$var.resd
      }
      reg2[k, (res_k$idx+1):N] = res_mean[k, (res_k$idx+1):N]
    }

    if (iter_k==1){
      d_current = mean(init.d)
      m_current = init.m
      sd_current = init.sd
      df_current = init.df
      tau = init.tau
      idx = init.idx
    }

    if (p > 1) {
      x = reg2
    # optimizing over d, m (fixed taus, I think)
    if (!use_t) {
    optim_params = optim(par = d_current,
                         fn = loglik_res_step_mv, method = "Brent", lower = -10,
                         upper = 10, x = x, idx = idx, ll.1 = ll.1,
                         use_t = use_t, control = list("fnscale" = -1))
    d_current = optim_params$par
    df_current <- sd_current <- m_current <- rep(NA, nrow(res))
    for (l in 1:length(m_current)) {
      m_current[l] <-  get.m(x.2 = x[l, (idx[l] + 1):N], dfrac = d_current)
      sd_current[l] <-  get.sd(x.2 = x[l, (idx[l] + 1):N], dfrac = d_current,
                               mean = m_current[l])
    }
    neglogL_current = -optim_params$value
    } else {
      optim_params = optim(par = c(d_current, m_current, sd_current, df_current),
                           fn = loglik_res_step_mv, method = "L-BFGS-B",
                           lower = c(-Inf, rep(c(-Inf, -Inf, 2 + 10^(-14)), each = length(m_current))),
                           upper = c(Inf, rep(c(Inf, Inf, 300), each = length(m_current))),
                           x = x, idx = idx, ll.1 = ll.1,
                           use_t = use_t, control = list("fnscale" = -1))
      d_current = optim_params$par[1]
      m_current = optim_params$par[1 + 1:length(m_current)]
      sd_current = optim_params$par[(1 + length(m_current)) + 1:length(m_current)]
      df_current = optim_params$par[(1 + 2*length(m_current)) + 1:length(m_current)]
      neglogL_current = -optim_params$value
    }

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
    } else {
      iter_flag = FALSE
    }

  }

  return(list(res = res, res_mean = res_mean, tim = tim, tau.idx = tau.idx,
              tau = tau, idx = idx, d = d_current,
              m = m_current, sd = sd_current, df = df_current,
              univ_tau = init.tau, univ_idx = init.idx, univ_d = init.d,
              hist.d = hist.d, hist.neglogL = hist.neglogL,
              iter_k = iter_k,
              fit.vals = fit.vals, var.resd = var.resd))
}

# plot sequences and fitted lines
#' @export
plot_t2cd_step_mv = function(results, tau.range = c(10, 50),
                             use_scale = TRUE, return_plot = FALSE) {
  res = results$res
  res_mean = results$res_mean
  tim = results$tim
  tau.idx = results$tau.idx
  N = ncol(res)

  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_idx = results$idx
  m <- results$m

  fit.vals <- var.resd1 <- var.resd1.1 <- mu <- fit1 <-
    matrix(nrow = nrow(res_mean), ncol = ncol(res_mean))

  for (k in 1:nrow(res)) {
    ### fitted values for first regime
    fit1[k, 1:opt_idx[k]] = res_step$fit.vals[k, 1:opt_idx[k]]
    var.resd1.1[k, 1:opt_idx[k]] = res_step$var.resd[k, 1:opt_idx[k]]
    x.2 = res_mean[k, (opt_idx[k]+1):N]
    mu[k, ] = c(fit1[k, 1:opt_idx[k]],
                trend_fi(s = x.2, eff.d = opt_d,  mu = m[k]))

    # Back on original scale
    if (use_scale){
      fit.vals[k, ] = mu[k, ]*(attributes(res_mean)$'scaled:scale'[k])
      var.resd1[k, 1:opt_idx[k]] =
        var.resd1.1[k, 1:opt_idx[k]]*(attributes(res_mean)$'scaled:scale'[k])^2
    }else{
      fit.vals[k, ] = mu[k, ]
      var.resd1[k, 1:opt_idx[k]] = var.resd1.1[k, 1:opt_idx[k]]
    }
  }

  # plotting
  if (return_plot){
    plot(tim[1, ], res[1, ],
         ylim = c(min(c(res[1, ], fit.vals[1, ])),
                  max(c(res[1, ], fit.vals[1, ]))), type = 'l',
         main = paste('Values fitted with d: ', round(opt_d,3), ' tau: ', round(opt_tau,3)),
         xlab = 'Time (hour)', ylab = 'Resistance (ohm)')
    if (is.na(opt_tau)){
      lines(tim[1, ], fit.vals[1, ], col = "blue", lwd = 1)
    }else{
      opt_tau.idx = which(tim[1, ] == opt_tau)
      lines(tim[1, 1:opt_idx], fit.vals[1, 1:opt_idx], col = "blue", lwd = 1)
      lines(tim[1, (opt_idx+1):N], fit.vals[1, (opt_idx+1):N], col = "green", lwd = 1)
      abline(v = opt_tau, lty = 2, col = "red")
    }
    abline(v = tau.range, lty = 1, col = "red")
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
    sd = results$sd[k]
    df <- results$df[k]
    sim = sim_fi(N-opt_idx[k], opt_d, mu = m, sig = sd, df = df)
    seq_fi = sim$s

    samp[k, ] = c(fit.vals1 + noise1, seq_fi*plot_results$scaling[k])
  }

  return(list(res=samp, tim=tim))
}
#' @export
boot_k_mv = function(k, results, plot_results, methodname){

  if (methodname=='step'){
    samp = bootstrap_sample_step_mv(results, plot_results, seed = k)
    use_t <- ifelse(any(is.na(results$df)), FALSE, TRUE)
    res = t2cd_step_mv(samp, use_t = use_t)
  }else if (methodname=='sigmoid'){
    samp = bootstrap_sample_sigmoid_mv(results, plot_results, seed = k)
    use_t <- ifelse(any(is.na(results$df)), FALSE, TRUE)
    res = t2cd_sigmoid_mv(samp, use_t = use_t)
  }
  tau = res$tau
  d = res$d

  return(list(tau=tau, d=d))
}
