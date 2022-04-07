# grid search for tau given FI parameters
# return the likelihood for each combination
#' @export
search_tau_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                           seqby = 1, resd.seqby = 5,
                           use_scale = TRUE,
                           fix.d = NULL, fix.m = NULL, fix.sd = NULL, fix.df = NA, use_t = FALSE){

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
  M <- c()
  d <- c()
  m <- c()
  sd <- c()
  df <- c()

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

    x.2 = res_mean[(tau_j+1):N]
    # MCG Change: Make sure d and m aren't re-estimated if they are provided
    d <- c(d, fix.d)
    m <- c(m, fix.m)
    sd <- c(sd, fix.sd)
    df <- c(df, fix.df)
    if (!use_t) {
      ll.2 <- loglik_res_step(fix.d, x.2 = x.2)
    } else {
      ll.2 <- loglik_t_res_step(c(fix.d, fix.m, fix.sd, fix.df), x.2 = x.2)
    }

    M = c(M, ll.1 + ll.2)
  }

  # tau and d at maximum log-likelihood
  M_df = data.frame(tau = tim[tau.idx], M = M, d = d, m = m, sd = sd, df = df)
  max.idx = which.max(M)
  max.tau = tim[tau.idx[max.idx]]
  max.d = d[max.idx]
  max.m = m[max.idx]
  max.sd = sd[max.idx]
  max.df = df[max.idx]

  return(list(M_df = M_df, res = res, tim = tim, tau.idx = tau.idx,
              tau = max.tau, d = max.d, m = max.m,  sd = max.sd, df = max.df,
              idx = tau.idx[max.idx], logL = M[max.idx]))

}

# helper functions for loglikelihood
# MCG Change: Compute profile likelihood, profiling out means
#' @export
negloglik_res_step_mv = function(param, x.2, use_t = FALSE){
  dfrac = param[1]
  if (use_t) {
    p <- nrow(x.2)
    mean <- param[1 + 1:p]
    sd <- param[1 + p + 1:p]
    df <- param[1 + 2*p + 1:p]
  }

  neglogL <- 0
  for (k in 1:nrow(x.2)) {
    if (!use_t) {
      nll <- -loglik_res_step(par = dfrac, x.2 = na.omit(x.2[k, ]))
    } else {
      nll <- -loglik_t_res_step(par = c(dfrac, mean[k], sd[k], df[k]), x.2 = na.omit(x.2[k, ]))
    }
    neglogL <- neglogL + nll
  }
  # n.2 = rowSums(!is.na(x.2))
  # m <- rep(NA, nrow(x.2))
  # for (k in 1:length(m)) {
  #   m[k] <- get.m(x.2 = na.omit(x.2[k, ]), dfrac = dfrac)
  # }
  # diff_p = t(diffseries_keepmean(t(x.2-m), dfrac))
  #
  # diff_p[diff_p==0] = NA
  # neglogL = sum(n.2*log(rowMeans(diff_p^2, na.rm = T)))

  return(neglogL)
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
  init.d = init.m = init.sd = init.df =
    fix.d = fix.m = fix.sd = fix.df =
    init.tau = init.idx = tau = idx = rep(NA, p)
  hist.d = hist.neglogL = c()

  while (iter_flag){
    # retain all data in regime 2
    reg2 = matrix(NA, nrow = p, ncol = N)
    for (k in 1:p){
      cat("k=", k, "\n")
      if (iter_k==1){
        res_k = t2cd_step(list(res=res[k,], tim=tim[k,]), t.max, tau.range, deg,
                          seqby = seqby, resd.seqby = resd.seqby,
                          use_scale = use_scale, use_t = use_t)
        init.d[k] = res_k$d
        init.m[k] = res_k$m
        init.sd[k] = res_k$sd
        init.df[k] = res_k$df
        init.tau[k] = res_k$tau
        init.idx[k] = res_k$idx
      }else{
        res_k = search_tau_step(list(res=res[k,], tim=tim[k,]), t.max, tau.range, deg,
                                seqby = seqby, resd.seqby = resd.seqby,
                                use_scale = use_scale,
                                fix.d = d_current, fix.m = m_current[k],
                                fix.sd = sd_current[k], fix.df = df_current[k],
                                use_t = use_t)
        fix.d[k] = res_k$d
        fix.m[k] = res_k$m
        fix.sd[k] = res_k$sd
        fix.df[k] = res_k$df
        tau[k] = res_k$tau
        idx[k] = res_k$idx
      }
      reg2[k, (res_k$idx+1):N] = res[k, (res_k$idx+1):N]
    }
    print(res_k$idx)

    # preprocessing
    if (use_scale){
      reg2 = t(scale(t(reg2), center = F))
    }
    x.2 = reg2
    if (iter_k==1){
      d_current = mean(init.d)
      m_current = init.m
      sd_current = init.sd
      df_current = init.df
    }

    # optimizing over d, m (fixed taus, I think)
    if (!use_t) {
    optim_params = optim(par = d_current,
                         fn = negloglik_res_step_mv, method = "Brent", lower = -10,
                         upper = 10, x.2 = x.2)
    d_current = optim_params$par
    df_current <- sd_current <- m_current <- rep(NA, nrow(res))
    for (l in 1:length(m_current)) {
      m_current[l] <-  get.m(x.2 = na.omit(x.2[l, ]), dfrac = d_current)
      sd_current[l] <-  get.sd(x.2 = na.omit(x.2[l, ]), dfrac = d_current, mean = m_current[l])
    }
    neglogL_current = optim_params$value
    } else {
      optim_params = optim(par = c(d_current, m_current, sd_current, df_current),
                           fn = negloglik_res_step_mv, method = "L-BFGS-B",
                           lower = c(-Inf, rep(c(-Inf, -Inf, 2 + 10^(-14)), each = length(m_current))),
                           upper = c(Inf, rep(c(Inf, Inf, 300), each = length(m_current))),
                           x.2 = x.2, use_t = use_t)
      d_current = optim_params$par[1]
      m_current = optim_params$par[1 + 1:length(m_current)]
      sd_current = optim_params$par[(1 + length(m_current)) + 1:length(m_current)]
      df_current = optim_params$par[(1 + 2*length(m_current)) + 1:length(m_current)]
      neglogL_current = optim_params$value
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

  }

  return(list(res = res, tim = tim, tau.idx = tau.idx,
              tau = tau, idx = idx, d = d_current,
              m = m_current, sd = sd_current, df = df_current,
              univ_tau = init.tau, univ_idx = init.idx, univ_d = init.d,
              hist.d = hist.d, hist.neglogL = hist.neglogL,
              iter_k = iter_k))
}

# plot sequences and fitted lines
#' @export
plot_t2cd_step_mv = function(results, tau.range = c(10, 50), deg = 3,
                             use_scale = TRUE, return_plot = TRUE){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
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
    m = results$m[k]
    x.2 = res_mean[k, (opt_idx[k]+1):N]
    diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), opt_d))

    mu.2 = (x.2 - m) - diff_p
    mu[k, (opt_idx[k]+1):N] = mu.2 + m

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
    res = t2cd_step.mv(samp)
  }else if (methodname=='sigmoid'){
    samp = bootstrap_sample_sigmoid.mv(results, plot_results, seed = k)
    res = t2cd_sigmoid.mv(samp)
  }
  tau_step = res$tau
  d_step = res$d

  return(list(tau=tau_step, d=d_step))
}
