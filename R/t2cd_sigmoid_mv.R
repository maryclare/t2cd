# multivariate implementation for T2CD-sigmoid
# optimize the likelihood for d and tau
# option to initialize tau at multiple indices
# loglikelihood, penalty to enforce tau within tau.range
#' @export
negloglik_partial_pen_res_sigmoid_mv = function(param, tim_cp, tau.idx, N, p, x.2,
                                                ll.1.mat, C, alpha0, alpha1, use_t = FALSE){
  dfrac = param

  neglogL = 0
  for (k in 1:p) {
    neglogL = neglogL -
      loglik_res_sigmoid(param = c(alpha0[k], alpha1[k], dfrac),
                         tim_cp = tim_cp[k, ], tau.idx = tau.idx,
                         x.2 = x.2[k, ], ll.1 = ll.1.mat[k, ],
                         pen = TRUE, C = C, use_t = use_t)
  }

  return(neglogL)
}

#' @export
# loglikelihood
loglik_res_sigmoid_mv = function(param, tim_cp, tau.idx, N, p, x.2,
                                 ll.1.mat, alpha0, alpha1, C = NULL, pen = FALSE, use_t = FALSE){
  dfrac = param[1]
  if (use_t) {
    m <- param[1 + 1:p]
    sd <- param[1 + p + 1:p]
    df <- param[1 + 2*p + 1:p]
  }

  logL = 0
  for (k in 1:p) {
    parev <- c(alpha0[k], alpha1[k], dfrac)
    if (use_t) {
      parev <- c(parev, m[k], sd[k], df[k])
    }
    logL = logL +
      loglik_res_sigmoid(param = parev,
                         tim_cp = tim_cp[k, ], tau.idx = tau.idx,
                         x.2 = x.2[k, ], ll.1 = ll.1.mat[k, ],
                         pen = pen, use_t = use_t, C = C)
  }

  return(logL)
}
#' @export
t2cd_sigmoid_mv = function(dat, t.max = 72, tau.range = c(10, 50),
                           init.tau = c(15, 30, 45), deg = 3, C = 1000,
                           seqby = 1, resd.seqby = 5, use_scale = TRUE, use_t = FALSE){
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
  if (!use_t) {
    init.param = matrix(NA, p, 3)
  } else {
    init.param = matrix(NA, p, 6)
  }
  init.d = rep(NA, p)
  for (k in 1:p){
    res_k = t2cd_sigmoid(list(res=matrix(res[k,], 1), tim=matrix(tim[k,], 1)),
                         t.max, tau.range, init.tau, deg, C,
                         seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale, use_t = use_t)
    init.param[k,] = res_k$par
    init.d[k] = res_k$d
  }

  x.2 = res_mean
  init.d.2 = init.d

  alpha0 = init.param[,1]
  alpha1 = init.param[,2]

  start <- c(mean(init.d.2))
  lower <- -Inf
  upper <- Inf
  if (use_t) {
    start <- c(start, init.param[, 4], init.param[, 5], init.param[, 6])
    lower <- c(lower, rep(-Inf, 2*p), rep(2 + 10^(-14), p))
    upper <- c(upper, rep(Inf, 2*p), rep(300, p))
  }

  optim_params = optim(par = start,
                       fn = loglik_res_sigmoid_mv, method = "L-BFGS-B",
                       lower = lower, upper = upper,
                       tim_cp = tim_cp, tau.idx = tau.idx, N = N, p = p, x.2 = x.2,
                       ll.1.mat = ll.1.mat, C = C, alpha0 = alpha0, alpha1 = alpha1, pen = TRUE,
                       control = list("fnscale" = -1), use_t = use_t)
  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(ifelse(wt_cp[,ncol(wt_cp)] != 0, wt_cp[,ncol(wt_cp)], 1),
                        N-tau.idx[length(tau.idx)]), p, byrow = F))

  opt_param = c(alpha0, alpha1, optim_params$par)
  opt_logL = loglik_res_sigmoid_mv(optim_params$par,
                                   tim_cp = tim_cp, tau.idx = tau.idx, N = N,
                                   p = p, x.2 = x.2,
                                   ll.1.mat = ll.1.mat, alpha0 = alpha0, alpha1 = alpha1, use_t = use_t)
  opt_d = optim_params$par[1]
  univ_d = init.param[,3]

  df <- sd <- m <- rep(NA, nrow(x.2))

  for (k in 1:length(m)) {
    if (!use_t) {
      m[k] <- get.m(x.2 = x.2[k, ], dfrac = opt_d, wt = wt[k, ])
      sd[k] <- get.sd(x.2 = x.2[k, ], dfrac = opt_d, mean = m[k], wt = wt[, ])
    } else {
      m[k] <- optim_params$par[1 + k]
      sd[k] <- abs(optim_params$par[1 + p + k])
      df[k] <- optim_params$par[1 + 2*p + k]
    }
  }


  opt_tau.idx = apply(wt, 1, function(x){return(which(x>=0.5, arr.ind = TRUE)[1]-1)})
  opt_tau = c()
  for (k in 1:p){
    opt_tau = c(opt_tau, tim[k,opt_tau.idx[k]])
  }

  return(list(res = res, tim = tim, tau.idx = tau.idx, idx = opt_tau.idx,
              d = opt_d, univ_d = univ_d, tau = opt_tau, param = opt_param, logL = opt_logL,
              m = m, sd = sd, df = df))
}

# plot sequences and fitted lines
#' @export
plot_t2cd_sigmoid_mv = function(results, tau.range = c(10, 50), deg = 3,
                                seqby = 1, resd.seqby = 5, return_plot = TRUE){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
  res_mean = t(scale(t(res), center = F)) # scaling
  N = ncol(res_mean)
  p = nrow(res_mean)

  fit1 = var.resd1 = matrix(nrow = nrow(res), ncol = ncol(res))
  for (k in 1:p){
    fitwls = refitWLS(tim[k,], res_mean[k,], deg = deg,
                      seqby = seqby, resd.seqby = resd.seqby)
    fit1[k, ] <- fitwls$fit.vals
    var.resd1[k, ] <- fitwls$var.resd
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
  for (k in 1:p) {
    wt_cp[k, ] = sigmoid(alpha0[k]+alpha1[k]*tim_cp[k, ])
  }
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(ifelse(wt_cp[,ncol(wt_cp)] != 0, wt_cp[,ncol(wt_cp)], 1),
                        N-tau.idx[length(tau.idx)]), p, byrow = F))

  # update variables if using original or first difference
  fit.vals <- mu.2 <- diff_p <- matrix(nrow = nrow(res), ncol = ncol(res))
  for (k in 1:p) {
    x.2 = res_mean[k, ]
    d = opt_d
    diff_p[k, ] = t(diffseries_keepmean(t(wt[k, ]*(x.2-m[k])), d))

    mu.2[k, ] = wt[k, ]*(x.2-m[k]) - diff_p[k, ]
    fit.vals[k, ] = (mu.2[k, ] + wt[k, ]*m[k])*(attributes(res_mean)$'scaled:scale'[k]) +
    (1-wt[k, ])*fit1[k, ]*(attributes(res_mean)$'scaled:scale'[k])

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
bootstrap_sample_sigmoid_mv = function(results, plot_results, seed = 0){

  set.seed(seed)
  res = results$res
  tim = results$tim
  N <- ncol(res)
  p <- nrow(res)

  opt_d = results$d

  samp <- matrix(nrow = p, ncol = N)

  for (k in 1:p) {
    # regime 1
    fit.vals1 = plot_results$fit.vals[k, ]
    var.resd1 = plot_results$var.resd1[k, ]
    noise1 = rnorm(N, 0, sqrt(var.resd1))

    # regime 2

    m = results$m[k]
    wt = plot_results$wt[k, ]
    x.2 = t(scale(t(res), center = F))[k, ] # scaling
    diff_p = t(diffseries_keepmean(t(wt*(x.2-m)), opt_d))
    sd.resd2 = sqrt(sum(wt*diff_p^2)/sum(wt))
    sim = sim.fi(N, opt_d, sd.resd2, mu = m)
    seq_fi = sim$s

    samp[k, ] = c((1-wt)*(fit.vals1 + noise1) + wt*(seq_fi)*plot_results$scaling[k])
  }

  return(list(res=samp, tim=tim))
}
