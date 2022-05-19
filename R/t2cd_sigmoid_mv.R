# multivariate implementation for T2CD-sigmoid
# optimize the likelihood for d and tau
# option to initialize tau at multiple indices
# loglikelihood, penalty to enforce tau within tau.range
#' @export
# loglikelihood
loglik_res_sigmoid_mv = function(param, tim_cp, tau.idx, p, x,
                                 ll.1.mat, alpha0, alpha1,
                                 C = NULL, pen = FALSE, use_t = FALSE){
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
      loglik_res(param = parev, sigmoid = TRUE,
                 tim_cp = na.omit(tim_cp[k, ]), tau.idx = na.omit(tau.idx[k, ]),
                 x = na.omit(x[k, ]), ll.1 = na.omit(ll.1.mat[k, ]),
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

  max.tim.ind =
    max(apply(dat$tim, 1, function(x) {max(which(!is.na(x) & x <= t.max))}))

  res = matrix(dat$res[, 1:max.tim.ind], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, 1:max.tim.ind], nrow = nrow(dat$tim))

  res[!is.na(tim) & tim > t.max] <- NA
  tim[!is.na(tim) & tim > t.max] <- NA
  if (use_scale){
    res_mean = t(scale(t(res), center = F)) # scaling
  }else{
    res_mean = res
  }
  N = ncol(res)
  p = nrow(res)
  x = res_mean

  # points in change range
  tim_cp = tim
  tim_cp[!is.na(tim) & (tim < tau.range[1] | tim > tau.range[2])] <- NA

  min.tau.idx <- apply(tim_cp, 1, function(x) {min(which(!is.na(x)))})
  max.tau.idx <- apply(tim_cp, 1, function(x) {max(which(!is.na(x)))})

  tau.idx <- matrix(NA, nrow = nrow(tim), ncol = max(max.tau.idx - min.tau.idx) + 1)
  for (k in 1:p) {
    tau.idx[k, 1:(max.tau.idx[k] - min.tau.idx[k] + 1)] <- min.tau.idx[k]:max.tau.idx[k]
  }

  tim_cp = tim_cp[, !apply(tim_cp, 2, function(x) {sum(is.na(x))}) == nrow(tim_cp), drop = FALSE]

  # determine range of d to search and find initializing parameters
  if (!use_t) {
    init.param = matrix(NA, p, 3)
  } else {
    init.param = matrix(NA, p, 6)
  }
  init.d = rep(NA, p)
  fit.vals <- var.resd <- ll.1.mat <- matrix(NA, ncol = N, nrow = p)
  for (k in 1:p){
    res_k = t2cd_sigmoid(list(res=matrix(x[k,], 1), tim=matrix(tim[k,], 1)),
                         t.max, tau.range, init.tau, deg, C,
                         seqby = seqby, resd.seqby = resd.seqby,
                         use_scale = FALSE, use_t = use_t)
    init.param[k,] = res_k$par
    init.d[k] = res_k$d
    ll.1.mat[k, 1:length(res_k$ll.1)] <- res_k$ll.1
    fit.vals[k, 1:length(res_k$fit.vals)] <- res_k$fit.vals
    var.resd[k, 1:length(res_k$var.resd)] <- res_k$var.resd
  }

  alpha0 = init.param[,1]
  alpha1 = init.param[,2]
  univ_d = init.param[,3]
  opt_d = mean(univ_d)

  if (p > 1) {
  start <- c(mean(init.d))
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
                       tim_cp = tim_cp, tau.idx = tau.idx, p = p, x = x,
                       ll.1.mat = ll.1.mat, C = C,
                       alpha0 = alpha0, alpha1 = alpha1, pen = TRUE,
                       control = list("fnscale" = -1), use_t = use_t)
  opt_param = c(alpha0, alpha1, optim_params$par)
  opt_logL = loglik_res_sigmoid_mv(optim_params$par,
                                   tim_cp = tim_cp, tau.idx = tau.idx,
                                   p = p, x = x,
                                   ll.1.mat = ll.1.mat, alpha0 = alpha0,
                                   alpha1 = alpha1, use_t = use_t)
  opt_d = optim_params$par[1]
  } else {
    opt_param = res_k$param
    opt_logL = res_k$logL
  }
  # weights
  wt <- matrix(NA, nrow = p, ncol = N)
  for (k in 1:p) {
    wt[k, 1:length(na.omit(x[k, ]))] <-
      get.wt(alpha = c(alpha0[k], alpha1[k]),
             tim_cp = na.omit(tim_cp[k, ]),
             tau.idx = na.omit(tau.idx[k, ]),
             N = length(na.omit(x[k, ])))
  }

  df <- sd <- m <- rep(NA, nrow(x))

  for (k in 1:length(m)) {
    if (!use_t) {
      m[k] <- get.m(x.2 = na.omit(x[k, ]), dfrac = opt_d, wt = na.omit(wt[k, ]))
      sd[k] <- get.sd(x.2 = na.omit(x[k, ]), dfrac = opt_d, mean = m[k],
                      wt = na.omit(wt[k, ]))
    } else {
      m[k] <- opt_param[2*p + 1 + k]
      sd[k] <- abs(opt_param[2*p + 1 + p + k])
      df[k] <- opt_param[2*p + 1 + 2*p + k]
    }
  }


  opt_tau.idx = apply(wt, 1, function(x){return(which(x>=0.5, arr.ind = TRUE)[1]-1)})
  opt_tau = c()
  for (k in 1:p){
    opt_tau = c(opt_tau, tim[k,opt_tau.idx[k]])
  }

  return(list(res = res, res_mean = res_mean, tim = tim,
              tim_cp = tim_cp,
              tau.idx = tau.idx,
              idx = opt_tau.idx,
              d = opt_d, univ_d = univ_d, tau = opt_tau, param = opt_param,
              logL = opt_logL,
              m = m, sd = sd, df = df,
              fit.vals = fit.vals, var.resd = var.resd))
}

# plot sequences and fitted lines
#' @export
plot_t2cd_sigmoid_mv = function(results, tau.range = c(10, 50),
                                return_plot = TRUE, use_scale = TRUE){
  res = results$res
  tim = results$tim
  tim_cp = results$tim_cp
  tau.idx = results$tau.idx
  res_mean = results$res_mean # scaling
  N = ncol(res_mean)
  p = nrow(res_mean)

  fit1 <- results$fit.vals
  var.resd1 <- results$var.resd

  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_param = results$param

  # fitted values
  alpha0 = opt_param[1:p]
  alpha1 = opt_param[p + 1:p]
  m = results$m

  # weights
  wt <- matrix(nrow = p, ncol = N)
  for (k in 1:p) {
    N.k <- length(na.omit(res[k, ]))
    wt[k, 1:N.k] <-
      get.wt(alpha = c(alpha0[k], alpha1[k]),
             tim_cp = na.omit(tim_cp[k, ]),
             tau.idx = na.omit(tau.idx[k, ]),
             N = N.k)
  }

  # update variables if using original or first difference
  fit.vals <- mu.2 <- diff_p <- matrix(nrow = nrow(res), ncol = ncol(res))
  for (k in 1:p) {
    N.k <- length(na.omit(res[k, ]))
    x.2 = res_mean[k, ]
    d = opt_d
    diff_p[k, 1:N.k] = t(diffseries_keepmean(t(na.omit(wt[k, ]*(x.2-m[k]))), d))

    mu.2[k, ] = wt[k, ]*(x.2-m[k]) - diff_p[k, ]
    fit.vals[k, ] = (mu.2[k, ] + wt[k, ]*m[k]) +(1-wt[k, ])*fit1[k, ]
    if (use_scale) {
      fit.vals[k, ] = fit.vals[k, ]*(attributes(res_mean)$'scaled:scale'[k])
      var.resd1[k, ] = var.resd1[k, ]*attributes(res_mean)$'scaled:scale'[k]^2
    }

  }

  # plotting
  if (return_plot){
    plot(tim[1, ], res[1, ],
         ylim = c(min(c(res[1, ], fit.vals[1, ])), max(c(res[1, ], fit.vals[1, ]))), type = 'l',
         main = paste('Values fitted with d: ', round(opt_d,3), ' tau: ', round(mean(opt_tau),3)),
         xlab = 'Time (hour)', ylab = 'Resistance (ohm)')

    if (is.na(opt_tau[1])){
      lines(tim[1, ], fit.vals[1, ], col = "blue", lwd = 1)
    }else{
      opt_tau.idx = which(tim == opt_tau[1])
      lines(tim[1, 1:opt_tau.idx],
            fit.vals[1, 1:opt_tau.idx], col = "blue", lwd = 1)
      lines(tim[1, (opt_tau.idx):ncol(fit.vals)],
            fit.vals[1, (opt_tau.idx):ncol(fit.vals)], col = "green", lwd = 1)
      abline(v = opt_tau[1], lty = 2, col = "red")
    }
    abline(v = tau.range, lty = 1, col = "red")
  }

  return(list(fit.vals = fit.vals,
                var.resd1 = var.resd1,
                wt = wt, scaling = attributes(res_mean)$'scaled:scale'))
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
    sd = results$sd[k]
    df = results$df[k]
    wt = plot_results$wt[k, ]
    sim = sim_fi(N, opt_d, mu = m, sig = sd, df = df)
    seq_fi = sim$s

    samp[k, ] = c((1-wt)*(fit.vals1 + noise1) + wt*(seq_fi)*plot_results$scaling[k])
  }

  return(list(res=samp, tim=tim))
}
