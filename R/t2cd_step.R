# helper functions for loglikelihood
#' @export
get.m <- function(x.2, dfrac, wt = rep(1, length(x.2))) {
  wt <- c(wt)
  diff.wt <- diffseries_keepmean(wt, dfrac)
  diff.wx <- diffseries_keepmean(wt*x.2, dfrac)
  sum(wt*diff.wx*diff.wt)/sum(wt*diff.wt^2)
}
#' @export
get.sd <- function(x.2, dfrac, mean, wt = rep(1, length(x.2))) {
  diff_p = t(diffseries_keepmean(t(wt*(x.2 - mean)), dfrac))
  sqrt(sum(wt*diff_p^2)/sum(wt))
}

#' @export
loglik_res_step = function(par, x.2, wt = rep(1, length(x.2))){

  dfrac = par[1]

  if (length(par) == 1) {
    mean <- get.m(x.2 = x.2, dfrac = dfrac, wt = wt)
    sd <- get.sd(x.2 = x.2, dfrac = dfrac, mean = mean, wt = wt)
  } else if (length(par) == 2) {
    mean <- par[2]
    sd <- get.sd(x.2 = x.2, dfrac = dfrac, mean = mean, wt = wt)
  } else if (length(par) == 3) {
    mean <- par[2]
    sd <- par[3]
  }

  logL = loglik_norm(x.2 = x.2, dfrac = dfrac, mean = mean, sd = sd, wt = wt)

  return(logL)
}
#' @export
loglik_t_res_step = function(par, x.2, wt = rep(1, length(x.2))) {
  dfrac = par[1]
  mean = par[2]
  sd = par[3]
  df = par[4]
  logL = loglik_t(x.2 = x.2, dfrac = dfrac, mean = mean, sd = sd, df = df, wt = wt)

  return(logL)
}
#' @export
loglik_t = function(x.2, dfrac, mean, sd, df, wt = rep(1, length(x.2))) {
  diff_p = c(diffseries_keepmean(wt*(x.2-mean), dfrac))
  sum(wt*ldt_mc(x = diff_p, mean = 0, sd = sd, df = df))
}
#' @export
loglik_norm = function(x.2, dfrac, mean, sd, wt = rep(1, length(x.2))) {
  diff_p = c(diffseries_keepmean(wt*(x.2-mean), dfrac))
  sum(wt*ldnorm_mc(x = diff_p, mean = 0, sd = sd))
}
#' @export
ldt_mc <- function(x, mean, sd, df) {
  # scale <- sd*sqrt((df - 2)/df)
  # dt((x - mean)/scale, df = df, log = TRUE) - log(scale)
  -log(pi)/2 - log(df - 2)/2 + log((gamma((df + 1)/2)/gamma(df/2))) +
    -(df + 1)*log(1 + (x - mean)^2/(sd^2*(df - 2)))/2 - log(sd^2)/2
}
#' @export
ldnorm_mc <- function(x, mean, sd) {
  -log(2*pi*sd^2)/2 - (x - mean)^2/(2*sd^2)
}
#' @export
fit_res_step = function(t.1, x.1, x.2, deg, seqby, resd.seqby, tau_j, N,
                        use_t = FALSE, start.2 = NULL, prof = TRUE){

  fit1 = refitWLS(t.1, x.1, deg = deg, seqby = seqby, resd.seqby = resd.seqby)
  resd1.1 = x.1 - fit1$fit.vals
  var.resd1.1 = fit1$var.resd
  ll.1 = sum(dnorm(resd1.1, log = TRUE, sd = sqrt(var.resd1.1)))

  if (!use_t | (use_t & is.null(start.2))) {

    if (is.null(start.2)) {
      start.norm <- 0
      if (!prof) {
        start.norm <- c(start.norm, mean(x.2), sd(x.2))
      }
    } else {
      start.norm <- start.2$d
      if (!prof) {
        start.norm <- c(start.norm, start.2$m, start.2$sd)
      }
    }

    optim.2 = optim(par = start.norm,
                    fn = loglik_res_step, method = "L-BFGS-B",
                    x.2 = x.2,
                    control = list("fnscale" = -1))
    fit_d = optim.2$par[1]
    if (length(optim.2$par) == 1) {
      fit_m <- get.m(x.2 = x.2, dfrac = fit_d)
      fit_sd <- get.sd(x.2 = x.2, dfrac = fit_d, mean = fit_m)
    } else if (length(optim.2$par) == 2) {
      fit_m <- optim.2$par[2]
      fit_sd <- get.sd(x.2 = x.2, dfrac = fit_d, mean = fit_m)
    } else {
      fit_m <- optim.2$par[2]
      fit_sd <- abs(optim.2$par[3])
    }
    ll.2 = optim.2$value
    fit_df <- NA
  }

  if (use_t) {

    if (is.null(start.2)) {
      start.t <- c(fit_d, fit_m, fit_sd, 3)
    } else {
      start.t <- c(start.2$d, start.2$m, start.2$sd, start.2$df)
    }

    opt <- optim(start.t, loglik_t_res_step, x.2 = x.2,
                 lower = c(-Inf, -Inf, -Inf, 2 + 10^(-14)), upper = c(Inf, Inf, Inf, 300),
                 control = list("fnscale" = -1), method = "L-BFGS-B")
    ll.2 = opt$value
    fit_d <- opt$par[1]
    fit_m <- opt$par[2]
    fit_sd <- abs(opt$par[3])
    fit_df <- opt$par[4]
  }
  return(list(fit_M=ll.1 + ll.2,fit_d=fit_d,fit_m=fit_m, fit_sd=fit_sd, fit_df = fit_df))
}

# fit first regime with model that accounts to heteroscedastic noise
# fit second regime with ARFIMA as in cpsearch2.R
# uses loglikelihood of original series

# T2CD-step method
# wrapper around search_dtau_step to fit with both ranges of d
# and pick the fit with highest likelihood
#' @export
t2cd_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                     seqby = 1, resd.seqby = 5,
                     use_scale = TRUE, use_t = FALSE, prof = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # deg: degree for B-spline
  # segby, resd.seqby: interval between knots
  # use_scale: if true, scale time series
  res1 = search_dtau_step(dat, t.max, tau.range, deg,
                          seqby = seqby, resd.seqby = resd.seqby,
                          use_scale = use_scale, use_t = use_t, prof = prof)
  return(res1)
}

### helper and plotting functions

# grid search for d and tau
# return the likelihood for each combination
#' @export
search_dtau_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                            seqby = 1, resd.seqby = 5,
                            use_scale = TRUE, use_t = FALSE, prof = TRUE){
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

  # iterate through each tau, return log-likelihood
  foreach.tau <- vector("list", length = length(tau.idx))
  for (j in 1:length(tau.idx)) {
    cat("j=", j, "\n")
    tau_j = tau.idx[j]

    # optimize for polynomial component
    x.1 = res_mean[1:tau_j]
    t.1 = tim[1:tau_j]
    n.1 = length(x.1)
    # optimize for ARFIMA
    x.2 = res_mean[(tau_j+1):N]

    if (j == 1) {
      start.2 <- NULL
    } else {
      start.2 <- list("d" = foreach.tau[[j-1]]$d,
                      "m" = foreach.tau[[j-1]]$m,
                      "sd" = foreach.tau[[j-1]]$sd,
                      "df" = foreach.tau[[j-1]]$df)
    }
    fit_res = tryCatch(fit_res_step(t.1 = t.1,  x.1 = x.1, x.2 = x.2, deg = deg,
                                    seqby = seqby, resd.seqby = resd.seqby, tau_j = tau_j,
                                    N = N, use_t = use_t,
                                    start.2 = start.2, prof = prof),
                       error = function(e) {return(NA)})
    if ((!use_t & any(is.na(fit_res[1:4]))) | (use_t & any(is.na(fit_res)))) {
      foreach.tau[[j]] <- list("M" = -Inf, "d" = NA, "m" = NA, "df" = NA, "sd" = NA)
    }else{
      foreach.tau[[j]] <- list("M" = fit_res$fit_M, "d" = fit_res$fit_d, "m" = fit_res$fit_m,
                  "sd" = fit_res$fit_sd,  "df" = fit_res$fit_df)
    }
  }

  M <- unlist(lapply(foreach.tau, function(x) {x$M}))
  d <- unlist(lapply(foreach.tau, function(x) {x$d}))
  m <- unlist(lapply(foreach.tau, function(x) {x$m}))
  sd <- unlist(lapply(foreach.tau, function(x) {x$sd}))
  df <- unlist(lapply(foreach.tau, function(x) {x$df}))

  # tau and d at maximum log-likelihood
  M_df = data.frame(tau = tim[tau.idx], M = M, d = d, m = m, sd = sd, df = df)
  max.idx = which.max(M)
  max.tau = tim[tau.idx[max.idx]]
  max.d = d[max.idx]
  max.m = m[max.idx]
  max.sd = sd[max.idx]
  max.df = df[max.idx]

  return(list(M_df = M_df, res = res, tim = tim, tau.idx = tau.idx,
              tau = max.tau, d = max.d, m = max.m, sd = max.sd, df = max.df,
              idx = tau.idx[max.idx], logL = M[max.idx]))
}
# plot sequences and fitted lines
#' @export
plot_t2cd_step = function(results, tau.range = c(10, 50), deg = 3,
                          use_scale = TRUE, return_plot = TRUE){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
  if (use_scale){
    res_mean = scale(res, center = F) # scaling
  }else{
    res_mean = matrix(res, ncol = 1)
  }
  N = length(res)

  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_idx = results$idx

  ### fitted values for first regime
  fit1 = refitWLS(tim[1:opt_idx], res_mean[1:opt_idx], deg = deg)
  var.resd1.1 = fit1$var.resd
  mu = c(fit1$fit.vals, rep(NA, N-opt_idx))

  # fitted values for second regime
    m = results$m
    x.2 = res_mean[(opt_idx+1):N]
    diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), opt_d))

    mu.2 = (x.2 - m) - diff_p
    mu[(opt_idx+1):N] = mu.2 + m

  if (use_scale){
    fit.vals = mu*attributes(res_mean)$'scaled:scale'
    var.resd1 = fit1$var.resd*attributes(res_mean)$'scaled:scale'^2
  }else{
    fit.vals = mu
    var.resd1 = fit1$var.resd
  }

  # plotting
  if (return_plot){
    plot(tim, res, ylim = c(min(c(res, fit.vals)), max(c(res, fit.vals))), type = 'l',
         main = paste('Values fitted with d: ', round(opt_d,3), ' tau: ', round(opt_tau,3)),
         xlab = 'Time (hour)', ylab = 'Resistance (ohm)')
    if (is.na(opt_tau)){
      lines(tim, fit.vals, col = "blue", lwd = 1)
    }else{
      opt_tau.idx = which(tim == opt_tau)
      lines(tim[1:opt_idx], fit.vals[1:opt_idx], col = "blue", lwd = 1)
      lines(tim[(opt_idx+1):N], fit.vals[(opt_idx+1):N], col = "green", lwd = 1)
      abline(v = opt_tau, lty = 2, col = "red")
    }
    abline(v = tau.range, lty = 1, col = "red")
  }

  return(list(fit.vals1 = mu[1:opt_idx], fit.vals2 = mu[(opt_idx+1):N],
              opt_idx = opt_idx, N = N, scaling = attributes(res_mean)$'scaled:scale',
              var.resd1 = fit1$var.resd))
}

# parametric bootstrap using outputs from t2cd_step and plot.t2cd_step
# Need to update for t-distribution
#' @export
bootstrap_sample_step = function(results, plot_results, seed = 0){

  set.seed(seed)
  res = results$res
  tim = results$tim
  N = length(res)
  opt_idx = results$idx

  # regime 1
  fit.vals1 = plot_results$fit.vals1*plot_results$scaling
  var.resd1 = plot_results$var.resd1*plot_results$scaling^2
  noise1 = rnorm(opt_idx, 0, sqrt(var.resd1))

  # regime 2
  opt_d = results$d
  m = results$m
  res_mean = scale(res, center = F) # scaling
  diff_p = c(diffseries_keepmean(matrix(res_mean[(opt_idx+1):N]-m, ncol = 1), opt_d))
  sd.resd2 = sqrt(mean(diff_p^2))
  sim = sim.fi(N-opt_idx, opt_d, sd.resd2, mu = m)
  seq_fi = sim$s

  samp = c(fit.vals1 + noise1, seq_fi*plot_results$scaling)

  return(list(res=matrix(samp, nrow=1), tim=tim))
}
#' @export
boot_k = function(k, results, plot_results, methodname){

  if (methodname=='step'){
    samp = bootstrap_sample_step(results, plot_results, seed = k)
    res = t2cd_step(samp)
  }else if (methodname=='sigmoid'){
    samp = bootstrap_sample_sigmoid(results, plot_results, seed = k)
    res = t2cd_sigmoid(list(res=matrix(samp$res,1), tim=matrix(samp$tim,1)))
  }
  tau_step = res$tau
  d_step = res$d

  return(list(tau=tau_step, d=d_step))
}

