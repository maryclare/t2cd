# weighted least squares (WLS) for fitting regression b-spline (cubic) to data
# 1. fit spline by LS
# 2. estimate noise variance
# 3. refit spline by WLS weighted by 1/variance
# 4. reiterate till convergence with difference between consecutive iterations below tol

# refit WLS till mean squared differences between consecutive fits is below tol or maxint is reached
refitWLS = function(tim, res_mean, deg = 3, seqby = 1, resd.seqby = 5,
                    resd.fit = 'pen', tol = 1e-10, maxint = 2){
  # tim, res_mean: time indices and observations
  # deg: degree of spline
  # seqby: interval of knots by units in tim
  # resd.seqby: interval of knots for residuals by units in tim
  # resd.fit: if 'pen', fit with penalty on second derivative, else 'bs'
  # tol, maxint: stopping criteria
  fitcurr = fitWLS1(tim, res_mean, weights = NULL, deg = deg, 
                    seqby = seqby, resd.seqby = resd.seqby,
                    resd.fit = resd.fit)
  numint = 2
  msdiff = Inf
  while(numint <= maxint & msdiff > tol){
    fitnext = fitWLS1(tim, res_mean, weights = fitcurr$var.resd, deg = deg,
                      seqby = seqby, resd.seqby = resd.seqby,
                      resd.fit = resd.fit)
    # update
    numint = numint + 1
    msdiff = mean((fitnext$fit.vals - fitcurr$fit.vals)^2)
    var.resd = fitcurr$var.resd
    fitcurr = fitnext
  }
  # smooth variance fit is preferred
  return(list(fit.vals = fitcurr$fit.vals, var.resd = var.resd,
              model.vals = fitcurr$model.vals, model.resd = fitcurr$model.resd))
}

### helper function for 1 iteration of fit

fitWLS1 = function(tim, res_mean, weights = NULL, deg = 3, 
                   seqby = 1, resd.seqby = 5, resd.fit = 'pen'){
  if (is.null(weights)){
    weights = rep(1, length(tim))
  }
  knots = seq(0, ceiling(max(tim)), by = seqby)
  if (resd.fit == 'bs'){
    bsdeg = splines.bs(tim, knots = knots[2:(length(knots)-1)], degree = deg)
    model = lm(res_mean ~ bsdeg, weights = weights)
    fit.vals = model$fitted.values
    resd = model$residuals 
  }else if (resd.fit == 'pen'){ # penalty on second derivative
    model = smooth.spline(tim, res_mean, w = 1/weights, cv = TRUE,
                          all.knots = FALSE, nknots = length(knots))
    fit.vals = model$y
    resd = res_mean - model$y
  }
  # model the residuals with heteroscedastic noise
  if (resd.fit == 'bs'){
    # spline regression on log squared residuals
    resd.knots = seq(0, ceiling(max(tim)), by = resd.seqby)
    resd.bsdeg = splines.bs(tim, knots = resd.knots[2:(length(resd.knots)-1)], degree = deg)
    resd.model = lm(log(resd^2) ~ resd.bsdeg)
    var.resd = exp(resd.model$fitted.values)
  }else if (resd.fit == 'pen'){
    resd.knots = seq(0, ceiling(max(tim)), by = seqby)
    resd.model = smooth.spline(tim, log(resd^2), cv = TRUE,
                          all.knots = FALSE, nknots = length(resd.knots))
    var.resd = exp(resd.model$y)
  }
  return(list(fit.vals = fit.vals, var.resd = var.resd, 
              model.vals = model, model.resd = resd.model))
}


