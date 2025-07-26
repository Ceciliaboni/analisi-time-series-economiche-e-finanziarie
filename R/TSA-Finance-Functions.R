################################################################################
##
## File:    TSA-Finance-Functions.R
## 
## Purpose: Useful functions for TSA of financial data.
##
## Created: 2020.11.05
##
## Version: 2023.01.06
## 
################################################################################


################################################################################
## Functions 
################################################################################

.dist <- 
function(fit, x, type)
{
  #### Values for the pdf
  dist <- fit@model$modeldesc$distribution
  x1 <- fit@fit$coef["lambda"]
  lambda <- ifelse(is.na(x1), -0.5, x1)
  x1 <- fit@fit$coef["skew"]
  skew <- ifelse(is.na(x1), 1, x1)
  x1 <- fit@fit$coef["shape"]
  shape <- ifelse(is.na(x1), 1, x1)
  #### Compute
  if (type == "d")
  {
    ddist(distribution = dist, y = x, mu = 0, sigma = 1, 
      lambda = lambda, skew = skew, shape = shape)
  }
  else if (type == "q")
  {
    qdist(distribution = dist, p = x, mu = 0, sigma = 1, 
      lambda = lambda, skew = skew, shape = shape)
  }
  else
  {
    stop("Argument type must be 'd' or 'q'")
  }
}
# ------------------------------------------------------------------------------


.hist <- 
function(x, xlim, n = 200, breaks = 100, xlab = "", main = "")
{
  #### Histogram
  hist(x = x, breaks = breaks, freq = FALSE, xlim = xlim, 
    xlab = xlab, main = main)
  #### Plot of the pdf
  x1 <- seq(from = xlim[1], to = xlim[2], length.out = n)
  pdf1 <- dnorm(x = x1, mean = mean(x), sd = sd(x))
  lines(x = x1, y = pdf1, col = "red", lwd = 2) 
}
# ------------------------------------------------------------------------------


.hist.fit <- 
function(fit, xlim, ylim = NULL, n = 200, breaks = 100, plot.norm = FALSE,
  main = "")
{
  #### Settings
  col <- c("red", "blue")
  #### Histogram
  hist(x = fit@fit$z, breaks = breaks, freq = FALSE, xlim = xlim, xlab = "z",
    ylim = ylim, main = main)
  #### Plot of the selected pdf
  x1 <- seq(from = xlim[1], to = xlim[2], length.out = n)
  pdf1 <- .dist(fit = fit, x = x1, type = "d")
  lines(x = x1, y = pdf1, col = col[1], lwd = 2)
  #### Add the Normal pdf
  if (plot.norm[1] & fit@model$modeldesc$distribution != "norm")
  {
    pdf2 <- dnorm(x = x1)
    lines(x = x1, y = pdf2, col = col[2], lwd = 2)
    #### Add legend
    legend <- c(fit@model$modeldesc$distribution, "norm")
    legend(x = "topright", legend = legend, 
      fill = NULL, col = col, border = "white", bty = "n", lty = 1)
  }
}
# ------------------------------------------------------------------------------


.qqplot.fit <- 
function(fit)
{
  #### Compute
  zemp <- sort(fit@fit$z)
  distr <- fit@model$modeldesc$distribution
  ####
  main <- paste0(distr, " Q-Q Plot")
  n <- NROW(zemp)
  p <- seq(from = 1 / (n + 1), by = 1 / (n + 1), length.out = n)  
  zth  <- .dist(fit = fit, x = p, type = "q")
  #### Plot
  qqplot(x = zth, y = zemp, plot.it = TRUE, xlab = "Theoretical quantiles",
    ylab = "Empirical quantiles", main = main)
  abline(a = 0, b = 1, col = "red")
}
# ------------------------------------------------------------------------------

.Eabsz <- function(fit)
{
  #### Distribution type
  dist <- fit@model$modeldesc$distribution
  #### Cases
  ans <- if (dist == "norm")
  {
    sqrt(2 / pi)
  }
  else if (dist == "std")
  {
    df <- as.numeric(fit5@fit$coef["shape"])
    2 * sqrt(df - 2) / ((df - 1) * beta(0.5 * df, 0.5))
  }
  else
  {
    paste0("Not implemented for distribution '", dist, "'") 
  }
  #### Answer
  ans
}
# ------------------------------------------------------------------------------


################################################################################
## Family-GARCH
################################################################################

.fgarch.2.gjr <-
function(fit)
{
  #### Extract
  vmodel <- fit@model$modeldesc$vmodel
  vsubmodel <- fit@model$modeldesc$vsubmodel

  ####
  if ( vmodel != "fGARCH" )
  {
    stop(".fgarch.2.gjr() needs an fGARCH model")
  }
  else
  {
    #### Extract
    est.all <- fit@fit$coef
    se.all  <- fit@fit$matcoef[, " Std. Error"]
    name.all <- names(est.all)
    name <- name.all[!is.na(se.all)]
    np.all <- NROW(name.all)
    np  <- NROW(name)
    est <- est.all[name]
    vcov <- fit@fit$cvar
    vcovR <- fit@fit$robust.cvar
    matcoef <- fit@fit$matcoef
    matcoefR <- fit@fit$robust.matcoef

    #### Extract
    ## alpha
    pattern <- "alpha"
    ind  <- which( substr(x = name, 1, nchar(pattern)) == pattern )
    alpha <- est[ind]
    inda <- ind
    ## eta1
    pattern <- "eta1"
    ind  <- which( substr(x = name.all, 1, nchar(pattern)) == pattern )
    inde <- ind
    eta1 <- est.all[ind]
    if (NROW(ind) > 0) { name.gamma <- paste0("gamma", 1 : NROW(ind)) } 
    else { name.gamma <- NULL}

    ####
    if (vsubmodel == "GJRGARCH")
    {
      #### Extract
      alpha.s <- alpha * (1 - eta1)^2
      gamma.s <- 4 * alpha * eta1
      ####
      D <- diag(np)
      D[cbind(inda, inda)] <- (1 - eta1)^2
      if (NROW(inde) > 0)
      {  
        D[cbind(inda, inde)] <- -2 * alpha * (1 - eta1)
        D[cbind(inde, inda)] <- 4 * eta1
        D[cbind(inde, inde)] <- 4 * alpha
      }
    }
    else if (vsubmodel == "TGARCH")
    {
      #### Extract
      alpha.s <- alpha * (1 - eta1)
      gamma.s <- 2 * alpha * eta1
      ####
      D <- diag(np)
      D[cbind(inda, inda)] <- 1 - eta1
      if (NROW(inde) > 0)
      {  
        D[cbind(inda, inde)] <- - alpha
        D[cbind(inde, inda)] <- 2 * eta1
        D[cbind(inde, inde)] <- 2 * alpha
      }
    }  
    ####
    est.all[inda] <- alpha.s
    est.all[inde] <- gamma.s
    names(est.all)[inde] <- name.gamma
    vcov  <- tcrossprod(D %*% vcov, D)
    vcovR <- tcrossprod(D %*% vcovR, D)

    ####
    vcov.all <- vcovR.all <- matrix(NA, np.all, np.all)
    ind <- which(name.all %in% name)
    vcov.all[ind, ind] <- vcov
    vcovR.all[ind, ind] <- vcovR

    ####
    x1 <- matcoef
    rownames(x1) <- names(est.all)
    se.all <- sqrt( abs( diag(vcov.all) ) )
    t.all  <- est.all / se.all
    pv.all <- 2 * ( 1 - pnorm( abs(t.all) ) )
    x1[, " Estimate"] <- est.all
    x1[, " Std. Error"] <- se.all
    x1[, " t value"] <- t.all
    x1[, "Pr(>|t|)"] <- pv.all
    matcoef <- x1
    ####
    x1 <- matcoefR
    rownames(x1) <- names(est.all)
    se.all <- sqrt( abs( diag(vcovR.all) ) )
    t.all  <- est.all / se.all
    pv.all <- 2 * ( 1 - pnorm( abs(t.all) ) )
    x1[, " Estimate"] <- est.all
    x1[, " Std. Error"] <- se.all
    x1[, " t value"] <- t.all
    x1[, "Pr(>|t|)"] <- pv.all
    matcoefR <- x1
  }

  #### Answer
  list(coef = matcoef[, " Estimate"],
    se.coef = matcoef[, " Std. Error"], tval = matcoef[, " t value"],
    robust.se.coef = matcoefR[, " Std. Error"],
    robust.tval = matcoefR[, " t value"],
    matcoef = matcoef, robust.matcoef = matcoefR,
    cvar = vcov, robust.cvar = vcovR)
}
# ------------------------------------------------------------------------------


################################################################################
## Garman-Klass
################################################################################

.garmanklass <-
function(data, 
  sd = TRUE, currency = FALSE)
{
  #### Auxiliary
  currency <- currency[1]
  nobs  <- NROW(data)

  #### Fix columns
  ind <- which(colnames(data) == "Price")
  if ( NROW(ind) > 0 )
  {
    colnames(data)[ind] <- "Close"
  }
  ind <- which(colnames(data) == "Adjusted")
  if ( NROW(ind) == 0 )
  {
    data$Adjusted <- data$Close
  }
  
  
  #### Intradaily
  ## Extract
  coef <- if (!currency)
  { 
    data$Adjusted / data$Close
  }
  else 
  { 
    1
  }
  H1 <- log( data$High * coef )
  L1 <- log( data$Low * coef )
  O1 <- log( data$Open * coef )
  C1 <- log( data$Close * coef )
  u1 <- H1 - O1
  d1 <- L1 - O1
  c1 <- C1 - O1
  ## Values
  x <- 0.511 * (u1 - d1)^2 +
    (-0.019) * (c1 * (u1 + d1) - 2 * u1 * d1) +
    (-0.383) * c1^2
  # x <- 0.5 * (H1 - L1)^2 - (2 * log(2) - 1) * (C1 - O1)^2
  #### Overnight adjustment
  if ( !currency )
  {
    retco <- c(NA, log( data$Open[-1] / data$Close[-nobs] ) )
    retoc <- log( data$Close / data$Open )
    x1 <- sum( retco^2, na.rm = TRUE); x2 <- sum( retoc^2, na.rm = TRUE )  
    f  <- x1 / (x1 + x2)
    f[f < 0.01] <- 0.01; f[f > 0.99] <- 0.99
    a <- 0.12
    x <- a * retco^2 / f + ( (1 - a) / (1 - f) ) * x
  }
  
  #### Answer
  if ( sd ) { 1.034 * sqrt( x ) } else { x }
}
# ------------------------------------------------------------------------------

################################################################################
## Predict
################################################################################

########
## Arguments:
## object:        an object created by ugarchfit()
## n.ahead:       number of steps ahead
## t:             time from which to start forecasts (it must be <= length(y))
## fixed.n.ahead: whether the number of steps ahead is retained fixed
## data:          Data used to compute predictions. 
##                Use this argument if we want to get forecasts beyond the 
##                estimation period (in this case it must be 
##                NROW(data) > NROW(object@model$modeldata$data)).  
##
########
.GARCH.predict <- 
function(object, n.ahead, t, 
  data = NULL, fixed.n.ahead = TRUE, alpha = NULL)
{
  #### Extract
  y <- if (NROW(data) == 0)
  {
    xts(x = object@model$modeldata$data, 
      order.by = object@model$modeldata$index)
  }
  else
  {
    data
  }
  nobs <- NROW(y)
  tStart <- t
  alpha <- alpha[1]
  
  #### Check
  if ( n.ahead <= 0 )
  {
    stop("Argument 'n.ahead' must be a positive integer")
  } 
  if ( t > nobs )
  {
    stop("Argument 't' must be <= NROW(y)")
  }
  if ( NROW(alpha) > 0 && alpha >= 0.5 )
  {
    stop("Argument 'alpha' must be lower than 0.5")
  }

  #### Adjust
  n.ahead <- round(n.ahead[1])

  #### Manage the model to rebuild the specification
  spec <- getspec(object)
  setfixed(spec) <- as.list(coef(object))  
  
  #### n.ahead not fixed
  if (!fixed.n.ahead)
  {
    #### Forecast
    forc <- ugarchforecast(fitORspec = spec, n.ahead = n.ahead, 
      n.roll = 0, data = y, out.sample = nobs - t)
    #### Extract
    time1 <- tStart : (tStart + n.ahead - 1) 
    mu    <- forc@forecast$seriesFor
    sigma <- forc@forecast$sigmaFor
    if ( NROW(alpha) > 0 )
    {
      q1 <- quantile(forc, probs = alpha / 2)
      q2 <- quantile(forc, probs = 1 - alpha / 2)
    }
  }
  else
  {
    #### Forecast
    forc <- ugarchforecast(fitORspec = spec, n.ahead = n.ahead, 
      n.roll = nobs - t, data = y, out.sample = nobs - t)
    #### Extract
    time1 <- (tStart - 1 + n.ahead) : (nobs - 1 + n.ahead)
    mu    <- forc@forecast$seriesFor[n.ahead, ]
    sigma <- forc@forecast$sigmaFor[n.ahead, ]
    if ( NROW(alpha) > 0 )
    {
      q1 <- quantile(forc, probs = alpha / 2)[n.ahead, ]
      q2 <- quantile(forc, probs = 1 - alpha / 2)[n.ahead, ]
    }
  }

  #### Join
  pred1 <- cbind(t = time1, pred = as.numeric(mu), 
    se = as.numeric(sigma))
  if ( NROW(alpha) > 0 )
  {
    pred1 <- cbind(pred1, left = q1, right = q2)
  }
  
  #### Answer
  list( n.ahead = n.ahead, fixed.n.ahead = fixed.n.ahead, alpha = alpha, 
    pred = as.data.frame(pred1) )
}
# ------------------------------------------------------------------------------
