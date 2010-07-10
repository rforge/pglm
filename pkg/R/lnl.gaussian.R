lnl.gaussian <- function(param, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE,
                         opposite = FALSE, direction = rep(0, length(param)),
                         initial.value = NULL, steptol = 1E-01){
  opposite <- ifelse(opposite, -1, +1)
  Ti <- as.numeric(table(id)[as.character(id)])
  N <- length(y)
  n <- length(unique(id))
  names.id <- as.character(unique(id))
  step <- 2
  Xb <- X$Between
  yb <- y$Between
  Xw <- X$within
  yw <- y$within
  Xp <- Xw + Xb
  yp <- yw + yb
  K <- ncol(Xb)
  repeat{
    step <- step / 2
    if (step < steptol) break
    beta <- param[1L:K] + step * direction[1L:K]
    gamma <- param[K+1L] + step * direction[K+1L]
    sig2 <- param[K+2L] + step * direction[K+2L]
    ep <- yp - as.numeric(crossprod(t(Xp), beta))
    eb <- yb - as.numeric(crossprod(t(Xb), beta))
    ew <- yw - as.numeric(crossprod(t(Xw), beta))
    SSRB <- sum(eb^2)
    SSRW <- sum(ew^2)
    SSRP <- sum(ep^2)
    lnL <- - 1 / 2 * log(2 * pi) - 1 / 2 * log(sig2) - 1 / (2 * Ti) * log(1 + Ti * gamma) -
      1 / 2 * ep^2 / sig2 + 1 / 2 * (Ti * gamma) / (1 + Ti * gamma) * eb^2
    lnL <- opposite * sum(lnL)

    lnL <- -N/2*log(2*pi)-N/2*log(sig2)-1/2*sum(log(1+Ti*gamma))-1/(2*sig2)*SSRP
    +1/2*sum(gamma/(1+Ti*gamma)*Ti^2*unique(eb)^2)
    
    if (is.null(initial.value) || lnL <= initial.value) break
  }
  if (gradient){
    gb <- ep / sig2 * Xp - Ti * gamma / (1 + Ti * gamma) * eb * Xb
    gg <- - 1 / (2  * (1 + Ti * gamma)) + Ti * eb^2 / (2 * (1 + gamma * Ti)^2)
    gs <- ep^2 / (2 * sig2^2) - 1 / (2 * sig2)
    gradi <- cbind(gb, gg, gs)
    attr(lnL, "gradient") <- opposite * gradi
  }
  lnL
}


lnl.gaussian <- function(param, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE,
                         opposite = FALSE, direction = rep(0, length(param)),
                         initial.value = NULL, steptol = 1E-01){
  type <- "sd"
  opposite <- ifelse(opposite, -1, +1)
  Ti <- as.numeric(table(id)[as.character(id)])
  T <- Ti[1]
  N <- length(y)
  n <- length(unique(id))
  names.id <- as.character(unique(id))
  Xb <- apply(X, 2, tapply, id, mean)[as.character(id), ]
  yb <- tapply(y, id, mean)[as.character(id)]
  step <- 2
  K <- ncol(X)
  repeat{
    step <- step / 2
    if (step < steptol) break
    beta <- param[1L:K] + step * direction[1L:K]
    if (type == "sd"){
      sigmu <- param[K+1L] + step * direction[K+1L]
      sigeps <- param[K+2L] + step * direction[K+2L]
      gamma <- sigmu^2 / sigeps^2
      sig2 <- sigeps^2
    }
    else{
      gamma <- param[K+1L] + step * direction[K+1L]
      sig2 <- param[K+2L] + step * direction[K+2L]
    }
    eb <- as.numeric(yb) - as.numeric(crossprod(t(Xb), beta))
    e <- y - as.numeric(crossprod(t(X), beta))
    SSRB <- T * sum(eb^2)
    SSRP <- sum(e^2)
    lnL <- -N/2*log(2*pi)-N/2*log(sig2)-n/2*log(1+T*gamma)-
      1/2*SSRP/sig2+1/2*n*T^2*gamma/(1+T*gamma)*SSRB
    lnL <- - 1 / 2 * (log(2 * pi) + log(sig2) + 1 / T * log(1 + gamma * T) +
                     e^2 / sig2 - gamma * T / (1 + gamma * T) * eb^2 / sig2)
    lnL <- sum(lnL)
    
    if (is.null(initial.value) || lnL <= initial.value) break
  }
  if (gradient){
    gb <- e / sig2 * X - T * gamma / (1 + T * gamma) * eb * Xb /sig2
    gg <- - 1 / (2  * (1 + T * gamma)) + T * eb^2 / (2 * (1 + gamma * Ti)^2) / sig2
    gs <- - 1 / (2 * sig2) + e^2 / (2 * sig2^2) -
      gamma * T / (1 + gamma * T) * eb^2 / (2 * sig2^2)
    gradi <- cbind(gb, gg, gs)
    if (type == "sd"){
      gmu <- 2 * sigmu / sigeps^2 * gg
      geps <- 2 * sigeps * gs - 2 * sigmu^2 / sigeps^3 * gg
      gradi <- cbind(gb, gmu, geps)
    }
    attr(lnL, "gradient") <- opposite * gradi
  }
  lnL
}
