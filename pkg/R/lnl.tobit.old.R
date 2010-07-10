lnl.tobit <- function(param, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE,
                         opposite = FALSE, direction = rep(0, length(param)),
                         initial.value = NULL, steptol = 1E-01){
  mills <- function(x) dnorm(x) / pnorm(x)
  opposite <- ifelse(opposite, -1, +1)
  Ti <- table(id)
  N <- length(y)
  n <- length(unique(id))
  names.id <- as.character(unique(id))
  step <- 2
  K <- ncol(X)
  repeat{
    step <- step / 2
    if (step < steptol) break
    beta <- param[1L:K] + step * direction[1L:K]
    sigma <- param[K+1L] + step * direction[K+1L]
    sigma <- sqrt(sigma)
    Xb <- as.numeric(crossprod(t(X), beta))
    if (model == "pooling"){
      lnL <- numeric(length = N)
      ez <- - Xb[y == 0] /sigma
      ep <- (y - Xb)[y > 0] / sigma
      mz <- mills(ez)
      lnL[y == 0] <- log(pnorm(ez))
      lnL[y  > 0] <- - 0.5 * log(2 * pi) - log(sigma) - 0.5 *ep^2
      lnl <- opposite * sum(lnL)
    }
    if (model == "random"){
      smu <- param[K + 2L] + step * direction[K + 2L]
      Pitr <- lapply(rn$nodes,
                     function(z){
                       result <- numeric(length = N)
                       ez <- - (Xb[y == 0] + sqrt(2) * smu * z) /sigma
                       ep <- (y - Xb - sqrt(2) * smu * z)[y > 0] / sigma
                       result[y == 0] <- pnorm(ez)
                       result[y  > 0] <- dnorm(ep) / sigma
                       result
                     }
                     )
      Pir <- lapply(Pitr, function(x) tapply(x, id, prod))
      Li <- Reduce("+", mapply("*", Pir, rn$weights, SIMPLIFY = FALSE)) / sqrt(pi)
      lnl <- opposite * sum(log(Li))
    }
    if (is.null(initial.value) || lnl <= initial.value) break
  }
  if (gradient){
    if (model == "pooling"){
      gradi <- matrix(0, nrow = nrow(X), ncol = ncol(X) + 1)
      gradi[y == 0, 1L:K] <- - mz * X[y == 0, , drop = FALSE] / sigma
      gradi[y == 0, K+1L] <- - ez * mz  / (2 * sigma^2)
      gradi[y  > 0, 1L:K] <- ep * X[y  > 0, , drop = FALSE] / sigma
      gradi[y  > 0, K+1L] <- - (1 - ep^2) / (2 * sigma^2)
      gradi <- opposite * gradi
    }
    if (model == "random"){
      g <- Reduce("+",
                  mapply(
                         function(w, x, p){
                           ez <- - (Xb[y == 0] + sqrt(2) * smu * x) /sigma
                           ep <- (y - Xb - sqrt(2) * smu * x)[y > 0] / sigma
                           mz <- mills(ez)
                           gradi <- matrix(0, nrow = N, ncol = 2)
                           gradi[y == 0, 1] <- - mz / sigma
                           gradi[y == 0, 2] <- - ez * mz  / (2 * sigma^2)
                           gradi[y  > 0, 1] <- ep  / sigma
                           gradi[y  > 0, 2] <- - (1 - ep^2) / (2 * sigma^2)
                           gradi <- cbind(gradi, gradi[, 1] * sqrt(2) * x)
                           w * as.numeric(p[as.character(id)]) * gradi
                         },
                         rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE
                         )
                  )
      gradi <- opposite * cbind(g[, 1] * X, g[, 2:3]) /
        as.numeric(Li[as.character(id)])/ sqrt(pi)
    }
    attr(lnl, 'gradi') <- gradi
    attr(lnl, 'gradient') <- apply(gradi, 2, sum)
  }
  if (hessian){
    if (model == "pooling"){
      hbb <- hbs <- hss <- numeric(length = N)
      hbb[y == 0] <- - (ez + mz) * mz / sigma^2
      hbs[y == 0] <- mz * (1 - (ez + mz) * ez)/(2 * sigma^3)
      hss[y == 0] <- ez * mz * (3 - (ez + mz) * ez) / (4 * sigma^4)
      hbb[y  > 0] <- - 1 / sigma^2
      hbs[y  > 0] <- - ep / sigma^3
      hss[y  > 0] <- (1 - 2 * ep^2) / (2 * sigma^4)
      hbb <- crossprod(hbb * X, X)
      hbs <- apply(hbs * X, 2, sum)
      hss <- sum(hss)
      H <- opposite * rbind(cbind(hbb, hbs), c(hbs, hss))
    }
    if (model == "random"){
      H <- mapply(
                  function(w, x, p){
                    P <- as.numeric((p/Li)[as.character(id)])
                    sp <- as.numeric(p / Li)
                    ez <- - (Xb[y == 0] + sqrt(2) * smu * x) /sigma
                    ep <- (y - Xb - sqrt(2) * smu * x)[y > 0] / sigma
                    mz <- mills(ez)
                    gradi <- matrix(0, nrow = N, ncol = 2)
                    gradi[y == 0, 1] <- - mz / sigma
                    gradi[y == 0, 2] <- - ez * mz  / (2 * sigma^2)
                    gradi[y  > 0, 1] <- ep  / sigma
                    gradi[y  > 0, 2] <- - (1 - ep^2) / (2 * sigma^2)
                    gradi <- cbind(gradi[, 1] * X, gradi[, 2], gradi[, 1] * sqrt(2) * x)
                    gradi <- apply(gradi, 2, tapply, id, sum)
                    H1 <- crossprod(sqrt(sp) * gradi)
                    hbb <- hbs <- hss <- numeric(length = N)
                    hbb[y == 0] <- - (ez + mz) * mz / sigma^2
                    hbs[y == 0] <- mz * (1 - (ez + mz) * ez)/(2 * sigma^3)
                    hss[y == 0] <- ez * mz * (3 - (ez + mz) * ez) / (4 * sigma^4)
                    hbb[y  > 0] <- - 1 / sigma^2
                    hbs[y  > 0] <- - ep / sigma^3
                    hss[y  > 0] <- (1 - 2 * ep^2) / (2 * sigma^4)
                    hbb <- crossprod(hbb * cbind(X, sqrt(2) * x) * P,
                                     cbind(X, sqrt(2) * x))
                    hbs <- apply(hbs * cbind(X, sqrt(2)* x) * P, 2, sum)
                    hss <- sum(hss * P)
                    H2 <- rbind(cbind(hbb, hbs), c(hbs, hss))
                    mX <- H2[1L:K, K+1L]
                    sX <- H2[1L:K, K+2L]
                    mm <- H2[K+1L, K+1L]
                    ss <- H2[K+2L, K+2L]
                    H2[K+1L, K+1L] <- ss
                    H2[K+2L, K+2L] <- mm
                    H2[1L:K, K+1L] <- H2[K+1L, 1L:K] <- sX
                    H2[1L:K, K+2L] <- H2[K+2L, 1L:K] <- mX
                    (H1 + H2) * w / sqrt(pi)
                  },
                  rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE
                  )
      H <- Reduce("+", H) - crossprod(apply(gradi, 2, tapply, id, sum))
    }
    attr(lnl, "hessian") <- opposite * H
  }
  lnl
}
