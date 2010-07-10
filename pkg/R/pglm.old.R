ordinal <- function(link = c('probit', 'logit')){
  link <- match.arg(link)
  list(family = 'ordinal', link = link)
}

negbin <- function(link = c('log'), vlink = c('nb1', 'nb2')){
  link <- match.arg(link)
  vlink <- match.arg(vlink)
  list(family = 'negbin', link = link, vlink = vlink)
}

pglm <-  function(formula, data, subset, na.action,
                  effect = c('individual','time','twoways'),
                  model  = c('random', 'pooling', 'within', 'between'),
                  family, index  = NULL, start = NULL, R = 20, ...){
  dots <- list(...)
  nframe <- length(sys.calls())
  args <- list(model = model, effect = effect)

  if (is.character(family)){
    if (family %in% c("ordprobit", "ordlogit", "tobit")){
      if (family == "ordprobit") family <- list(family = "ordinal", link = "probit")
      if (family == "ordlogit") family <- list(family = "ordinal", link = "logit")
      if (family == "tobit") family <- list(family = "tobit", link = NULL)
    }
    else  family <- get(family, mode = "function")
  }
  if (is.function(family)) family <- family()
  
  link <- thelink <- family$link
  if (family$family == "negbin") vlink <- family$vlink
  family <- family$family
  
  # check and match the arguments
  effect <- match.arg(effect)
  if (!any(is.na(model))) model <- match.arg(model)

  # Check whether data is a pdata.frame and if not create it ; ignore
  # this step if model = pooling and if index = NULL so that ordinary
  # glm models can be fitted
  if (model != "pooling" | !is.null(index)){
    if (inherits(data, "pdata.frame") && !is.null(index))
      warning("the index argument is ignored because data is a pdata.frame")
    if (!inherits(data, "pdata.frame")) data <- pdata.frame(data, index)
    # Create a Formula object if necessary
    if (!inherits(formula, "pFormula")) formula <- pFormula(formula)
  }
  
  # eval the model.frame
  cl <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset", "na.action"),names(mf),0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf$formula <- formula
  mf$data <- data
  data <- eval(mf, parent.frame())
  # return the model.frame or estimate the model
  if (is.na(model)){
    attr(data, "formula") <- formula
    return(data)
  }
  Kw <- NULL
  if (model != "pooling"){
    X <- model.matrix(formula, data, rhs = 1, model = "pooling", effect = effect)
    if (model == "within" && family == "within"){
      Xw <- model.matrix(formula, data, rhs = 1, model = "within", effect = effect)
      Kw <- colnames(Xw)
      X <- X[, Kw, drop = FALSE]
    }
    if (ncol(X) == 0) stop("empty model")
    y <- pmodel.response(formula, data, model = "pooling", effect = effect)
    id <- attr(data, "index")[[1]]
  }
  else{
    X <- model.matrix(formula, data)
    y <- model.response(data)
    if (inherits(data, "pdata.frame")) id <- attr(data, "index")[[1]]
    else id <- NULL
  }
  # compute the starting vgalues
  start <- starting.values(family, model, Kw, X, y, cl, start)
  if (model == "random") rn <- gauss.quad(R, kind = 'hermite')

  # call to maxLik with the relevant arguments
  ml <- cl
  m <- match(c("print.level", "ftol", "tol", "reltol", "gradtol", "steptol",
               "lambdatol", "qrtol", "iterlim", "fixed", "activePar", "method"),names(ml),0)
  ml <- ml[c(1, m)]

  if (family == "poisson"){
    ml$logLik <- function(start)
      lnl.poisson(start, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE)
    ml$grad <- function(start) 
      attr(lnl.poisson(start, y, X, id, model, link, rn, gradient = TRUE, hessian = FALSE),
           'gradi')
    ml$hess <- function(start) 
      attr(lnl.poisson(start, y, X, id, model, link, rn, gradient = TRUE, hessian = TRUE),
           'hessian')
  }
  if (family == "gaussian"){
    ml$logLik <- function(start)
      lnl.gaussian(start, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE)
##     ml$grad <- function(start) 
##       attr(lnl.poisson(start, y, X, id, model, link, rn, gradient = TRUE, hessian = FALSE),
##            'gradi')
##     ml$hess <- function(start) 
##       attr(lnl.poisson(start, y, X, id, model, link, rn, gradient = TRUE, hessian = TRUE),
##            'hessian')
  }
  if (family == "tobit"){
    ml$logLik <- function(start)
      lnl.tobit(start, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE)
    ml$grad <- function(start) 
      attr(lnl.tobit(start, y, X, id, model, link, rn, gradient = TRUE, hessian = FALSE),
           'gradi')
    ml$hess <- function(start) 
      attr(lnl.tobit(start, y, X, id, model, link, rn, gradient = TRUE, hessian = TRUE),
           'hessian')
  }
  if (family == "negbin"){
    ml$logLik <- function(start)
      lnl.negbin(start, y, X, id, model, link, vlink, rn, gradient = FALSE, hessian = FALSE)
    ml$grad <- function(start) 
      attr(lnl.negbin(start, y, X, id, model, link, vlink, rn,
                      gradient = TRUE, hessian = FALSE), 'gradi')
    ml$hess <- function(start) 
      attr(lnl.negbin(start, y, X, id, model, link, vlink, rn,
                      gradient = TRUE, hessian = TRUE), 'hessian')
  }

  argschar <- function(args){
    paste(as.character(names(args)), as.character(args),
          sep="=", collapse=",")
  }
  args <- list(param = "start",
               y = "y", X = "X", id = "id", model = "model", family = "link",
               rn = "rn", gradient = FALSE, hessian = FALSE)

  thefunc <- paste("function(start) lnl.", family, "(", argschar(args), ")", sep = "")
  ml$logLik <- eval(parse(text = thefunc))
  args$gradient <- TRUE
  thefunc <- paste("function(start) attr(lnl.", family, "(", argschar(args), "), \"gradient\")", sep = "")
  ml$grad <- eval(parse(text = thefunc))
  args$hessian <- TRUE
  thefunc <- paste("function(start) attr(lnl.", family, "(", argschar(args), "), \"hessian\")", sep = "")
  ml$hess <- eval(parse(text = thefunc))
  
  
##   if (family == "binomial"){
##     ml$logLik <- function(start)
##       lnl.binomial(start, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE)
##     ml$grad <- function(start)
##       attr(lnl.binomial(start, y, X, id, model, link, rn, gradient = TRUE, hessian = FALSE),
##            'gradi')
##     ml$hess <- function(start)
##       attr(lnl.binomial(start, y, X, id, model, link, rn, gradient = TRUE, hessian = TRUE),
##            'hessian')
##   }

  
  if (family == "ordinal"){
    ml$logLik <- function(start)
      lnl.ordinal(start, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE, start.sigma = FALSE)
    ml$grad <- function(start)
      attr(lnl.ordinal(start, y, X, id, model, link, rn, gradient = TRUE, hessian = FALSE, start.sigma = FALSE),
           'gradi')
    ml$hess <- function(start)
      attr(lnl.ordinal(start, y, X, id, model, link, rn, gradient = TRUE, hessian = TRUE, start.sigma = FALSE),
           'hessian')
  }

  ml$start <- start
  ml[[1]] <- as.name('maxLik')
  result <- eval(ml, sys.frame(which = nframe))
  ml$start.sigma <- TRUE
  ml$iterlim <- 0
  print(ml);stop()
  sigma <-  eval(ml, sys.frame(which = nframe))
  print(sigma);stop()
  result[c('call', 'args', 'model')] <- list(cl, args, data)
  result
}

