\name{pglm}
\alias{pglm}
\alias{ordinal}
\title{Panel Estimators for Generalized Linear Models}
\description{
  Estimation by maximum likelihood of glm (binomial and Poisson) and
  'glm-like' models (Negbin and ordered) on longitudinal data
}
\usage{
pglm(formula, data, subset, na.action,
     effect = c("individual", "time", "twoways"),
     model = c("random", "pooling", "within", "between"),
     family, index = NULL, start = NULL, R = 20,  ...) 
}
\arguments{
  \item{formula}{a symbolic description of the model to be estimated,}
  \item{data}{the data: a \code{pdata.frame} object or an ordinary
    \code{data.frame},}
  \item{subset}{an optional vector specifying a subset of observations,}
  \item{na.action}{a function which indicates what should happen when
    the data contains '\code{NA}'s,}
  \item{effect}{the effects introduced in the model, one of
    \code{"individual"}, \code{"time"} or \code{"twoways"},}
  \item{model}{one of \code{"pooling"}, \code{"within"},
    \code{"between"}, \code{"random",},}
  \item{family}{the distribution to be used,}
  \item{index}{the index,}
  \item{start}{a vector of starting values,}
  \item{R}{the number of function evaluation for the gaussian quadrature
  method used,}
  \item{...}{further arguments.}
}
\value{
  An object of class \code{"pglm"}, a list with elements:
  \item{coefficients}{the named vector of coefficients,}
  \item{logLik}{the value of the log-likelihood,}
  \item{hessian}{the hessian of the log-likelihood at convergence,}
  \item{gradient}{the gradient of the log-likelihood at convergence,}
  \item{call}{the matched call,}
  \item{est.stat}{some information about the estimation (time used,
    optimisation method),}
  \item{freq}{the frequency of choice,}
  \item{residuals}{the residuals,}
  \item{fitted.values}{the fitted values,}
  \item{formula}{the formula (a \code{mFormula} object),}
  \item{expanded.formula}{the formula (a \code{formula} object),}
  \item{model}{the model frame used,}
  \item{index}{the index of the choice and of the alternatives.}
}

% \details{
% }

\examples{
## an ordered probit example
#data('Fairness', package = 'pglm')
#op <- pglm(as.numeric(answer)~recurring+driving+education+rule+
#          education:rule, Fairness[Fairness$good == 'parking',],
#          family=ordinal('probit'), R=5, print.level=0,
#          method='bfgs', index='id', iterlim=2)

## a binomial (logit) example
data('Unions', package = 'pglm')
anb <- pglm(union~wage+exper+rural, Unions, family=binomial('probit'),
model="pooling",  method = "bfgs", print.level=0, R = 5, iterlim=2)
}

% \seealso{
% }

\author{Yves Croissant}

% \references{
% }

\keyword{regression}
