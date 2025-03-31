\name{predict_Lasso_SIR}
\alias{predict_Lasso_SIR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction based on Lasso SIR
}
\description{
This function calculates the value of \eqn{X\beta} for a new data set.
}
\usage{
predict_Lasso_SIR( lassosirobj, newdata )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ lassosirobj}{ LassoSIR object when running the function LassoSIR.
%%     ~~Describe \code{x} here~~
}
  \item{ newdata}{ A data frame consisting of the values of the predictors.
%%     ~~Describe \code{x} here~~
}

}
\details{
Based on the estiamted \eqn{\beta}, this function provides the value \eqn{X\beta} for any new input.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
predict_value = predict_value, beta = lassosirobj$beta, no.dim = lassosirobj$no.dim
  \item{predict_value }{the value of \eqn{X\beta}.}
  \item{beta }{the estiamted value of the \eqn{\beta}.}
  \item{no.dim}{the dimension of the central space.}
%% ...
}
\references{
  Lin, Q., Zhao, Z. , and Liu, J. (2017) On consistency and sparsity for
  sliced inverse regression in high dimension.   Annals of
  Statistics.

  Lin, Q., Zhao, Z. , and Liu, J. (2016) Sparse Sliced Inverse
  Regression for High Dimensional Data.
}
\author{
  Zhigen Zhao, Qian Lin, Jun S. Liu
}


\examples{

p <- 10
n <- 200


H <- 20
m <- n/H

beta <- array(0, c(p, 1) )
beta[1:3,1] <- rnorm(3, 0, 1)

X <- array(0, c(n, p ) )

rho <- 0.3
Sigma <- diag(p)
elements <- rho^(c((p-1):0,1:(p-1) ) )
for(i in 1:p)
    Sigma[i,] <- elements[(p+1-i):(2*p-i) ]


X <- matrix( rnorm(p*n), c(n, p) )
X <- X\%*\% chol(Sigma)

Y <-  ( X\%*\% beta )^3/2 + rnorm(n,0,1)
sir.lasso <- LassoSIR( X, Y, H, choosing.d="automatic",
          solution.path=FALSE, categorical=FALSE, nfolds=10,
          screening=FALSE)

res = predict_Lasso_SIR( sir.lasso, newdata=data.frame( matrix( rnorm(5*p), c(5, p) )  ) )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
