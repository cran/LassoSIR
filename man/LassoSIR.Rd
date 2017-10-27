\name{LassoSIR}
\alias{LassoSIR}
\title{
LassoSIR
}
\description{
  This function calculates the sufficient dimension reduction (SDR) space
  using the Sparse Sliced Inverse Regression  Via Lasso (Lasso-SIR).

  The input is a continuous design matrix X and a response vector Y
  which can be either continuous or categorical. X is arranged such that
  each column corresponds to one variable and each row corresponds
  to one subject.
  
  The function gives users options to choose (i) the dimension of the
  SDR space, (ii) screening based on the diagonal thresholding, (iii)
  the number of slices (H), and many others.
}
\usage{
LassoSIR(X, Y, H = 0, choosing.d = "automatic", solution.path = FALSE,
        categorical = FALSE, nfolds = 10, screening = TRUE, no.dim = 0)
}

\arguments{
  \item{X}{
    This argument is the continuous design matrix X. X is arranged such that
    each column corresponds to one variable and each row corresponds
    to one subject.
  }
  \item{Y}{
    The response vector Y, which can be either continuous or categorical.
  }
  \item{H}{
    The number of slices.
    
    (i) If the boolean variable "categorical" is true, H   is chosen as
    the number of categories automatically.
    
    (ii) If the response variable is continuous, namely, "categorical" is
    false, user need to specify the number of slices. If H is set as 0,
    the code will ask the user to enter the number of slices interactively;
    
    (iii) the default choice of H is zero.
  }
  \item{choosing.d}{
    This argument asks for the method of choosing the dimension of SDR. If
    no.dim is non zero, then choosing.d is set as "given". Otherwise,
    choosing.d can be set as "automatic" or "manual".
    
    When choosing.d is set as "manual", this function will calculate the
    eigenvalues of var(EX|Y) and plot these eigenvalues. After that, the user will be asked to enter
    the dimension interactively.
    
    When choosing.d is set as "automatic", the dimension will be
    determined automatically according to Algorithm 5 from the original
    paper.  
    
    The default option is "automatic".
  }
  \item{solution.path}{
    When setting this boolean variable as TRUE, a plot with solution path
    based on the final proposed model will be plotted.
    
    The default option is FALSE.
  }
  \item{categorical}{
    When setting this boolean variable as TRUE, the response variable is
    categorical; otherwise, the response variable is continuous.
    
    The default option is FALSE.
  }
  \item{nfolds}{
    This argument set the number of folds in the cross validation. The
    default option is 10.
  }
  \item{screening}{
    When setting this boolean variable as TRUE, a diagonal thresholding (DT-SIR)
    step is applied to reduce the dimension before applying Lasso-SIR.
  }
  \item{no.dim}{
    This argument specifies the dimension of SDR. The default option is 0
    and this dimension is chosen manually or automatically based on the
    choice of choosing.d.
  }
}
\details{
  This function estimates the sufficient dimension reduction space using
  the sparse sliced inverse regression for high dimensional data via
  Lasso (LassoSIR). 
}
\value{
  When solution.path is set as true, the function returns a glmnet
  object.
  
  When solution.path is set as false, the tuning parameter in Lasso is
  chosen by using the cross validation. The function returns the
  following values:
  
  \item{beta}{the estimated coefficient in SDR.}

  \item{eigen.value}{the eigen value of the estimator of \eqn{var(EY|X)}.}

  \item{no.dim}{the dimension of the central space.}

  \item{H}{the number of slices.}
  
  \item{categorical}{a boolean variable to indicate the type of the response.}
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
\note{
  NA
}
\seealso{
  NA
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

}
\keyword{ Lasso }% use one of  RShowDoc("KEYWORDS")
\keyword{ SIR }% __ONLY ONE__ keyword per line
\keyword{ Optimal}
\keyword{ Model free}
