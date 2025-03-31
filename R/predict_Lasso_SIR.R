predict_Lasso_SIR <- function( lassosirobj, newdata )
{
  n = dim(newdata)[1]
  p = dim(newdata)[2]
  if( is.data.frame(newdata) )
    newdata = array( unlist( data.frame(newdata)), c(n,p))

  predict_value = newdata %*% lassosirobj$beta

  list( predict_value = predict_value, beta = lassosirobj$beta, no.dim = lassosirobj$no.dim )
}
