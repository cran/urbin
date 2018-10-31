elaIntWeights <- function( xShares ) {
  nInt <- length( xShares )
  weights <- rep( NA, nInt - 1 )
  for( m in 1:(nInt-1) ){
    weights[m] <- ifelse( m == 1, 1, 0.5 ) * xShares[m] +
      ifelse( m+1 == nInt, 1, 0.5 ) * xShares[m+1]
  }
  if( abs( sum( weights ) - 1 ) > 1e-5 ) {
    stop( "internal error: weights do not sum up to one" )
  }
  return( weights )
}
