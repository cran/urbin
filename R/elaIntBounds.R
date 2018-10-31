elaIntBounds <- function( xBound, nInt, argName = "xBound" ) {
  if( length( xBound ) != nInt + 1 ) {
    stop( "argument '", argName, "' must be a vector with ", nInt + 1, 
      " elements" )
  }
  if( any( xBound != sort( xBound ) ) ) {
    stop( "the elements of the vector specified by argument '", argName, 
      "' must be in increasing order" )
  }
  if( max( table( xBound ) ) > 1 ) {
    stop( "the vector specified by argument '", argName, 
      "' may not contain two (or more) elements with the same value" )
  }
  if( is.infinite( xBound[ nInt + 1 ] & nInt > 1 ) ) {
    xBound[ nInt + 1 ] <- 3 * xBound[ nInt ] - 2 * xBound[ nInt - 1 ]
  }
  return( xBound )
}
