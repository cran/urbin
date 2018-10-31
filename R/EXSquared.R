EXSquared <- function( lowerBound, upperBound ) {
  result <- ( upperBound^3 - lowerBound^3 )/( 3 * ( upperBound - lowerBound ) )
  return( result )
}
