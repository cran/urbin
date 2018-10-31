checkXPos <- function( xPos, minLength, maxLength, minVal, maxVal,
  requiredVal = NA ) {
  if( any( xPos != round( xPos ) ) ) {
    stop( "argument 'xPos' must be a vector of integers" )
  }
  if( length( xPos ) < minLength ) {
    stop( "argument 'xPos' must have a length equal to or larger than ", 
      minLength )
  }
  if( length( xPos ) > maxLength ) {
    stop( "argument 'xPos' must have a length smaller than or equal to ", 
      maxLength )
  }
  if( any( xPos < minVal ) ) {
    stop( "all elements of argument 'xPos' must be equal to or larger than ",
      minVal )
  }
  if( any( xPos > maxVal ) ) {
    stop( "all elements of argument 'xPos' must be smaller than or equal to ",
      maxVal )
  }
  if( max( table( xPos ) ) > 1 ) {
    stop( "all elements of argument 'xPos' may only occur once" )
  }
  if( !is.na( requiredVal ) ) {
    if( sum( xPos == requiredVal ) != 1 ) {
      stop( "argument 'xPos' must have exactly one element that is ",
        requiredVal )
    }
  }
}
