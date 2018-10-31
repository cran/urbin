checkYCat <- function( yCat, nYCat, maxLength = 1 ) {
  if( is.null( yCat ) ) {
    stop( "argument 'yCat' must be be specified for 'mlogit' models" )
  }
  if( any( yCat != round( yCat ) ) ) {
    stop( "argument 'yCat' must be a vector of integers" )
  }
  if( length( yCat ) > maxLength ) {
    stop( "argument 'yCat' must have a length smaller than or equal to ", 
      maxLength )
  }
  if( any( yCat < 0 ) ) {
    stop( "all elements of argument 'yCat' must non-negative" )
  }
  if( any( yCat > nYCat ) ) {
    stop( "all elements of argument 'yCat' must be smaller than or equal to ",
      nYCat )
  }
  if( max( table( yCat ) ) > 1 ) {
    stop( "all elements of argument 'yCat' may only occur once" )
  }
}
