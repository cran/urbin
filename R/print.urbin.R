print.urbin <- function( x, ... ) {
  if( "semEla" %in% names( x ) ) {
    printVec <- c( semEla = x$semEla, stdEr = x$stdEr )
  } else if( "effect" %in% names( x ) ) {
    printVec <- c( effect = x$effect, stdEr = x$stdEr )
  } else {
    stop( "internal error: object of class 'urbin' does not include",
      " expected components" )
  }
  print( printVec )
  invisible( x )
}
