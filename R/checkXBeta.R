checkXBeta <- function( xBeta ) {
  if( any( abs( xBeta ) > 3.5 ) ) {
    warning( "At least one x'beta has an implausible value: ", 
      paste( xBeta, collapse = ", " ) )
  }  
}
