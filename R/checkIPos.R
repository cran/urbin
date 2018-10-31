checkIPos <- function( iPos, xPos, allXVal, model ) {
  if( length( iPos ) != 1 ) {
    stop( "argument 'iPos' must be a single integer value" )
  }
  if( iPos != round( iPos ) ) {
    stop( "argument 'iPos' must be an integer" )
  }
  if( iPos < 0 ) {
    stop( "argument 'iPos' must be a non-negative value" )
  }
  if( iPos != 0 && iPos %in% xPos ) {
    stop( "argument 'xPos' must not indicate the intercept",
      " (as indicated by argument 'iPos')" )
  }
  if( iPos > length( allXVal ) ) {
    stop( "argument 'iPos' must be smaller than or equal to ",
      length( allXVal ) )
  }
  if( iPos > 0 & !is.null( allXVal ) & !is.list( allXVal ) ) {
    if( model == "oprobit" ) {
      if( allXVal[ iPos ] != -1 ) {
        stop( "the value of argument 'allXVal' indicated by argument 'iPos',",
          " i.e., allXVal[iPos], must be minus one (-1)",
          " for ordered probit models" )
      }
    } else {
      if( allXVal[ iPos ] != 1 ) {
        stop( "the value of argument 'allXVal' indicated by argument 'iPos',",
          " i.e., allXVal[iPos], must be one (1)" )
      }
    }
  }
}
