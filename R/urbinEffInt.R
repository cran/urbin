urbinEffInt <- function( allCoef, allXVal = NULL, xPos, refBound, intBound, model,
  allCoefVcov = NULL, xMeanSd = NULL, iPos = 1, yCat = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # check allXVal and allCoef
  if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
    if( model == "lpm" ) {
      if( !is.null( allXVal ) ) {
        warning( "argument allXVal is ignored for lpm models",
          " (set this argument to 'NULL' to avoid this warning)" )
      }
      allXVal <- rep( 0, nCoef )
      if( iPos != 0 ){
        allXVal[ iPos ] <- 1
      }
      allXVal[ xPos ] <- NA
      nXVal <- length( allXVal )
    }  
    if( nXVal != nCoef ){
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }  
    # if( any( !is.na( allXVal[ xPos ] ) ) ) {
    #   allXVal[ xPos ] <- NA
    #   warning( "values of argument 'allXVal[ xPos ]' are ignored",
    #     " (set these values to 'NA' to avoid this warning)" )
    # }
  } else if( model == "mlogit" ){
    # number of alternative categories of the dependent variable
    nYCat <- round( nCoef / nXVal )
    if( nCoef != nXVal * nYCat ) {
      stop( "length of argument 'allCoef' must be a multiple",
        " of the length of argument 'allXVal'" )
    } 
    # create matrix of coefficients
    mCoef <- matrix( allCoef, nrow = nXVal, ncol = nYCat )
    # add column for coefficients of the reference category
    mCoef <- cbind( mCoef, 0 )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # check argument yCat
  if( model == "mlogit" ) {
    checkYCat( yCat, nYCat, maxLength = nYCat + 1 ) 
    yCat[ yCat == 0 ] <- nYCat + 1
  } else if( !is.null( yCat ) ) {
    warning( "argument 'yCat' is ignored" )
  }
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, 
    maxVal = ifelse( model == "mlogit", nXVal, nCoef ) )
  # Check value of quadratic term in argument allXVal
  if( length( xPos ) == 2 ){
    if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
      stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
        " to the squared value of 'allXVal[ xPos[1] ]'" )
    }
  } 
  # check position of the intercept
  checkIPos( iPos, xPos, allXVal, model ) 
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd,
    nXVal = nXVal, iPos = iPos, pCall = match.call() )
  # check the boundaries of the intervals
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )

  # calculate xBars
  intX <- mean( intBound )
  refX <- mean( refBound ) 
  if( length( xPos ) == 2 ) {
    intX <- c( intX, EXSquared( intBound[1], intBound[2] ) )
    refX <- c( refX, EXSquared( refBound[1], refBound[2] ) )
  }
  if( length( intX ) != length( xPos ) || 
      length( refX ) != length( xPos ) ) {
    stop( "internal error: 'intX' or 'refX' does not have the expected length" )
  }
  # define X' * beta 
  if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
    intXbeta <- sum( allCoef * replace( allXVal, xPos, intX ) )
    refXbeta <- sum( allCoef * replace( allXVal, xPos, refX ) )
    if( model != "lpm" ) {
      checkXBeta( c( intXbeta, refXbeta ) )
    }
  } else if( model == "mlogit" ){
    intXbeta <- replace( allXVal, xPos, intX ) %*% mCoef
    refXbeta <- replace( allXVal, xPos, refX ) %*% mCoef
    checkXBeta( c( intXbeta, refXbeta ) )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  
  # calculate the effect
  if( model == "lpm" ) {
    eff <- intXbeta - refXbeta
  } else if( model %in% c( "probit", "oprobit" ) ) {
    eff <- pnorm( intXbeta ) - pnorm( refXbeta )
  } else   if( model == "logit" ){  
    eff <- exp( intXbeta )/( 1 + exp( intXbeta ) ) - 
      exp( refXbeta )/( 1 + exp( refXbeta ) )
  } else if( model == "mlogit" ){
    pFunRef <- exp( refXbeta ) / sum( exp( refXbeta ) )
    pFunInt <- exp( intXbeta ) / sum( exp( intXbeta ) )
    eff <- sum( pFunInt[ yCat ] - pFunRef[ yCat ] )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }

  # partial derivative of the effect w.r.t. all estimated coefficients
  if( model == "lpm" ) {
    derivCoef <- rep( 0, nCoef ) 
    derivCoef[ xPos ] <- intX - refX
  } else if( model %in% c( "probit", "oprobit" ) ) {
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( dnorm( intXbeta ) - dnorm( refXbeta ) ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] = dnorm( intXbeta ) * intX - 
      dnorm( refXbeta ) * refX
  } else if( model == "logit" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( intXbeta )/( 1 + exp( intXbeta ) )^2 - 
        exp( refXbeta )/( 1 + exp( refXbeta ) )^2 ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] <- exp( intXbeta )/( 1 + exp( intXbeta ) )^2 * intX - 
      exp( refXbeta )/( 1 + exp( refXbeta ) )^2 * refX
  } else if( model == "mlogit" ){
    derivCoef <- matrix( 0, nrow = nXVal, ncol = nYCat )
    for( p in 1:nYCat ){
      for( yCati in yCat ) {
        derivCoef[ -xPos, p ] <- derivCoef[ -xPos, p ] + 
          ( pFunRef[ yCati ] * pFunRef[ p ] - pFunInt[ yCati ] * pFunInt[ p ] 
            - ( p == yCati ) * ( pFunRef[ yCati ] - pFunInt[ yCati ] ) ) *
          allXVal[ - xPos ]
        derivCoef[ xPos, p ] <- derivCoef[ xPos, p ] + 
          ( pFunRef[ yCati ] * pFunRef[ p ] - 
              ( p == yCati ) * pFunRef[ yCati ] ) * refX -
          ( pFunInt[ yCati ] * pFunInt[ p ] -
              ( p == yCati ) * pFunInt[ yCati ] ) * intX
      }
    }
    derivCoef <- c( derivCoef )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  
  # object to be returned
  result <- list()
  result$call <- match.call()
  result$allCoefVcov <- allCoefVcov
  result$derivCoef <- derivCoef
  result$effect <- eff
  result$stdEr <- effSE
  class( result ) <- "urbin"
  return( result )
}
