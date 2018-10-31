urbinElaInt <- function( allCoef, allXVal, xPos, xBound, model, 
  allCoefVcov = NULL, iPos = 1, yCat = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # check allXVal and allCoef
  if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
    if( nXVal != nCoef ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    } 
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
  checkXPos( xPos, minLength = 2, maxLength = nCoef + 1, 
    minVal = 0, maxVal = ifelse( model == "mlogit", nXVal, nCoef ), 
    requiredVal = 0 )
  # check position of the intercept
  checkIPos( iPos, xPos, allXVal, model ) 
  # number of intervals
  nInt <- length( xPos ) 
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, length( allCoef ), xPos, 
    pCall = match.call() )
  # vector of shares of observations in each interval
  shareVec <- rep( NA, nInt )
  for( i in 1:nInt ){
    if( xPos[i] != 0 ) {
      shareVec[ i ] <- allXVal[ xPos[i] ] 
    }
  }
  if( any( shareVec[ xPos != 0 ] < 0 ) ) {
    stop( "all elements of argument 'allXVal'",
      " that are indicated by argument 'xPos'",
      " (i.e., the shares of observations in each interval)",
      " must be non-negative" )
  }
  if( sum( shareVec[ xPos != 0 ] ) > 1 + sqrt(.Machine$double.eps) ) {
    stop( "the sum of the elements of argument 'allXVal'",
      " that are indicated by argument 'xPos'",
      " (i.e., the shares of observations in each interval)",
      " must not be larger than one" )
  }
  shareVec[ xPos == 0 ] <- 1 - sum( shareVec[ xPos != 0 ] )
  # weights
  weights <- elaIntWeights( shareVec )
  # prepare calculation of semi-elasticity 
  if( model == "lpm" ) {
    pFun <- rep( 0, length( xPos ) )
    for( i in 1:length( xPos ) ) {
      if( xPos[i] != 0 ) {
        pFun[i] <- allCoef[ xPos[i] ]
      }
    }
  } else if( model %in% c( "probit", "oprobit", "logit" ) ){
    xBeta <- rep( NA, nInt )
    for( i in 1:nInt ){
      allXValTemp <- replace( allXVal, xPos, 0 )
      if( xPos[i] != 0 ) {
        allXValTemp[ xPos[i] ] <- 1
      }
      xBeta[ i ] <- sum( allCoef * allXValTemp )
    }
    checkXBeta( xBeta )
  } else if( model == "mlogit" ){
    xBeta <- matrix( NA, nrow = nInt, ncol = nYCat + 1 ) 
    for( p in 1:( nYCat + 1 ) ){
      for( i in 1:nInt ){
        allXValTemp <- replace( allXVal, xPos, 0 )
        if( xPos[i] != 0 ) {
          allXValTemp[ xPos[i] ] <- 1
        }
        xBeta[i,p] <- sum( mCoef[ ,p] * allXValTemp )
      }
      checkXBeta( xBeta[,p] )
    }
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # vector of probabilities of y=1 for each interval
  if( model %in% c( "probit", "oprobit" ) ){
    pFun <- pnorm( xBeta )
    dFun <- dnorm( xBeta )
  } else if( model == "logit" ){
    pFun <- exp( xBeta ) / ( 1 + exp( xBeta ) )
    dFun <- exp( xBeta ) / ( 1 + exp( xBeta ) )^2
  } else if( model == "mlogit" ){
    pFunMat <- exp( xBeta ) / rowSums( exp( xBeta ) )
    pFun <- rowSums( pFunMat[ , yCat, drop = FALSE ] )
  } else if( model != "lpm" ) {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # effect of increasing to the next higher interval
  effNextInt <- pFun[ -1 ] - pFun[ -nInt ]
  # expected proportions of observations that increase to the next interval
  shareNextInt <- 0.5 * xBound[ 2:nInt ] * 
    ( shareVec[ -nInt ] / ( xBound[ 2:nInt ] - xBound[ 1:(nInt-1) ] ) +
        shareVec[ -1 ] / ( xBound[ 3:(nInt+1 ) ] - xBound[ 2:nInt ] ) )
  # (approximate) semi-elasticity
  semEla <- sum( effNextInt * shareNextInt )
  # partial derivatives of semi-elasticities wrt coefficients
  if( model == "lpm" ) {
    derivCoef <- rep( 0, nCoef )
    derivCoef[ xPos[1] ] <- - shareNextInt[1]
    derivCoef[ xPos[nInt] ] <- shareNextInt[nInt-1]
    if( nInt > 2 ) {
      for( n in 2:( nInt-1 ) ) {
        derivCoef[ xPos[n] ] <- 
          shareNextInt[n-1] - shareNextInt[n]
      }
    }
  } else if( model %in% c( "probit", "oprobit", "logit" ) ){
    derivCoef <- rep( 0, nCoef )
    derivCoef[ -xPos ] <- 
      sum( ( dFun[ -1 ] - dFun[ -nInt ] ) * shareNextInt ) * 
      allXVal[ -xPos ]
    derivCoef[ xPos[1] ] <- - dFun[1] * shareNextInt[1]
    derivCoef[ xPos[nInt] ] <- dFun[nInt] * shareNextInt[nInt-1]
    if( nInt > 2 ) {
      for( n in 2:( nInt-1 ) ) {
        derivCoef[ xPos[n] ] <- dFun[n] * 
          ( shareNextInt[n-1] - shareNextInt[n] )
      }
    }
  } else if( model == "mlogit" ){
    derivCoef <- matrix( 0, nrow = nXVal, ncol = nYCat )
    for( p in 1:nYCat ){
      for( yCati in yCat ) {
        derivCoef[ -xPos, p ] <- derivCoef[ -xPos, p ] +
          sum( ( pFunMat[ -nInt, yCati ] * pFunMat[ -nInt, p ] -
              pFunMat[ -1, yCati ] * pFunMat[ -1, p ] -
              ( p == yCati ) * 
                ( pFunMat[ -nInt, yCati ] - pFunMat[ -1, yCati ] ) ) *
              shareNextInt ) * allXVal[ -xPos ]
        derivCoef[ xPos[1], p ] <- derivCoef[ xPos[1], p ] +
          ( pFunMat[ 1, yCati ] * pFunMat[ 1, p ] -
              ( p == yCati ) * pFunMat[ 1, yCati ] ) * shareNextInt[1]
        derivCoef[ xPos[nInt], p ] <- derivCoef[ xPos[nInt], p ] -
          ( pFunMat[ nInt, yCati ] * pFunMat[ nInt, p ] -
              ( p == yCati ) * pFunMat[ nInt, yCati ] ) * 
          shareNextInt[nInt-1]
        if( nInt > 2 ) {
          for( n in 2:( nInt-1 ) ) {
            derivCoef[ xPos[n], p ] <- derivCoef[ xPos[n], p ] +
              ( pFunMat[ n, yCati ] * pFunMat[ n, p ] -
                  ( p == yCati ) * pFunMat[ n, yCati ] ) * 
              ( shareNextInt[n] - shareNextInt[n-1] )
          }
        }
      }     
    }
    derivCoef <- c( derivCoef )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }

  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  
  # create object that will be returned
  result <- list()
  result$call <- match.call()
  result$allCoefVcov <- allCoefVcov
  result$derivCoef <- derivCoef
  result$semEla <- unname( semEla )
  result$stdEr <- semElaSE
  class( result ) <- "urbin"
  return( result )
}
