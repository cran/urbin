## -----------------------------------------------------------------------------
linEla <- function( xCoef, xVal, xCoefSE = rep( NA, length( xCoef ) ) ){
  if( ! length( xCoef ) %in% c( 1, 2 ) ) {
    stop( "argument 'xCoef' must be a scalar or vector with 2 elements" )
  }
  if( length( xCoef ) != length( xCoefSE ) ) {
    stop( "arguments 'xCoef' and 'xCoefSE' must have the same length" )
  }
  if( length( xVal ) != 1 || !is.numeric( xVal ) ) {
    stop( "argument 'xVal' must be a single numeric value" )
  }
  if( length( xCoef ) == 1 ) {
    xCoef <- c( xCoef, 0 )
    xCoefSE <- c( xCoefSE, 0 )
  }
  semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal ) * xVal
  derivCoef <- c( xVal, ifelse( xCoef[2] == 0, 0, 2 * xVal^2 ) )
  vcovCoef <- diag( xCoefSE^2 )
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
# model equation that is linear in x_k
ela1a <- linEla( 0.05, 23.4 )
ela1a
ela1b <- linEla( 0.05, 23.4, 0.001 )
ela1b
all.equal( ela1b, linEla( c( 0.05, 0 ), 23.4, c( 0.001, 0 ) ) )
# Example
# model equation that is quadratic in x_k
ela2a <- linEla( c( 0.05, -0.00002 ), 23.4 )
ela2a
ela2b <- linEla( c( 0.05, -0.00002 ), 23.4, c( 0.001, 0.00002 ) )
ela2b

## -----------------------------------------------------------------------------
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

## ----elaIntBounds-------------------------------------------------------------
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

## -----------------------------------------------------------------------------
linElaInt <- function( xCoef, xShares, xBound,
                       xCoefSE = rep( NA, length( xCoef ) ) ){
  nInt <- length( xCoef )
  if( nInt < 2 || !is.vector( xCoef ) ) {
    stop( "argument 'xCoef' must be a vector with at least two elements" )
  }
  if( length( xCoefSE ) != nInt ) {
    stop( "arguments 'xCoef' and 'xCoefSE' must be vectors of the same length" )
  }
  if( length( xShares ) != nInt ) {
    stop( "arguments 'xCoef' and 'xShares' must be vectors of the same length" )
  }
  if( any( xShares < 0 ) ) {
    stop( "all shares in argument 'xShares' must be non-negative" )
  }
  if( abs( sum( xShares ) - 1 ) > 0.015 ) {
    stop( "the shares in argument 'xShares' must sum to one" )
  }
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # weights
  weights <- elaIntWeights( xShares )
  # semi-elasticities 'around' each inner boundary and their weights
  semElas <- rep( NA, nInt - 1 )
  for( m in 1:(nInt-1) ){
    semElas[m] <- 2 * ( xCoef[ m+1 ] - xCoef[ m ] ) * xBound[ m+1 ] /
      ( xBound[m+2] - xBound[m] )
  }
  # (average) semi-elasticity
  semElaAvg <- sum( semElas * weights )
  # derivatives of the (average) semi-elasticity wrt the coefficients
  derivCoef <- rep( NA, nInt )
  derivCoef[1] <-
   -2 * weights[1] * xBound[2] / ( xBound[3] - xBound[1] )
  derivCoef[nInt] <-
    2 * weights[nInt-1] * xBound[nInt] / ( xBound[nInt+1] - xBound[nInt-1] )
  if( nInt > 2 ) {
    for( n in 2:( nInt-1 ) ) {
      derivCoef[n] <-
        2 * weights[n-1] * xBound[n] / ( xBound[n+1] - xBound[n-1] ) -
        2 * weights[n]   * xBound[n+1] / ( xBound[n+2] - xBound[n] )
    }
  }
  # variance-covariance matrix of the coefficiencts
  vcovCoef <- diag( xCoefSE^2 )
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla = semElaAvg, stdEr = semElaSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
ela3a <- linElaInt( xCoef = c( 0, 0.22, 0.05, 0.6 ),
  xShares = c( 0.35, 0.4, 0.12, 0.13 ),
  xBound = c( 0, 500, 1000, 1500, Inf ) )
ela3a
# Example
ela3b <- linElaInt( xCoef = c( 0, 0.22, 0.05, 0.6 ),
  xShares = c( 0.35, 0.4, 0.12, 0.13 ),
  xBound = c( 0, 500, 1000, 1500, Inf ),
  xCoefSE = c( 0, 0.002, 0.005, 0.001 ) )
ela3b

## -----------------------------------------------------------------------------
EXSquared <- function( lowerBound, upperBound ) {
  result <- ( upperBound^3 - lowerBound^3 )/( 3 * ( upperBound - lowerBound ) )
  return( result )
}

## -----------------------------------------------------------------------------
linEffInt <- function( xCoef, refBound, intBound,
                       xCoefSE = rep( NA, length( xCoef ) ) ){
  if( ! length( xCoef ) %in% c( 1, 2 ) ) {
    stop( "argument 'xCoef' must be a scalar or vector with 2 elements" )
  }
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )
  if( length( xCoef ) != length( xCoefSE ) ) {
    stop( "arguments 'xCoef' and 'xCoefSE' must have the same length" )
  }
  if( length( xCoef ) == 1 ) {
      xCoef <- c( xCoef, 0 )
      xCoefSE <- c( xCoefSE, 0 )
  }
  # difference between the xBars of the two intervals
  xDiff <- mean( intBound ) - mean( refBound )
  # difference between the xSquareBars of the two intervals
  xSquaredDiff <-
    EXSquared( intBound[1], intBound[2] ) -
    EXSquared( refBound[1], refBound[2] )
  # effect E_{k,ml}
  eff <-  xCoef[1] * xDiff + xCoef[2] * xSquaredDiff
  # partial derivative of E_{k,ml} w.r.t. the beta_k and beta_{k+1}
  derivCoef <- c( xDiff, ifelse( xCoef[2] == 0, 0, xSquaredDiff ) )
  # variance covariance of the coefficients (covariances set to zero)
  vcovCoef <- diag( xCoefSE^2 )
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # object to be returned
  result <- c( effect = eff, stdEr = effSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
# model equation that is linear in x_k
eff1a <- linEffInt( 0.4, refBound = c( 19, 34 ), intBound = c( 35, 55 ) )
eff1a
eff1b <- linEffInt( 0.4, refBound = c( 19, 34 ), intBound = c( 35, 55 ),
                 xCoefSE = 0.03 )
eff1b
# Example
# model equation that is quadratic in x_k
eff2a <- linEffInt( c( 0.4, -0.0003 ),
                 refBound = c( 19, 34 ), intBound = c( 35, 55 ) )
eff2a
eff2b <- linEffInt( c( 0.4, -0.0003 ),
                 refBound = c( 19, 34 ), intBound = c( 35, 55 ),
                 xCoefSE = c( 0.002, 0.000001 ) )
eff2b

## -----------------------------------------------------------------------------
linEffGr <- function( xCoef, xShares, Group,
                      xCoefSE = rep( NA, length( xCoef ) ) ){
  if( sum( xShares ) > 1 ){
    stop( "the shares in argument 'xShares' sum up to a value larger than 1" )
  }
  if( length( xCoef ) != length( xShares ) ){
    stop( "arguments 'xCoef' and 'xShares' must have the same length" )
  }
  if( length( xCoef ) != length( Group ) ){
    stop( "arguments 'xCoef' and 'Group' must have the same length" )
  }
  if( ! all( Group %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'Group' must be -1, 0, or 1" )
  }
  # D_mr
  DRef <- xShares * ( Group == -1 ) / sum( xShares[ Group == -1 ] )
  # D_ml
  DEffect <- xShares * ( Group == 1 ) / sum( xShares[ Group == 1 ] )
  # effect: sum of delta_m * ( D_ml - D_mr )
  effeG <- sum( xCoef * ( DEffect - DRef ) )
  # approximate standard error
  effeGSE <- sqrt( sum( ( xCoefSE * ( DEffect - DRef ) )^2 ) )
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example:
eff3a <- linEffGr( xCoef = c( 0, 0.2, 0.04, 0.06, 0.3 ),
                   xShares = c( 0.14, 0.35, 0.3, 0.01, 0.2 ),
                   Group = c( -1, 1, -1, -1, 0 ) )
eff3a
# Example:
eff3b <- linEffGr( xCoef = c( 0, 0.2, 0.04, 0.06, 0.3 ),
                   xShares = c( 0.14, 0.35, 0.3, 0.01, 0.2 ),
                   Group = c( -1, 1, -1, -1, 0 ),
                   xCoefSE = c( 0, 0.0001, 0.002, 0.05, 0.09 ))
eff3b
# Example:
eff3c <- linEffGr( xCoef = c( 0, 0.2, 0.04, 0.06 ),
                   xShares = c( 0.14, 0.35, 0.3, 0.01 ),
                   Group = c( -1, 1, -1, -1 ))
eff3c
# Example:
eff3d <- linEffGr( xCoef = c( 0, 0.2, 0.04, 0.06 ),
                   xShares = c( 0.14, 0.35, 0.3, 0.01 ),
                   Group = c( -1, 1, -1, -1 ),
                   xCoefSE = c( 0, 0.0001, 0.002, 0.05 ))
eff3d

## ----checkXPos,echo=FALSE-----------------------------------------------------
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

## ----checkXBeta,echo=FALSE----------------------------------------------------
checkXBeta <- function( xBeta ) {
  if( any( abs( xBeta ) > 3.5 ) ) {
    warning( "At least one x'beta has an implausible value: ",
      paste( xBeta, collapse = ", " ) )
  }
}

## -----------------------------------------------------------------------------
probEla <- function( allCoef, allXVal,
  allCoefSE = rep( NA, length( allCoef )), xPos ){
  nCoef <- length( allCoef )
  if( nCoef != length( allCoefSE ) ) {
    stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
  }
  if( length( allCoef ) != length( allXVal ) ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  if( length( xPos ) == 2 ){
    xCoef <- allCoef[ xPos ]
    xCoefSE <- allCoefSE[ xPos ]
    if( !isTRUE( all.equal( allXVal[ xPos[2] ], allXVal[ xPos[1] ]^2 ) ) ) {
      stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
        "to the squared value of 'allXVal[ xPos[1] ]' " )
    }
  } else if( length( xPos ) == 1 ) {
    xCoef <- c( allCoef[ xPos ], 0 )
    xCoefSE <- c( allCoefSE[ xPos ], 0 )
  } else {
    stop( "argument 'xPos' must be a scalar or a vector with two elements" )
  }
  xVal <- allXVal[ xPos[ 1 ] ]
  xBeta <- sum( allCoef * allXVal )
  checkXBeta( xBeta )
  dfun <- pnorm( xBeta )
  semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  derivCoef <- c( dfun * xVal,
    ifelse( length( xPos ) == 1, 0, dfun * 2 * xVal^2 ) )
  vcovCoef <- diag( xCoefSE^2 )
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
ela4a <- probEla( c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
                 c( 1, 2.34, 3.3, 3.3^2, 0.0, 0.987 ),
                 xPos = 2 )
ela4a
# Example
ela4b <- probEla( c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
                 c( 1, 2.34, 3.3, 3.3^2, 0.0, 0.987 ),
                 c( 0.032, 0.004, 0.00001, 0.034, 0.0009, 0.056 ),
                 xPos = 2 )
ela4b
# Example
ela4c <- probEla( c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
                 c( 1, 2.34, 3.3, 3.3^2, 0.0, 0.987 ),
                 c( 0.032, 0.004, 0.00001, 0.034, 0.0009, 0.056 ),
                 xPos = c( 3, 4 ))
ela4c

## -----------------------------------------------------------------------------
probElaInt <- function( allCoef, allXVal, xPos, xBound,
                        allCoefSE = rep( NA, length( allCoef ) ) ){
  # number of coefficients
  nCoef <- length( allCoef )
  # number of intervals
  nInt <- length( xPos )
  # checking arguments
  if( length( allXVal ) != nCoef ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
  if( length( allCoefSE ) != nCoef ) {
    stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
  }
  checkXPos( xPos, minLength = 2, maxLength = nCoef,
    minVal = 0, maxVal = nCoef, requiredVal = 0 )
  if( any( allXVal[ xPos ] < 0 ) ) {
    stop( "all elements of argument 'allXVal'",
      " that are indicated by argument 'xPos'",
      " (i.e., the shares of observations in each interval)",
      " must be non-negative" )
  }
  if( sum( allXVal[ xPos ] > 1 ) ) {
    stop( "the sum of the elements of argument 'allXVal'",
      " that are indicated by argument 'xPos'",
      " (i.e., the shares of observations in each interval)",
      " must not be larger than one" )
  }
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # vector of probabilities of y=1 for each interval and
  # vector of shares of observations in each interval
  xBeta <- shareVec <- rep( NA, nInt )
  for( i in 1:nInt ){
    allXValTemp <- replace( allXVal, xPos, 0 )
    if( xPos[i] != 0 ) {
      allXValTemp[ xPos[i] ] <- 1
      shareVec[ i ] <- allXVal[ xPos[ i ] ]
    }
    xBeta[ i ] <- sum( allCoef * allXValTemp )
  }
  shareVec[ xPos == 0 ] <- 1 - sum( shareVec[ xPos != 0 ] )
  checkXBeta( xBeta )
  phiVec <- pnorm( xBeta )
  # weights
  weights <- elaIntWeights( shareVec )
  # calculation of the semi-elasticity
  semEla <- linElaInt( phiVec, shareVec, xBound )
  ### calculation of its standard error
  # partial derivatives of each semi-elasticity around each boundary
  # w.r.t. all estimated coefficients
  gradM <- matrix( 0, nCoef, nInt - 1 )
  gradPhiVec <- dnorm( xBeta )
  for( m in 1:( nInt - 1 ) ) {
    gradM[ -xPos, m ] <- 2 * ( gradPhiVec[m+1] - gradPhiVec[m] ) *
        allXVal[ -xPos ] * xBound[m+1] / ( xBound[m+2] - xBound[m] )
    gradM[ xPos[m], m ] <- - 2 * gradPhiVec[m] * xBound[m+1] /
      ( xBound[m+2] - xBound[m] )
    gradM[ xPos[m+1], m ] <- 2 * gradPhiVec[m+1] * xBound[m+1] /
      ( xBound[m+2] - xBound[m] )
  }
  # partial derivative of the semi-elasticity
  # w.r.t. all estimated coefficients
  derivCoef <- rep( 0, nCoef )
  for( m in 1:( nInt - 1 ) ){
    derivCoef <- derivCoef + weights[m] * gradM[ , m ]
  }
  # variance-covariance matrix of the coefficiencts
  vcovCoef <- diag( allCoefSE^2 )
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla[1], stdEr = semElaSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
ela5a <- probElaInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                    allXVal = c( 1, 0.4, 0.12, 0.13 ),
                    xPos = c( 2, 0, 3, 4 ),
                    xBound = c( 0, 500, 1000, 1500, Inf ))
ela5a
# Example
ela5b <- probElaInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                    allXVal = c( 1, 0.4, 0.12, 0.13 ),
                    xPos = c( 2, 0, 3, 4 ),
                    xBound = c( 0, 500, 1000, 1500, Inf ),
                    allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ))
ela5b

## -----------------------------------------------------------------------------
probEffInt <- function( allCoef, allXVal, xPos, refBound, intBound,
                        allCoefSE = rep( NA, length( allCoef ) ) ){
  # number of coefficients
  nCoef <- length( allCoef )
  # check arguments
  if( length( allXVal ) != nCoef ){
    stop( "argument 'allCoef' and 'allXVal' must have the same length" )
  }
  if( length( allCoefSE ) != nCoef ){
    stop( "argument 'allCoef' and 'allCoefSE' must have the same length" )
  }
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )
  if( any( !is.na( allXVal[ xPos ] ) ) ) {
    allXVal[ xPos ] <- NA
    warning( "values of argument 'allXVal[ xPos ]' are ignored",
      " (set these values to 'NA' to avoid this warning)" )
  }

  # calculate xBars
  intX <- mean( intBound )
  refX <- mean( refBound )
  if( length( xPos ) == 2 ) {
    intX <- c( intX, EXSquared( intBound[1], intBound[2] ) )
    refX <- c( refX, EXSquared( refBound[1], refBound[2] ) )
  }
  if( length( intX ) != length( xPos ) || length( refX ) != length( xPos ) ) {
    stop( "internal error: 'intX' or 'refX' does not have the expected length" )
  }
  # define X' * beta
  intXbeta <- sum( allCoef * replace( allXVal, xPos, intX ) )
  refXbeta <- sum( allCoef * replace( allXVal, xPos, refX ) )
  checkXBeta( c( intXbeta, refXbeta ) )
  # effect E_{k,ml}
  eff <- pnorm( intXbeta ) - pnorm( refXbeta )
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  derivCoef <- rep( NA, nCoef )
  derivCoef[ -xPos ] = ( dnorm( intXbeta ) - dnorm( refXbeta ) ) *
    allXVal[ -xPos ]
  derivCoef[ xPos ] = dnorm( intXbeta ) * intX - dnorm( refXbeta ) * refX
  # variance covariance of the coefficients (covariances set to zero)
  vcovCoef <- diag( allCoefSE^2 )
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # object to be returned
  result <- c( effect = eff, stdEr = effSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
eff4a <- probEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                     allXVal = c( 1, NA, 0.16, 0.13 ),
                     xPos = 2,
                     refBound = c( 8, 12 ), intBound = c( 13, 15 ))
eff4a

eff4b <- probEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                     allXVal = c( 1, NA, 0.16, 0.13 ),
                     xPos = 2,
                     refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                     allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ) )
eff4b

# Example
eff5a <- probEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.006 ),
                     allXVal = c( 1, NA, NA, 0.13 ),
                     xPos = c( 2, 3 ),
                     refBound = c( 8, 12 ), intBound = c( 13, 15 ))
eff5a

eff5b <- probEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.006 ),
                     allXVal = c( 1, NA, NA, 0.13 ),
                     xPos = c( 2, 3 ),
                     refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                     allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ) )
eff5b

## -----------------------------------------------------------------------------
probEffGr <- function( allCoef, allXVal, xPos, Group,
                       allCoefSE = rep( NA, length( allCoef ) ) ){
  nCoef <- length( allCoef )
  xShares <- allXVal[ xPos ]
  xCoef <- allCoef[ xPos ]
  if( sum( xShares ) > 1 ){
    stop( "the shares in argument 'xShares' sum up to a value larger than 1" )
  }
  if( length( xCoef ) != length( xShares ) ){
    stop( "arguments 'xCoef' and 'xShares' must have the same length" )
  }
  if( length( xCoef ) != length( Group ) ){
    stop( "arguments 'xCoef' and 'Group' must have the same length" )
  }
  if( ! all( Group %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'Group' must be -1, 0, or 1" )
  }
  # D_mr
  DRef <- sum( xCoef[ Group == -1 ] * xShares[ Group == -1 ]) /
          sum( xShares[ Group == -1 ] )
  XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DRef
  # D_ml
  DEffect <- sum( xCoef[ Group == 1 ] * xShares[ Group == 1 ]) /
             sum( xShares[ Group == 1 ] )
  XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DEffect
  # effect
  effeG <- pnorm( XBetaEffect ) - pnorm( XBetaRef )
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  derivCoef <- rep( NA, nCoef )
  derivCoef[ -xPos ] = ( dnorm( XBetaEffect ) - dnorm( XBetaRef ) ) *
                         allXVal[ -xPos ]
  derivCoef[ xPos ] = dnorm( XBetaEffect ) * DEffect - dnorm( XBetaRef ) * DRef
  # variance covariance of the coefficients (covariances set to zero)
  vcovCoef <- diag( allCoefSE^2 )
  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
Eff6a <- probEffGr( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034, -0.005, 0.89, -1.2 ),
                    allXVal = c( 1, 0.1, 0.3, 0.25, 0.15, 0.2, 2.34, 10.8 ),
                    xPos = c( 2:6 ), Group = c( 0, -1, -1, 1, 1 ) )
Eff6a
# Example
Eff6b <- probEffGr( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034, -0.005, 0.89, -1.2 ),
                    allXVal = c( 1, 0.1, 0.3, 0.25, 0.15, 0.2, 2.34, 10.8 ),
                    xPos = c( 2:6 ), Group = c( 0, -1, -1, 1, 1 ),
                    allCoefSE = c( 0.03, 0.0001, 0.005, 0, 0.01, 0.004, 0.05, 0.8 ))
Eff6b
# the same examples but without categories that are neither
# in the new reference category nor in the new category of interest
all.equal( Eff6a,
  probEffGr( allCoef = c( 0.28, 0.0075, 0, -0.034, -0.005, 0.89, -1.2 ),
    allXVal = c( 1, 0.3, 0.25, 0.15, 0.2, 2.34, 10.8 ),
    xPos = c( 2:5 ), Group = c( -1, -1, 1, 1 ) ) )
all.equal( Eff6b,
  probEffGr( allCoef = c( 0.28, 0.0075, 0, -0.034, -0.005, 0.89, -1.2 ),
    allXVal = c( 1, 0.3, 0.25, 0.15, 0.2, 2.34, 10.8 ),
    xPos = c( 2:5 ), Group = c( -1, -1, 1, 1 ),
    allCoefSE = c( 0.03, 0.005, 0, 0.01, 0.004, 0.05, 0.8 ) ) )

## -----------------------------------------------------------------------------
logEla <- function( allCoef, allCoefBra, allXVal, allXValBra,
                    allCoefSE = rep( NA, length( allCoef ) ),
                    xPos, yCat, yCatBra, lambda, method =  "binary" ){
  if( method == "binary" || method == "CondL" ){
    nCoef <- length( allCoef )
    # Checking standard errors
    if( nCoef != length( allCoefSE ) ) {
      stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
    }
  } else if( method == "MNL" ){
    NCoef <- length( allCoef )
    mCoef <- matrix( allCoef, nrow = length( allXVal ) )
    nCoef <- dim( mCoef )[1]
    pCoef <- dim( mCoef )[2]
    # Checking standard errors
    if( NCoef != length( allCoefSE ) ) {
      stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
    }
  } else{
    nCoef <- length( allCoef )
    NCoef <- length( allCoefBra )
    # Checking standard errors
    if( nCoef != length( allCoefSE ) ){
      stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
    }
  }
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # Check x values
  if( method == "binary" || method == "MNL" ){
    if( nCoef != length( allXVal ) ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  } else if( method == "CondL" ){
    mXVal <- matrix( allXVal, nrow = nCoef )
    nXVal <- dim( mXVal )[1]
    pXVal <- dim( mXVal )[2]
    if( nCoef != dim( mXVal )[1] ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  } else{
    mXValBra <- matrix( allXValBra, nrow = NCoef )
    nXValBra <- dim( mXValBra )[1]
    pXValBra <- dim( mXValBra )[2]
    if( NCoef != nXValBra ) {
      stop( "arguments 'allCoefBra' and 'allXValBra' must have the same length" )
    }
    O <- length( allXVal )
    mXVal <- matrix( unlist( allXVal[ yCatBra ] ), nrow = nCoef )
    nXVal <- dim( mXVal )[1]
    pXVal <- dim( mXVal )[2]
    if( nCoef != nXVal ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  }
  # Identify coefficients of interest (kth/tth covariate)
  if( length( xPos ) == 2 ){
    if( method == "binary" ){
      xCoef <- allCoef[ xPos ]
      if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
        stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
          "to the squared value of 'allXVal[ xPos[1] ]' " )
      }
    } else if( method == "MNL" ){
      xCoef <- mCoef[ xPos, ]
      if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
        stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
          "to the squared value of 'allXVal[ xPos[1] ]' " )
      }
    } else if( method == "CondL" ){
      xCoef <- allCoef[ xPos ]
      for( p in 1:pXVal ){
        if( !isTRUE( all.equal( mXVal[xPos[2], p], mXVal[xPos[1], p]^2 ) ) ) {
          stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
            "to the squared value of 'allXVal[ xPos[1] ]' " )
        }
      }
    } else{
      xCoef <- allCoef[ xPos ]
      for( p in 1:pXVal ){
        if( !isTRUE( all.equal( mXVal[xPos[2], p], mXVal[xPos[1], p]^2 ) ) ) {
          stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
            "to the squared value of 'allXVal[ xPos[1] ]' " )
        }
      }
    }
  } else if( length( xPos ) == 1 ) {
    if( method == "binary" || method == "CondL" ){
      xCoef <- c( allCoef[ xPos ], 0 )
    } else if( method == "MNL" ){
      xCoef <- matrix( c( mCoef[ xPos, ], rep( 0, dim( mCoef )[ 2 ] ) ),
        nrow = 2, byrow = TRUE  )
    } else{
      xCoef <- c( allCoef[ xPos ], 0 )
    }
  } else {
    stop( "argument 'xPos' must be a scalar or a vector with two elements" )
  }
  if( method == "binary" ){
    xVal <- allXVal[ xPos[1] ]
    xBeta <- sum( allCoef * allXVal )
    checkXBeta( xBeta )
    dfun <- exp( xBeta )/( 1 + exp( xBeta ) )^2
    semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  } else if( method == "MNL" ){     #checkXBeta missing
    xVal <- allXVal[ xPos[1] ]
    xBeta <- colSums( mCoef * allXVal )
    pfun <- rep( NA, length( xBeta ))
    term <- 0
    for( i in 1:length( xBeta )){
      pfun[i] <- exp( xBeta[i] )/( 1 + sum( exp( xBeta ) ) )
      term <- term + ( ( xCoef[ 1, yCat ] + 2 * xCoef[ 2, yCat ] * xVal ) -
          ( xCoef[ 1, i ] + 2 * xCoef[ 2, i ] * xVal ) * pfun[i] )
    }
    semEla <- xVal * pfun[ yCat ] *
      ( ( xCoef[ 1, yCat ] + 2 * xCoef[ 2, yCat ] * xVal )/
          ( 1 + sum( exp( xBeta ))) + term )
    dfun <- pfun[ yCat ] * ( 1/( 1 + sum( exp( xBeta ) ) ) + term )
  } else if( method == "CondL" ){    #checkXBeta missing
    xVal <- rep( NA, pXVal )
    for( p in 1:pXVal ){
      xVal[p] <- mXVal[ xPos[ 1 ], p ]
    }
    xBeta <- colSums( allCoef * mXVal )
    pfun <- exp( xBeta[ yCat ] )/( sum( exp( xBeta ) ) )
    semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal[ yCat ] ) *
      xVal[ yCat ] * ( pfun - pfun^2 )
  } else{                            #checkXBeta missing
    xVal <- rep( NA, pXVal )
    for( p in 1:pXVal ){
      xVal[p] <- mXVal[ xPos[ 1 ], p ]
    }
    coef <- matrix( NA, nrow = O, ncol = nCoef )
    for( o in 1:O ){
      coef[o, ] <- allCoef/lambda[o]
    }
    xBeta <- lapply( 1:O, function( i, m, v ){ colSums( m[[i]] * v[[i]] ) },
      m=allXVal, v=coef )
    IV <- unlist( lapply( 1:O, function( i, m ){ log( sum( exp( m[[i]] ) ) ) },
      m=xBeta ) )
    pfun <- exp( xBeta[[ yCatBra ]][ yCat ] )/
      ( sum( exp( xBeta[[ yCatBra ]] ) ) )
    xBetaBra <- colSums( allCoefBra * mXValBra )
    pfunBra <- exp( xBetaBra[ yCatBra ] + lambda[ yCatBra ] * IV[ yCatBra ] )/
      ( sum( exp( xBetaBra + lambda * IV ) ) )
    semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal[ yCat ] ) * xVal[ yCat ] *
      ( pfunBra * ( pfun - pfun^2 ) * 1/lambda[ yCatBra ] +
          pfun^2 * ( pfunBra - pfunBra^2 ) * lambda[ yCatBra ] * IV[ yCatBra ] )
  }
  if( method == "binary" || method == "MNL" ){
    derivCoef <- rep( 0, length( allCoef ) )
    derivCoef[ xPos[1] ] <- dfun * xVal
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- dfun * 2 * xVal^2
    }
  } else if( method == "CondL" ){
    derivCoef <- rep( 0, length( allCoef ) )
    derivCoef[ xPos[1] ] <- ( pfun - pfun^2 ) * xVal[ yCat ]
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- ( pfun - pfun^2 ) * 2 * xVal[ yCat ]^2
    }
  } else{
    derivCoef <- rep( 0, length( allCoef ) )
    derivCoef[ xPos[1] ] <- (
      pfunBra * ( pfun - pfun^2 ) / lambda[ yCatBra ] +
        pfun^2 * ( pfunBra - pfunBra^2 ) * lambda[ yCatBra ] *
        IV[ yCatBra ] ) *
      xVal[ yCat ]
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- (
        pfunBra * ( pfun - pfun^2 ) / lambda[ yCatBra ] +
          pfun^2 * ( pfunBra - pfunBra^2 ) * lambda[ yCatBra ] *
          IV[ yCatBra ] ) *
        2 * xVal[ yCat ]^2
    }
  }
  vcovCoef <- diag( allCoefSE^2 )
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
ela6a <- logEla( allCoef = c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
                 allXVal = c( 1, 3.3, 4.5, 2.34, 0.1, 0.987 ), xPos = 2 )
ela6a

ela6b <- logEla( allCoef = c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
                 allXVal = c( 1, 3.3, 4.5, 2.24, 0.1, 0.987 ),
                 allCoefSE = c( 0.001, 0.02, 0.000002, 0.05, 1.2, 0.03 ),
                 xPos = 2 )
ela6b

# Example
ela7a <- logEla( allCoef = c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
                 allXVal = c( 1, 3.3, 3.3^2, 2.34, 0.1, 0.987 ),
                 xPos = c( 2, 3 ) )
ela7a

ela7b <- logEla( allCoef = c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
                 allXVal = c( 1, 3.3, 3.3^2, 2.34, 0.1, 0.987 ),
                 allCoefSE = c( 0.001, 0.02, 0.000002, 0.05, 1.2, 0.03 ),
                 xPos = c( 2, 3 ) )
ela7b

# Example
ela8a <- logEla( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                 allXVal = c( 1, 8.4, 0.06 ), xPos = 3,
                 method = "MNL", yCat = 2 )
ela8a

ela8b <- logEla( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                 allXVal = c( 1, 8.4, 0.06 ),
                 allCoefSE = c( 0.002, 0.003, 0.004, 0.006, 0.00001, 0.08 ),
                 xPos = 3,
                 method = "MNL", yCat = 2 )
ela8b

# Example
ela9a <- logEla( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                 allXVal = c( 1, 0.04, 0.0016 ), xPos = c( 2, 3 ),
                 method = "MNL", yCat = 2 )
ela9a

ela9b <- logEla( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                 allXVal = c( 1, 0.04, 0.0016 ),
                 allCoefSE = c( 0.002, 0.003, 0.004, 0.006, 0.00001, 0.08 ),
                 xPos = c( 2, 3 ),
                 method = "MNL", yCat = 2 )
ela9b

# Example
ela10a <- logEla( allCoef = c( 0.445, 0.03, 0.00002 ),
                  allXVal = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
                  xPos = 2,
                  method = "CondL", yCat = 2 )
ela10a

ela10b <- logEla( allCoef = c( 0.445, 0.03, -0.002 ),
                  allXVal = c( 1, 0.3, 0.09, 1, 0.1, 0.01 ),
                  xPos = c( 2, 3 ),
                  method = "CondL", yCat = 2 )
ela10b

# Example
ela11a <- logEla( allCoef = c( 0.445, 0.03, 0.00002 ),
                  allXVal = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
                  allCoefSE = c( 0.002, 0.003, 0.004 ),
                  xPos = 2,
                  method = "CondL", yCat = 2 )
ela11a

ela11b <- logEla( allCoef = c( 0.445, 0.03, -0.002 ),
                  allXVal = c( 1, 0.3, 0.09, 1, 0.1, 0.01 ),
                  allCoefSE = c( 0.002, 0.003, 0.004 ),
                  xPos = c( 2, 3 ),
                  method = "CondL", yCat = 2 )
ela11b

# Example
matrix1 <- matrix( c( 1, 2.5, 0.3, 0.09, 1, 0.33, 0.9, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, 2.8, 0.099, 0.211 ), nrow = 4 )
ela12a <- logEla( allCoefBra = c( 0.445, 0.03, -0.002 ),
                  allCoef = c( 1.8, 0.005, -0.12, 0.8 ),
                  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
                  allXVal = list( matrix1, matrix2 ),
                  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ),
                  method = "NestedL" )
ela12a

matrix1 <- matrix( c( 1, 0.3, 0.09, 0.09, 1, 0.33, 0.1089, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, 0.31, 0.099, 0.211 ), nrow = 4 )
ela12b <- logEla( allCoefBra = c( 0.445, 0.03, -0.002 ),
                  allCoef = c( 1.8, 0.005, -0.12, 0.8 ),
                  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
                  allXVal = list( matrix1, matrix2 ),
                  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ),
                  method = "NestedL" )
ela12b

# Example
matrix1 <- matrix( c( 1, 2.5, 0.3, 0.09, 1, 0.33, 0.9, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, 2.8, 0.099, 0.211 ), nrow = 4 )
ela13a <- logEla( allCoefBra = c( 0.445, 0.03, -0.002 ),
                  allCoef = c( 1.8, 0.005, -0.12, 0.8 ),
                  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
                  allXVal = list( matrix1, matrix2 ),
                  allCoefSE = c( 0.001, 0.089, 0.0003, 0.12 ),
                  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ),
                  method = "NestedL" )
ela13a

matrix1 <- matrix( c( 1, 0.3, 0.09, 0.09, 1, 0.33, 0.1089, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, 0.31, 0.099, 0.211 ), nrow = 4 )
ela13b <- logEla( allCoefBra = c( 0.445, 0.03, -0.002 ),
                  allCoef = c( 1.8, 0.005, -0.12, 0.8 ),
                  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
                  allXVal = list( matrix1, matrix2 ),
                  allCoefSE = c( 0.001, 0.089, 0.0003, 0.12 ),
                  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ),
                  method = "NestedL" )
ela13b

## -----------------------------------------------------------------------------
logElaInt <- function( allCoef, allXVal, xPos, xBound, yCat = NA,
                       allCoefSE = rep( NA, length( allCoef ) ),
                       method = "binary" ){
  # number of coefficients
  if( method == "binary" || method == "MNL" ){
    mCoef <- matrix( allCoef, nrow = length( allXVal ))
    nCoef <- dim( mCoef )[1]
    pCoef <- dim( mCoef )[2]
    # checking arguments
    if( length( allXVal ) != nCoef ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  } else{
    nCoef <- length( allCoef )
    mXVal <- matrix( allXVal, nrow = nCoef )
    pCoef <- dim( mXVal )[2]
    # checking arguments
    if( dim( mXVal )[1] != nCoef ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  }
  # number of intervals
  nInt <- length( xPos )
  checkXPos( xPos, minLength = 2, maxLength = nCoef,
             minVal = 0, maxVal = nCoef, requiredVal = 0 )
  if( method == "binary" || method == "MNL" ){
    if( any( allXVal[ xPos ] < 0 ) ) {
        stop( "all elements of argument 'allXVal'",
              " that are indicated by argument 'xPos'",
              " (i.e., the shares of observations in each interval)",
              " must be non-negative" )
      }
    if( sum( allXVal[ xPos ] > 1 ) ) {
        stop( "the sum of the elements of argument 'allXVal'",
              " that are indicated by argument 'xPos'",
              " (i.e., the shares of observations in each interval)",
              " must not be larger than one" )
    }
  } else{
    for( p in 1:pCoef ){
    if( any( mXVal[ xPos, p ] < 0 ) ) {
        stop( "all elements of argument 'allXVal'",
              " that are indicated by argument 'xPos'",
              " (i.e., the shares of observations in each interval)",
              " must be non-negative" )
      }
    if( sum( mXVal[ xPos, p ] > 1 ) ) {
        stop( "the sum of the elements of argument 'allXVal'",
              " that are indicated by argument 'xPos'",
              " (i.e., the shares of observations in each interval)",
              " must not be larger than one" )
      }
    }
  }
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # vector of probabilities of y=1 for each interval and
  # vector of shares of observations in each interval
  xBeta <- matrix( rep( rep( NA, nInt ), pCoef ), ncol = pCoef )
  if( method == "binary" || method == "MNL" ){
    shareVec <- rep( NA, nInt )
    for( p in 1:pCoef ){
          for( i in 1:nInt ){
              allXValTemp <- replace( allXVal, xPos, 0 )
              if( xPos[i] != 0 ) {
                allXValTemp[ xPos[i] ] <- 1
                shareVec[i] <- allXVal[ xPos[i] ]
              }
            xBeta[i,p] <- sum( mCoef[ ,p] * allXValTemp )
          }
    }
    shareVec[ xPos == 0 ] <- 1 - sum( shareVec[ xPos != 0 ] )
  } else{
    shareVec <- matrix( rep( rep( NA, nInt ), pCoef ), ncol = pCoef )
    for( p in 1:pCoef ){
         for( i in 1:nInt ){
             allXValTemp <- replace( mXVal[ ,p], xPos, 0 )
             if( xPos[i] != 0 ) {
                allXValTemp[ xPos[i] ] <- 1
                shareVec[i,p] <- mXVal[ xPos[i], p ]
             }
           xBeta[i,p] <- sum( allCoef * allXValTemp )
         }
      shareVec[ xPos == 0, p ] <- 1 - sum( shareVec[ xPos != 0, p ] )
    }
    shareVec <- shareVec[ , yCat ]
  }
  #checkXBeta( xBeta )  #Please check this one with a matrix
  if( method == "binary" ){
     expVec <- as.vector( exp( xBeta )/( 1 + exp( xBeta ) ) )
  } else if( method == "MNL" ){
     expVec <- as.vector( exp( xBeta[ , yCat ])/( 1 + rowSums( exp( xBeta ) ) ) )
  } else{
     expVec <- as.vector( exp( xBeta[ , yCat ])/( rowSums( exp( xBeta ) ) ) )
  }
  # weights
  weights <- elaIntWeights( shareVec )
  # calculation of the semi-elasticity
  semEla <- linElaInt( expVec, shareVec, xBound )
  ### calculation of its standard error
  # partial derivatives of each semi-elasticity around each boundary
  # w.r.t. all estimated coefficients
  if( method == "binary" ){
    gradM <- matrix( 0, nCoef, nInt - 1 )
    gradExpVec <- exp( xBeta )/( 1 + exp( xBeta ) )^2
    for( m in 1:( nInt - 1 ) ) {
      gradM[ -xPos, m ] <- 2 * ( gradExpVec[m+1] - gradExpVec[m] ) *
          allXVal[ -xPos ] * xBound[m+1] / ( xBound[m+2] - xBound[m] )
      gradM[ xPos[m], m ] <- - 2 * gradExpVec[m] * xBound[m+1] /
        ( xBound[m+2] - xBound[m] )
      gradM[ xPos[m+1], m ] <- 2 * gradExpVec[m+1] * xBound[m+1] /
        ( xBound[m+2] - xBound[m] )
    }
  } else if( method == "MNL" ){
    gradM <- array( 0, c( nCoef, nInt - 1, pCoef ) )
    gradExpVecP <- ( exp( xBeta[ , yCat ] ) *
         ( 1 + rowSums( exp( xBeta[ , -yCat, drop = FALSE ] ) ) ) )/
         ( 1 + rowSums( exp( xBeta ) ) )^2
    for( p in 1:pCoef ){
      gradExpVecO <- ( exp( xBeta[ , yCat ] ) * exp( xBeta[ , p] ) )/
         ( 1 + rowSums( exp( xBeta ) ) )^2
      for( m in 1:( nInt - 1 ) ) {
        if( p == yCat ){
          gradM[ -xPos, m, p ] <- 2 * ( gradExpVecP[m+1] - gradExpVecP[m] ) *
               allXVal[ -xPos ] * xBound[m+1] / ( xBound[m+2] - xBound[m] )
          gradM[ xPos[ m ], m, p ] <- - 2 * gradExpVecP[m] * xBound[m+1] /
                                                ( xBound[m+2] - xBound[m] )
          gradM[ xPos[ m + 1 ], m, p ] <- 2 * gradExpVecP[m+1] * xBound[m+1] /
                                                ( xBound[m+2] - xBound[m] )
        } else {
          gradM[ -xPos, m, p ] <- 2 * ( gradExpVecO[m] - gradExpVecO[m+1] ) *
               allXVal[ -xPos ] * xBound[m+1] / ( xBound[m+2] - xBound[m] )
          gradM[ xPos[ m ], m, p ] <- 2 * gradExpVecO[m] * xBound[m+1] /
                                                ( xBound[m+2] - xBound[m] )
          gradM[ xPos[ m + 1 ], m, p ] <- - 2 * gradExpVecO[m+1] * xBound[m+1] /
                                                ( xBound[m+2] - xBound[m] )
        }
      }
    }
    gradM <- apply( gradM, 2, function( x ) x )
  } else{
    gradM <- matrix( 0, nCoef, nInt - 1 )
    for( m in 1:( nInt - 1 ) ) {
      gradM[ -xPos, m ] <- 2 *
          ( ( exp( xBeta[ m+1, yCat ] ) * mXVal[ -xPos, yCat ] *
              sum( exp( xBeta[ m+1, ] ) ) -
              exp( xBeta[ m+1, yCat ] ) *
              rowSums( exp( xBeta[ m+1, ] ) * mXVal[ -xPos, , drop = FALSE ] ) )/
            ( sum( exp( xBeta[ m+1, ] ) ) )^2 -
            ( exp( xBeta[ m, yCat ] ) * mXVal[ -xPos, yCat ] *
              sum( exp( xBeta[ m, ] ) ) -
              exp( xBeta[ m, yCat ] ) *
              rowSums( exp( xBeta[ m,  ] ) * mXVal[ -xPos, , drop = FALSE ] ) )/
            ( sum( exp( xBeta[ m, ] ) ) )^2 ) *
              xBound[m+1] / ( xBound[m+2] - xBound[m] )
      gradM[ xPos[m], m ] <- 0
      gradM[ xPos[m+1], m ] <- 0
    }
  }
  # partial derivative of the semi-elasticity
  # w.r.t. all estimated coefficients
  derivCoef <- rep( 0, length( allCoef ) )
  for( m in 1:( nInt - 1 ) ){
    derivCoef <- derivCoef + weights[m] * gradM[,m]
  }
  # variance-covariance matrix of the coefficiencts
  vcovCoef <- diag( allCoefSE^2 )
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla[1], stdEr = semElaSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
ela8a <- logElaInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                    allXVal = c( 1, 0.4, 0.12, 0.13 ),
                    xPos = c( 2, 0, 3, 4 ),
                    xBound = c( 0, 500, 1000, 1500, Inf ) )
ela8a

ela8b <- logElaInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                    allXVal = c( 1, 0.4, 0.12, 0.13 ),
                    xPos = c( 2, 0, 3, 4 ),
                    xBound = c( 0, 500, 1000, 1500, Inf ),
                    allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ) )
ela8b

# Example
ela9a <- logElaInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                    allXVal = c( 1, 0.4, 0.12 ),
                    xPos = c( 2, 0, 3 ),
                    xBound = c( 0, 500, 1000, Inf ), yCat = 2,
                    method = "MNL" )
ela9a

ela9b <- logElaInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                    allXVal = c( 1, 0.4, 0.12 ),
                    xPos = c( 2, 0, 3 ),
                    xBound = c( 0, 500, 1000, Inf ), yCat = 2,
                    allCoefSE = c( 0.003, 0.045, 0.007, 0.009, 0.0008, 0.9 ),
                    method = "MNL" )
ela9b

# Example
ela10a <- logElaInt( allCoef = c( 1.33, 0.022, 0.58, 1.6 ),
                     allXVal = c( 1, 0.4, 0.12, 0.0002,
                                  1, 0.28, 0.01, 0.000013 ),
                     xPos = c( 2, 0, 3 ),
                     xBound = c( 0, 1000, 1500, Inf ), yCat = 2,
                     method = "CondL" )
ela10a

ela10b <- logElaInt( allCoef = c( 1.33, 0.022, 0.58, 1.6 ),
                     allXVal = c( 1, 0.4, 0.12, 0.0002,
                                  1, 0.28, 0.01, 0.000013 ),
                     xPos = c( 2, 0, 3 ),
                     xBound = c( 0, 1000, 1500, Inf ),
                     allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ), yCat = 2,
                     method = "CondL" )
ela10b

## -----------------------------------------------------------------------------
logEffInt <- function( allCoef, allCoefBra = NA, allXVal, allXValBra=NA,
                       xPos, refBound, intBound, yCat, yCatBra, lambda,
                       allCoefSE = rep( NA, length( allCoef ) ),
                       method = "binary" ){
if( method == "binary" ){
  # number of coefficients
  nCoef <- length( allCoef )
  # check arguments
  if( length( allXVal ) != nCoef ){
    stop( "argument 'allCoef' and 'allXVal' must have the same length" )
  }
  if( length( allCoefSE ) != nCoef ){
    stop( "argument 'allCoef' and 'allCoefSE' must have the same length" )
  }
} else if( method == "MNL" ){
  # number of coefficients
    NCoef <- length( allCoef )
    mCoef <- matrix( allCoef, nrow = length( allXVal ))
    nCoef <- dim( mCoef )[1]
    pCoef <- dim( mCoef )[2]
   # check arguments
    if( length( allXVal ) != nCoef ){
      stop( "argument 'allCoef' and 'allXVal' must have the same length" )
    }
    if( length( allCoefSE ) != NCoef ){
      stop( "argument 'allCoef' and 'allCoefSE' must have the same length" )
    }
} else if( method == "CondL"){
  # number of coefficients
    nCoef <- length( allCoef )
    mXVal <- matrix( allXVal, nrow = nCoef )
    pCoef <- dim( mXVal )[2]
  # check arguments
    if( dim( mXVal )[1] != nCoef ){
      stop( "argument 'allCoef' and 'allXVal' must have the same length" )
    }
    if( length( allCoefSE ) != nCoef ){
      stop( "argument 'allCoef' and 'allCoefSE' must have the same length" )
    }
} else{
   nCoef <- length( allCoef )
   NCoef <- length( allCoefBra )
   mXValBra <- matrix( allXValBra, nrow = NCoef )
   nXValBra <- dim( mXValBra )[1]
   pXValBra <- dim( mXValBra )[2]
 # check arguments
   if( NCoef != nXValBra ){
     stop( "arguments 'allCoefBra' and 'allXValBra' must have the same length")
   }
   O <- length( allXVal )
   nXVal <- unlist( lapply( allXVal, function(x) dim( x )[1] ) )
   pCoef <- unlist( lapply( allXVal, function(x) dim( x )[2] ) )
   if( nCoef != nXVal[ yCatBra ] ){
     stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
   }
   if( nCoef != length( allCoefSE) ){
     stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
   }
}
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )
if( method == "binary" || method == "MNL" ){
  if( any( !is.na( allXVal[ xPos ] ) ) ) {
    allXVal[ xPos ] <- NA
    warning( "values of argument 'allXVal[ xPos ]' are ignored",
      " (set these values to 'NA' to avoid this warning)" )
  }
}else if( method == "CondL" ){
  for( p in 1:pCoef ){
    if( any( !is.na( mXVal[ xPos, p ] ) ) ){
      mXVal[ xPos, p ] <- NA
      warning( "values of argument 'allXVal[ xPos ]' are ignored",
      " (set these values to 'NA' to avoid this warning)" )
    }
  }
}else{
  for( p in 1:pCoef[ yCatBra ] ){
    if( any( !is.na( allXVal[[ yCatBra ]][ xPos, p ] ) ) ){
      mXVal[ xPos, p ] <- NA
      warning( "values of argument 'allXVal[ xPos ]' are ignored",
      " (set these values to 'NA' to avoid this warning)" )
    }
  }
}
  # calculate xBars
  intX <- mean( intBound )
  refX <- mean( refBound )
  if( length( xPos ) == 2 ) {
    intX <- c( intX, EXSquared( intBound[1], intBound[2] ) )
    refX <- c( refX, EXSquared( refBound[1], refBound[2] ) )
  }
  if( length( intX ) != length( xPos ) || length( refX ) != length( xPos ) ) {
    stop( "internal error: 'intX' or 'refX' does not have the expected length" )
  }
  # define X' * beta
  if( method == "binary" ){
    intXbeta <- sum( allCoef * replace( allXVal, xPos, intX ) )
    refXbeta <- sum( allCoef * replace( allXVal, xPos, refX ) )
    checkXBeta( c( intXbeta, refXbeta ) )
  } else if( method == "MNL" ){
    intXbeta <- colSums( mCoef * replace( allXVal, xPos, intX ) )
    refXbeta <- colSums( mCoef * replace( allXVal, xPos, refX ) )
  } else if( method == "CondL" ){
    mXValint <- mXValref <- mXVal
    for( p in 1:pCoef ){
      mXValint[ ,p] <- replace( mXValint[ ,p], xPos, intX )
      mXValref[ ,p] <- replace( mXValref[ ,p], xPos, refX )
    }
    intXbeta <- colSums( allCoef * mXValint )
    refXbeta <- colSums( allCoef * mXValref )
  } else{
    mCoef <- matrix( rep( allCoef, O ), nrow = nCoef, O ) %*% diag( 1/ lambda )
    mXValint <- mXValref <- allXVal
    for( i in 1:O ){
      for( p in 1:pCoef[i] ){
        mXValint[[i]][ ,p] <- replace( mXValint[[i]][ ,p], xPos, intX )
        mXValref[[i]][ ,p] <- replace( mXValref[[i]][ ,p], xPos, refX )
      }
    }
    refXbeta <- intXbeta <- rep( list( NA ), O )
    for( l in 1:O ){
      intXbeta[[ l ]] <- colSums( mCoef[ ,l ] * mXValint[[ l ]] )
      refXbeta[[ l ]] <- colSums( mCoef[ ,l ] * mXValref[[ l ]] )
    }
    XbetaBra <- colSums( allCoefBra * mXValBra )
  }
  # effect E_{k,ml}
  if( method == "binary" ){
    eff <- exp( intXbeta )/( 1 + exp( intXbeta ) ) -
           exp( refXbeta )/( 1 + exp( refXbeta ) )
  } else if( method == "MNL" ){
    eff <- exp( intXbeta[ yCat ] )/( 1 + sum( exp( intXbeta ) ) ) -
           exp( refXbeta[ yCat ] )/( 1 + sum( exp( refXbeta ) ) )
  } else if( method == "CondL"){
    eff <- exp( intXbeta[ yCat ] )/( sum( exp( intXbeta ) ) ) -
           exp( refXbeta[ yCat ] )/( sum( exp( refXbeta ) ) )
  } else{
    intBranch <- refBranch <- rep( list( NA ), O )
    for( l in 1:O ){
      intBranch[[ l ]] <- exp( XbetaBra[ l ] + lambda[ l ] *
                          log( sum( exp( intXbeta[[ l ]] ) ) ) )
      refBranch[[ l ]] <- exp( XbetaBra[ l ] + lambda[ l ] *
                          log( sum( exp( refXbeta[[ l ]] ) ) ) )
    }
    intBranch <- unlist( intBranch )
    refBranch <- unlist( refBranch )
    eff <- exp( intXbeta[[ yCatBra ]][ yCat ] )/( sum( exp( intXbeta[[ yCatBra ]] ) ) ) *
           intBranch[ yCatBra ]/ sum( intBranch ) -
           exp( refXbeta[[ yCatBra ]][ yCat ] )/( sum( exp( refXbeta[[ yCatBra ]] ) ) ) *
           refBranch[ yCatBra ]/ sum( refBranch )
  }
  # calculating approximate standard error
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  if( method == "binary" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( intXbeta )/( 1 + exp( intXbeta ) )^2 -
                            exp( refXbeta )/( 1 + exp( refXbeta ) )^2 ) *
                            allXVal[ -xPos ]
    derivCoef[ xPos ] <- exp( intXbeta )/( 1 + exp( intXbeta ) )^2 * intX -
                         exp( refXbeta )/( 1 + exp( refXbeta ) )^2 * refX
  } else if( method == "MNL" ){
    derivCoef <- matrix( NA, nrow=nCoef, ncol=pCoef )
    for( p in 1:pCoef ){
      if( p == yCat ){
        derivCoef[ -xPos, p ] <-
          ( exp( intXbeta[ p ] ) *
          ( 1 + sum( exp( intXbeta[ -yCat ] ) ) )/
          ( 1 + sum( exp( intXbeta ) ) )^2 -
            exp( refXbeta[ p ] ) *
          ( 1 + sum( exp( refXbeta[ -yCat ] ) ) )/
          ( 1 + sum( exp( refXbeta ) ) )^2 ) * allXVal[ - xPos ]
        derivCoef[ xPos, p ] <-
          ( exp( intXbeta[ p ] ) *
          ( 1 + sum( exp( intXbeta[ -yCat ] ) ) )/
          ( 1 + sum( exp( intXbeta ) ) )^2 ) * intX -
          ( exp( refXbeta[ p ] ) *
          ( 1 + sum( exp( refXbeta[ -yCat ] ) ) )/
          ( 1 + sum( exp( refXbeta ) ) )^2 ) * refX
      } else{
        derivCoef[ -xPos, p ] <-
          ( ( exp( refXbeta[ yCat ] ) * exp( refXbeta[ p ] ) )/
            ( 1 + sum( exp( refXbeta ) ) )^2 -
            ( exp( intXbeta[ yCat ] ) * exp( intXbeta[ p ] ) )/
            ( 1 + sum( exp( intXbeta ) ) )^2 ) * allXVal[ -xPos ]
        derivCoef[ xPos, p ] <-
            ( ( exp( refXbeta[ yCat ] ) * exp( refXbeta[ p ] ) )/
              ( 1 + sum( exp( refXbeta ) ) )^2 ) * intX -
            ( ( exp( intXbeta[ yCat ] ) * exp( intXbeta[ p ] ) )/
              ( 1 + sum( exp( intXbeta ) ) )^2 ) * refX
      }
    }
    derivCoef <- c( derivCoef )
  } else if( method == "CondL" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( intXbeta[ yCat] ) * mXVal[ -xPos, yCat] *
                            sum( exp( intXbeta ) ) -
                            exp( intXbeta[ yCat] ) * rowSums( exp( intXbeta ) *
                            mXVal[ -xPos, ] ) )/
                          ( sum( exp( intXbeta ) ) )^2 -
                          ( exp( refXbeta[ yCat] ) * mXVal[ -xPos, yCat] *
                            sum( exp( refXbeta ) ) -
                            exp( refXbeta[ yCat] ) * rowSums( exp( refXbeta ) *
                            mXVal[ -xPos, ] ) )/
                          ( sum( exp( refXbeta ) ) )^2
    derivCoef[ xPos ] <-  ( exp( intXbeta[ yCat] ) * intX *
                            sum( exp( intXbeta ) ) -
                            exp( intXbeta[ yCat] ) * sum( exp( intXbeta ) * intX ) )/
                          ( sum( exp( intXbeta ) ) )^2 -
                          ( exp( refXbeta[ yCat] ) * refX *
                            sum( exp( refXbeta ) ) -
                            exp( refXbeta[ yCat] ) * sum( exp( refXbeta ) * refX ) )/
                          ( sum( exp( refXbeta ) ) )^2
  } else{
    derivCoef <- rep( NA, nCoef )
    PImp <- exp( intXbeta[[ yCatBra ]][ yCat ])/( sum( exp( intXbeta[[ yCatBra ]] ) ) )
    PIlp <- exp( refXbeta[[ yCatBra ]][ yCat ])/( sum( exp( refXbeta[[ yCatBra ]] ) ) )
    PImo <- intBranch[ yCatBra ]/ sum( intBranch )
    PIlo <- refBranch[ yCatBra ]/ sum( refBranch )
    Om <- matrix(
          unlist( lapply( allXVal, function(x) rowSums( x[ -xPos, , drop = FALSE ] ) ) ),
          ncol = O )
    derivCoef[ -xPos ] <- ( ( allXVal[[ yCatBra ]][ -xPos, yCat ]/lambda[ yCatBra ] -
                          ( rowSums(
                            ( allXVal[[ yCatBra ]][ -xPos, ]/lambda[ yCatBra ] ) %*%
                              diag( exp( intXbeta[[ yCatBra ]] ) ) ) )/
                          ( sum( exp( intXbeta[[ yCatBra ]] ) ) ) ) +
                          ( rowSums( allXVal[[ yCatBra ]][ -xPos, ] ) -
                          ( rowSums( Om %*% diag( exp( intBranch ) ) )/
                          ( sum( intBranch ) ) ) ) ) * PImp * PImo -
                          ( ( allXVal[[ yCatBra ]][ -xPos, yCat ]/lambda[ yCatBra ] -
                          ( rowSums(
                            ( allXVal[[ yCatBra ]][ -xPos, ]/lambda[ yCatBra ] ) %*%
                              diag( exp( refXbeta[[ yCatBra ]] ) ) ) )/
                          ( sum( exp( refXbeta[[ yCatBra ]] ) ) ) ) +
                          ( rowSums( allXVal[[ yCatBra ]][ -xPos, ] ) -
                          ( rowSums( Om %*% diag( exp( refBranch ) ) )/
                          ( sum( refBranch ) ) ) ) ) * PIlp * PIlo
    derivCoef[ xPos ] <-  ( ( intX/lambda[ yCatBra ] -
                          ( sum( intX/lambda[ yCatBra ]  *
                                 exp( intXbeta[[ yCatBra ]] ) ) )/
                          ( sum( exp( intXbeta[[ yCatBra ]] ) ) ) ) +
                          ( intX * pCoef[ yCatBra ] -
                          ( sum( intX * exp( intBranch ) )/
                          ( sum( intBranch ) ) ) ) ) * PImp * PImo -
                          ( ( refX/lambda[ yCatBra ] -
                          ( sum( refX/lambda[ yCatBra ]  *
                                 exp( refXbeta[[ yCatBra ]] ) ) )/
                          ( sum( exp( refXbeta[[ yCatBra ]] ) ) ) ) +
                          ( refX * pCoef[ yCatBra ] -
                          ( sum( refX * exp( refBranch ) )/
                          ( sum( refBranch ) ) ) ) ) * PImp * PImo
  }
  # variance covariance of the coefficients (covariances set to zero)
  vcovCoef <- diag( allCoefSE^2 )
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # object to be returned
  result <- c( effect = eff, stdEr = effSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
eff6a <- logEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                    allXVal = c( 1, NA, 0.16, 0.13 ),
                    xPos = 2,
                    refBound = c( 8, 12 ), intBound = c( 13, 15 ) )
eff6a

eff6b <- logEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                    allXVal = c( 1, NA, 0.16, 0.13 ),
                    xPos = 2,
                    refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                    allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ) )
eff6b

# Example
eff7a <- logEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                    allXVal = c( 1, NA, NA, 0.0004 ),
                    xPos = c( 2, 3 ),
                    refBound = c( 8, 12 ), intBound = c( 13, 15 ))
eff7a

eff7b <- logEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
                    allXVal = c( 1, NA, NA, 0.13 ),
                    xPos = c( 2, 3 ),
                    refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                    allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ) )
eff7b

#Example
eff8a <- logEffInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                    allXVal = c( 1, NA, 0.12 ),
                    xPos = 2,
                    refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                    yCat = 2, method = "MNL" )
eff8a

eff8b <- logEffInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                    allXVal = c( 1, NA, 0.12 ),
                    xPos = 2,
                    refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                    yCat = 2,
                    allCoefSE = c( 0.003, 0.045, 0.007, 0.009, 0.0008, 0.9 ),
                    method = "MNL" )
eff8b

#Example
eff9a <- logEffInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                    allXVal = c( 1, NA, NA ),
                    xPos = c( 2, 3 ),
                    refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                    yCat = 2, method = "MNL" )
eff9a

eff9b <- logEffInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ),
                    allXVal = c( 1, NA, NA ),
                    xPos = c( 2, 3 ),
                    refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                    yCat = 2,
                    allCoefSE = c( 0.003, 0.045, 0.007, 0.009, 0.0008, 0.9 ),
                    method = "MNL" )
eff9b

#Example
eff10a <- logEffInt( allCoef = c( 0.2, 0.3, 0.5, 0.091 ),
                     allXVal = c( 1, NA, NA, 2.45, 1, NA, NA, 0.79 ),
                     xPos = c( 2, 3 ),
                     refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                     yCat = 2, method = "CondL" )
eff10a

eff10b <- logEffInt( allCoef = c( 0.2, 0.3, 0.5, 0.091 ),
                     allXVal = c( 1, NA, NA, 2.45, 1, NA, NA, 0.79 ),
                     xPos = c( 2, 3 ),
                     refBound = c( 8, 12 ), intBound = c( 13, 15 ),
                     allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ),
                     yCat = 2, method = "CondL" )
eff10b

# Example
matrix1 <- matrix( c( 1, NA, 0.3, 0.09, 1, NA, 0.9, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, NA, 0.099, 0.211 ), nrow = 4 )
eff12a <- logEffInt( allCoefBra = c( 0.445, 0.03, -0.002 ),
                     allCoef = c( 1.8, 0.005, -0.12, 0.8 ),
                     allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
                     allXVal = list( matrix1, matrix2 ),
                     refBound = c( 0.5, 1.5 ), intBound = c( 2, 3.5 ),
                     xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ),
                     method = "NestedL" )
eff12a

matrix1 <- matrix( c( 1, NA, 0.3, 0.09, 1, NA, 0.9, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, NA, 0.099, 0.211 ), nrow = 4 )
eff12b <- logEffInt( allCoefBra = c( 0.445, 0.03, -0.002 ),
                     allCoef = c( 1.8, 0.005, -0.12, 0.8 ),
                     allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
                     allXVal = list( matrix1, matrix2 ),
                     allCoefSE = c( 0.003, 0.045, 0.007, 0.0032 ),
                     refBound = c( 0.5, 1.5 ), intBound = c( 2, 3.5 ),
                     xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ),
                     method = "NestedL" )
eff12b

## -----------------------------------------------------------------------------
logEffGr <- function( allCoef, allXVal, xPos, Group, yCat = NA,
                      allCoefSE = rep( NA, length( allCoef ) ),
                      method = "binary" ){
  if( method == "binary" ){
     nCoef <- length( allCoef )
     xCoef <- allCoef[ xPos ]
     xShares <- allXVal[ xPos ]
  } else if( method == "MNL" ){
     nCoef <- length( allCoef )
     mCoef <- matrix( allCoef, nrow = length( allXVal ) )
     NCoef <- dim( mCoef )[2]
     pCoef <- dim( mCoef )[1]
     xCoef <- mCoef[ xPos, ]
     xShares <- allXVal[ xPos ]
  } else{
     nCoef <- length( allCoef )
     xCoef <- allCoef[ xPos ]
     mXVal <- matrix( allXVal, nrow = nCoef )
     pCoef <- dim( mXVal )[2]
     xShares <- mXVal[ xPos, ]
  }
  if( method == "binary" || method == "MNL" ){
    if( sum( xShares ) > 1 ){
      stop( "the shares in argument 'xShares' sum up to a value larger than 1" )
    }
  } else{
    for( p in 1:pCoef ){
      if( sum( xShares[ , p ] ) > 1 ){
        stop( "the shares in argument 'xShares' sum up to a value larger than 1" )
      }
    }
  }
  if( method == "binary" ){
    if( length( xCoef ) != length( xShares ) ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( length( xCoef ) != length( Group ) ){
      stop( "arguments 'xCoef' and 'Group' must have the same length" )
    }
  } else if( method == "MNL" ){
    if( dim( xCoef )[1] != length( xShares ) ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( dim( xCoef )[1] != length( Group ) ){
      stop( "arguments 'xCoef' and 'Group' must have the same length" )
    }
  } else{
    if( length( xCoef ) != dim( xShares )[1] ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( length( xCoef ) != length( Group ) ){
      stop( "arguments 'xCoef' and 'Group' must have the same length" )
    }
  }
  if( !all( Group %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'Group' must be -1, 0, or 1" )
  }
  if( method == "binary" ){
    # D_mr
    DRef <- sum( xCoef[ Group == -1 ] * xShares[ Group == -1 ]) /
            sum( xShares[ Group == -1 ] )
    XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DRef
    # D_ml
    DEffect <- sum( xCoef[ Group == 1 ] * xShares[ Group == 1 ]) /
               sum( xShares[ Group == 1 ] )
    XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DEffect
    # effect
    effeG <- exp( XBetaEffect )/( 1 + exp( XBetaEffect ) ) -
             exp( XBetaRef )/( 1 + exp( XBetaEffect ) )
  } else if( method == "MNL" ){
    # D_mr
    DRef <- colSums( xCoef[ Group == -1, , drop = FALSE ] *
                     xShares[ Group == -1 ] )/
                sum( xShares[ Group == -1 ] )
    XBetaRef <- colSums( mCoef[ -xPos, , drop = FALSE ] *
                         allXVal[ -xPos ]) + DRef
    # D_ml
    DEffect <- colSums( xCoef[ Group == 1, , drop = FALSE ] *
                        xShares[ Group == 1 ] )/
                   sum( xShares[ Group == 1 ] )
    XBetaEffect <- colSums( mCoef[ -xPos, , drop = FALSE ] *
                            allXVal[ -xPos ]) + DEffect
    # effect
    effeG <- exp( XBetaEffect[ yCat ] )/( 1 + sum( exp( XBetaEffect ) ) ) -
             exp( XBetaRef[ yCat ] )/( 1 + sum( exp( XBetaRef ) ) )
  } else{
    # D_mr
    DRef <- colSums( xCoef[ Group == -1 ] *
                     xShares[ Group == -1, , drop = FALSE ] )/
                sum( xShares[ Group == -1, , drop = FALSE ] )
    XBetaRef <- colSums( allCoef[ -xPos ] *
                         mXVal[ -xPos, , drop = FALSE ] ) + DRef
    # D_ml
    DEffect <- colSums( xCoef[ Group == 1 ] *
                        xShares[ Group == 1, , drop = FALSE ] )/
                   sum( xShares[ Group == 1, , drop = FALSE ] )
    XBetaEffect <- colSums( allCoef[ -xPos ] *
                            mXVal[ -xPos, , drop = FALSE ] ) + DEffect
    # effect
    effeG <- exp( XBetaEffect[ yCat ] )/( sum( exp( XBetaEffect ) ) ) -
             exp( XBetaRef[ yCat ] )/( sum( exp( XBetaRef ) ) )
  }
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  if( method == "binary" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( exp( XBetaEffect )/( 1 + exp( XBetaEffect ))^2 -
                           exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 ) *
                           allXVal[ -xPos ]
    derivCoef[ xPos ] = exp( XBetaEffect )/( 1 + exp( XBetaEffect))^2 * DEffect -
                        exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 * DRef
  } else if( method == "MNL" ){
    derivCoef <- matrix( NA, nrow=pCoef, ncol=NCoef )
    for( p in 1:NCoef ){
      if( p == yCat ){
        derivCoef[ -xPos, p ] <-
          ( exp( XBetaEffect[ p ] ) *
          ( 1 + sum( exp( XBetaEffect[ -yCat ] ) ) )/
          ( 1 + sum( exp( XBetaEffect ) ) )^2 -
            exp( XBetaRef[ p ] ) *
          ( 1 + sum( exp( XBetaRef[ -yCat ] ) ) )/
          ( 1 + sum( exp( XBetaRef ) ) )^2 ) * allXVal[ -xPos ]
        derivCoef[ xPos, p ] <-
          ( exp( XBetaEffect[ p ] ) *
          ( 1 + sum( exp( XBetaEffect[ -yCat ] ) ) )/
          ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * DEffect -
          ( exp( XBetaRef[ p ] ) *
          ( 1 + sum( exp( XBetaRef[ -yCat ] ) ) )/
          ( 1 + sum( exp( XBetaRef ) ) )^2 ) * DRef
      } else{
        derivCoef[ -xPos, p ] <-
          ( ( exp( XBetaRef[ yCat ] ) * exp( XBetaRef[ p ] ) )/
            ( 1 + sum( exp( XBetaRef ) ) )^2 -
            ( exp( XBetaEffect[ yCat ] ) * exp( XBetaEffect[ p ] ) )/
            ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * allXVal[ -xPos ]
        derivCoef[ xPos, p ] <-
            ( ( exp( XBetaRef[ yCat ] ) * exp( XBetaRef[ p ] ) )/
              ( 1 + sum( exp( XBetaRef ) ) )^2 ) * DRef -
            ( ( exp( XBetaEffect[ yCat ] ) * exp( XBetaEffect[ p ] ) )/
              ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * DEffect
      }
    }
    derivCoef <- c( derivCoef )
  } else{
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( exp( XBetaEffect[ yCat ] ) * mXVal[ -xPos, yCat ] *
                           sum( exp( XBetaEffect ) ) -
                           exp( XBetaEffect[ yCat ] ) * sum( exp( XBetaEffect ) *
                           mXVal[ -xPos, ] ) )/
                         ( sum( exp( XBetaEffect ) ) )^2 -
                         ( exp( XBetaRef[ yCat ] ) * mXVal[ -xPos, yCat ] *
                           sum( exp( XBetaRef ) ) -
                           exp( XBetaRef[ yCat ] ) * sum( exp( XBetaRef ) *
                           mXVal[ -xPos, ] ) )/
                         ( sum( exp( XBetaRef ) ) )^2
    derivCoef[ xPos ] =  ( exp( XBetaEffect[ yCat ] ) * DEffect[ yCat ] *
                           sum( exp( XBetaEffect ) ) -
                           exp( XBetaEffect[ yCat ] ) *
                           sum( exp( XBetaEffect ) * DEffect[ yCat ] ) )/
                         ( sum( exp( XBetaEffect ) ) )^2 -
                         ( exp( XBetaRef[ yCat ] ) * DRef[ yCat ] *
                           sum( exp( XBetaRef ) ) -
                           exp( XBetaRef[ yCat ] ) *
                           sum( exp( XBetaRef ) * DRef[ yCat ] ) )/
                         ( sum( exp( XBetaRef ) ) )^2
  }
  # variance covariance of the coefficients (covariances set to zero)
  vcovCoef <- diag( allCoefSE^2 )
  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}

## -----------------------------------------------------------------------------
# Example
eff10a <- logEffGr( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034,
                                -0.005, 0.89, -1.2 ),
                    allXVal = c( 1, 0.1, 0.3, 0.25, 0.15, 0.2, 2.34, 10.8 ),
                    xPos = c( 2:6 ), Group = c( 0, -1, -1, 1, 1 ) )
eff10a

eff10b <- logEffGr( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034,
                                -0.005, 0.89, -1.2 ),
                    allXVal = c( 1, 0.1, 0.3, 0.25, 0.15, 0.2, 2.34, 10.8 ),
                    xPos = c( 2:6 ), Group = c( 0, -1, -1, 1, 1 ),
                    allCoefSE = c( 0.03, 0.0001, 0.005, 0, 0.01,
                                   0.004, 0.05, 0.8 ) )
eff10b

# Example
eff11a <- logEffGr( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034,
                                -0.005, 0.89, 0, 0.005, 0.06, 1.7, 0 ),
                    allXVal = c( 1, 0.5, 0.3, 0.2 ), xPos = c( 2:4 ),
                    Group = c( -1, -1, 1 ), yCat = 2, method = "MNL" )
eff11a

eff11b <- logEffGr( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034,
                                -0.005, 0.89, 0, 0.005, 0.06, 1.7, 0 ),
                    allXVal = c( 1, 0.5, 0.3, 0.2 ), xPos = c( 2:4 ),
                    Group = c( -1, -1, 1 ), yCat = 2, method = "MNL",
                    allCoefSE = c( 0.03, 0.0001, 0.005, 0, 0.01, 0.004,
                                   0.05, 0, 0.004, 0.5, 0.0078, 0 ) )
eff11b

# Example
eff12a <- logEffGr( allCoef = c( 0.28, 0.003, 0.0075, 0 ),
                    allXVal = c( 1, 0.5, 0.3, 0.2, 1, 0.4, 0.4, 0.1 ),
                    xPos = c( 2:4 ),
                    Group = c( -1, -1, 1 ), yCat = 2, method = "CondL" )
eff12a

eff12b <- logEffGr( allCoef = c( 0.28, 0.003, 0.0075, 0 ),
                    allXVal = c( 1, 0.5, 0.3, 0.2, 1, 0.4, 0.4, 0.1 ),
                    xPos = c( 2:4 ),
                    allCoefSE = c( 0.03, 0.0001, 0.005, 0 ),
                    Group = c( -1, -1, 1 ), yCat = 2, method = "CondL" )
eff12b

## ----checkXPos,eval=FALSE-----------------------------------------------------
# checkXPos <- function( xPos, minLength, maxLength, minVal, maxVal,
#   requiredVal = NA ) {
#   if( any( xPos != round( xPos ) ) ) {
#     stop( "argument 'xPos' must be a vector of integers" )
#   }
#   if( length( xPos ) < minLength ) {
#     stop( "argument 'xPos' must have a length equal to or larger than ",
#       minLength )
#   }
#   if( length( xPos ) > maxLength ) {
#     stop( "argument 'xPos' must have a length smaller than or equal to ",
#       maxLength )
#   }
#   if( any( xPos < minVal ) ) {
#     stop( "all elements of argument 'xPos' must be equal to or larger than ",
#       minVal )
#   }
#   if( any( xPos > maxVal ) ) {
#     stop( "all elements of argument 'xPos' must be smaller than or equal to ",
#       maxVal )
#   }
#   if( max( table( xPos ) ) > 1 ) {
#     stop( "all elements of argument 'xPos' may only occur once" )
#   }
#   if( !is.na( requiredVal ) ) {
#     if( sum( xPos == requiredVal ) != 1 ) {
#       stop( "argument 'xPos' must have exactly one element that is ",
#         requiredVal )
#     }
#   }
# }

## ----checkXBeta,eval=FALSE----------------------------------------------------
# checkXBeta <- function( xBeta ) {
#   if( any( abs( xBeta ) > 3.5 ) ) {
#     warning( "At least one x'beta has an implausible value: ",
#       paste( xBeta, collapse = ", " ) )
#   }
# }

