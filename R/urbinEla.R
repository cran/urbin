urbinEla <- function( allCoef, allXVal, xPos, model, 
  allCoefVcov = NULL, seSimplify = !is.matrix( allCoefVcov ), 
  xMeanSd = NULL, iPos = 1, yCat = NULL ){
  
  # check argument seSimplify
  if( length( seSimplify ) != 1 || !is.logical( seSimplify ) ) {
    stop( "argument 'seSimplify' must be TRUE or FALSE" )
  }
  if( model == "lpm" && !is.null( allCoefVcov ) && 
      ( seSimplify == is.matrix( allCoefVcov ) ) ) {
    warning( "argument 'seSimplify' is ignored in linear probability models" )    
  } else if( seSimplify && is.matrix( allCoefVcov ) ) {
    warning( "if the (full) covariance matrix is provided",
      " via argument 'allCoefCov'", 
      " it is NOT recommended to set argument 'seSimplify' to TRUE" )    
  } else if( !seSimplify && !is.null( allCoefVcov) && 
      !is.matrix( allCoefVcov ) ) {
    warning( "the returned standard error is likely very imprecise;",
      " you can provide the full covariance matrix", 
      " via argument 'allCoefVcov'",
      " or do NOT set argument 'seSimplify' to FALSE",
      " to obtain a more precise standard error" )
  }
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # check allXVal and allCoef
  if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
    # LPM model: allXVal can be a scalar even if there is a quadratic term
    if( model == "lpm" && length( xPos ) == 2 && length( allXVal ) == 1 ){
      allXVal <- c( allXVal, allXVal^2 )
      nXVal <- length( allXVal )
    }
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
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, 
    maxVal = ifelse( model == "mlogit", nXVal, nCoef ) )
  # check position of the intercept
  checkIPos( iPos, xPos, allXVal, model ) 
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, length( allCoef ), xPos, xMeanSd,
    nXVal = nXVal, iPos = iPos, pCall = match.call() )
  # Identify coefficients of interest (kth/tth covariate)
  if( length( xPos ) == 2 ){
    if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
      xCoef <- allCoef[ xPos ]
      if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
        stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
          " to the squared value of 'allXVal[ xPos[1] ]'" )
      }
    } else if( model == "mlogit" ){
      xCoef <- mCoef[ xPos, ]
      if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
        stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
          " to the squared value of 'allXVal[ xPos[1] ]'" ) 
      }    
    } else {
      stop( "argument 'model' specifies an unknown type of model" )
    }
  } else if( length( xPos ) == 1 ) {
    if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
      xCoef <- c( allCoef[ xPos ], 0 )
    } else if( model == "mlogit" ){
      xCoef <- rbind( mCoef[ xPos, ], 0 ) 
    } else {
      stop( "argument 'model' specifies an unknown type of model" )
    }    
  } else {
    stop( "argument 'xPos' must be a scalar or a vector with two elements" )
  }
  # prepare calculation of semi-elasticity 
  if( model == "lpm" ) {
    xVal <- allXVal[ xPos[ 1 ] ]
    semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal ) * xVal
  } else if( model %in% c( "probit", "oprobit" ) ){
    xVal <- allXVal[ xPos[ 1 ] ]
    xBeta <- sum( allCoef * allXVal )
    checkXBeta( xBeta )
    dfun <- dnorm( xBeta )
    semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  } else if( model == "logit" ){
    xVal <- allXVal[ xPos[1] ]
    xBeta <- sum( allCoef * allXVal )
    checkXBeta( xBeta )
    dfun <- exp( xBeta )/( 1 + exp( xBeta ) )^2
    semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  } else if( model == "mlogit" ){
    xVal <- allXVal[ xPos[1] ]
    xCoefLinQuad <- xCoef[ 1, ] + 2 * xCoef[ 2, ] * xVal
    xBeta <- allXVal %*% mCoef
    checkXBeta( xBeta )
    pfun <- exp( xBeta ) / sum( exp( xBeta ) )
    semElaCat <- xVal * pfun * 
      ( xCoefLinQuad - sum( xCoefLinQuad * pfun ) )
    semEla <- sum( semElaCat[ yCat ] )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  } 
  # partial derivatives of semi-elasticities wrt coefficients
  if( model == "lpm" ){
    derivCoef <- rep( 0, nCoef ) 
    derivCoef[ xPos[1] ] <- xVal
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- 2 * xVal^2
    }
  } else if( model %in% c( "probit", "oprobit" ) ){
    if( seSimplify ) {
      derivCoef <- rep( 0, length( allCoef ) )
    } else {
      derivCoef <- - semEla * xBeta * allXVal
    }
    derivCoef[ xPos[1] ] <- derivCoef[ xPos[1] ] + dnorm( xBeta ) * xVal
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- derivCoef[ xPos[2] ] + dnorm( xBeta ) * 2 * xVal^2
    }
  } else if( model == "logit" ){
    if( seSimplify ) {
      derivCoef <- rep( 0, length( allCoef ) )
    } else {
      derivCoef <- allXVal * semEla *
        ( 1 - 2 * exp( xBeta ) / ( 1 + exp( xBeta ) ) )
    }
    derivCoef[ xPos[1] ] <- derivCoef[ xPos[1] ] + dfun * xVal
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- derivCoef[ xPos[2] ] + dfun * 2 * xVal^2
    }
  } else if( model == "mlogit" ){
    derivCoef <- rep( 0, length( allCoef ) )
    for( p in 1:nYCat ) {
      coefNoYCat <- ( 1 + (p-1)*nXVal ):( p * nXVal ) 
      for( yCati in yCat ) {
        derivCoef[ coefNoYCat ][ xPos[1] ] <- 
          derivCoef[ coefNoYCat ][ xPos[1] ] -
          pfun[ yCati ] * xVal * pfun[ p ] + 
          ( p == yCati ) * pfun[ yCati ] * xVal
        if( length( xPos ) == 2 ) {
          derivCoef[ coefNoYCat ][ xPos[2] ] <- 
            derivCoef[ coefNoYCat ][ xPos[2] ] -
            pfun[ yCati ] * 2 * xVal^2 * pfun[ p ] +
            ( p == yCati ) * pfun[ yCati ] * 2 * xVal^2
        }
        if( !seSimplify ) {
          derivCoef[ coefNoYCat ][ -xPos ] <-
            derivCoef[ coefNoYCat ][ -xPos ] +
            ( - pfun[ yCati ] * semElaCat[ p ] -
                pfun[ p ] * semElaCat[ yCati ] +
                ( p == yCati ) * semElaCat[ p ] ) *
            allXVal[ -xPos ]
          derivCoef[ coefNoYCat ][ xPos[1] ] <-
            derivCoef[ coefNoYCat ][ xPos[1] ] +
            ( - pfun[ yCati ] * semElaCat[ p ] -
                pfun[ p ] * semElaCat[ yCati ] +
                ( p == yCati ) * semElaCat[ p ] ) * 
            xVal
          if( length( xPos ) == 2 ) {
            derivCoef[ coefNoYCat ][ xPos[2] ] <-
              derivCoef[ coefNoYCat ][ xPos[2] ] +
              ( - pfun[ yCati ] * semElaCat[ p ] -
                  pfun[ p ] * semElaCat[ yCati ] + 
                  ( p == yCati ) * semElaCat[ p ] ) * 
              xVal^2
          }
        }
      }
    }
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
