urbinEffCat <- function( allCoef, allXVal, xPos, xGroups, model,
  allCoefVcov = NULL, iPos = 1, yCat = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = nCoef, minVal = 1, 
    maxVal = ifelse( model == "mlogit", nXVal, nCoef ) )
  # check position of the intercept
  checkIPos( iPos, xPos, allXVal, model ) 
  # number of categories
  nCat <- length( xPos ) + 1
  # check allXVal and allCoef
  if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
    xCoef <- c( allCoef[ xPos ], 0 )
    if( nXVal != nCoef ){
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
    # coefficients of the dummy variables (variable of interest)
    xCoef <- rbind( mCoef[ xPos, ], 0 )
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
  # shares in each category
  xShares <- allXVal[ xPos ]
  if( any( xShares < 0 ) ){
    stop( "the share of the observations in at least one category",
      " is negative" )
  }
  if( sum( xShares ) > 1 ){
    stop( "the shares of the observations in the individual categories",
      " (without the reference category) sum up to a value larger than 1" )
  }
  xShares <- c( xShares, 1 - sum( xShares ) )
  if( length( xShares ) != nCat ) {
    stop( "internal error: length of shares not equal to number of categories" )
  }
  if( xShares[ nCat ] < 0.05  ) {
    warning( "there are only ", 100 * xShares[ length( xShares ) ],
      "% of the observations in the reference category --",
      " please check whether this is indeed the case" )
  }
  if( xShares[ nCat ] > 0.95  ) {
    warning( "there are ", 100 * xShares[ length( xShares ) ],
      "% of the observations in the reference category --",
      " please check whether this is indeed the case" )
  }
  # check argument xGroups
  if( length( xGroups ) != ( length( xPos ) + 1 ) ){
    stop( "the vector specified by argument 'xGroups' must have",
      " one more element that the vector specified by argument 'xPos'" )
  }
  if( ! all( xGroups %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'xGroups' must be -1, 0, or 1" )
  }
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos = NA, 
    pCall = match.call() )
  if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
    # D_mr
    DRef <- ifelse( xGroups == -1, xShares, 0 ) / 
      sum( xShares[ xGroups == -1 ] )
    # D_ml
    DEffect <- ifelse( xGroups == 1, xShares, 0 ) / 
      sum( xShares[ xGroups == 1 ] )
    if( model == "lpm" ) {
      # effect: sum of delta_m * ( D_ml - D_mr )
      effeG <- sum( xCoef * ( DEffect - DRef ) )
    } else {
      # linear predictors
      XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + 
        sum( DRef * xCoef )
      XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + 
        sum( DEffect * xCoef )
      checkXBeta( c( XBetaRef, XBetaEffect ) )
      if( model %in% c( "probit", "oprobit" ) ) {
        # effect
        effeG <- pnorm( XBetaEffect ) - pnorm( XBetaRef )
      } else if( model == "logit" ) {
        # effect
        effeG <- exp( XBetaEffect )/( 1 + exp( XBetaEffect ) ) - 
          exp( XBetaRef )/( 1 + exp( XBetaRef ) )
      }
    }
  } else if( model == "mlogit" ){
    # D_mr  
    DRef <- ifelse( xGroups == -1, xShares, 0 ) / 
      sum( xShares[ xGroups == -1 ] )
    XBetaRef <- allXVal[ -xPos ] %*% mCoef[ -xPos, , drop = FALSE ] + 
      DRef %*% xCoef
    # D_ml
    DEffect <- ifelse( xGroups == 1, xShares, 0 ) / 
      sum( xShares[ xGroups == 1 ] )
    XBetaEffect <- allXVal[ -xPos ] %*% mCoef[ -xPos, , drop = FALSE ] + 
      DEffect %*% xCoef
    checkXBeta( c( XBetaRef, XBetaEffect ) )
    # probabilities
    pFunRef <- exp( XBetaRef ) / sum( exp( XBetaRef ) )
    pFunEffect <- exp( XBetaEffect ) / sum( exp( XBetaEffect ) )
    # effect
    effeG <- sum( pFunEffect[ yCat ] - pFunRef[ yCat ] )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # partial derivative of the effect w.r.t. all estimated coefficients
  if( model == "lpm" ) {
    derivCoef <- rep( 0, nCoef )
    derivCoef[ xPos ] <- DEffect[ -nCat ] - DRef[ -nCat ]
  } else if( model %in% c( "probit", "oprobit" ) ) {
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( dnorm( XBetaEffect ) - dnorm( XBetaRef ) ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] = dnorm( XBetaEffect ) * DEffect[ -nCat ] - 
      dnorm( XBetaRef ) * DRef[ -nCat ]
  } else if( model == "logit" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( XBetaEffect )/( 1 + exp( XBetaEffect ))^2 - 
        exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] <- exp( XBetaEffect )/( 1 + exp( XBetaEffect))^2 * DEffect[ -nCat ] - 
      exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 * DRef[ -nCat ]
  } else if( model == "mlogit" ){
    derivCoef <- matrix( 0, nrow = nXVal, ncol = nYCat )
    for( p in 1:nYCat ){
      for( yCati in yCat ) {
        if( p == yCati ){
          derivCoef[ -xPos, p ] <- derivCoef[ -xPos, p ] +
            ( pFunEffect[ p ] - pFunEffect[ p ]^2 - 
                pFunRef[ p ] + pFunRef[ p ]^2 ) *
            allXVal[ -xPos ]
          derivCoef[ xPos, p ] <- derivCoef[ xPos, p ] +
            ( pFunEffect[ p ] - pFunEffect[ p ]^2 ) * DEffect[-nCat] -
            ( pFunRef[ p ] - pFunRef[ p ]^2 ) * DRef[-nCat]
        } else {  
          derivCoef[ -xPos, p ] <- derivCoef[ -xPos, p ] +
            ( pFunRef[ yCati ] * pFunRef[ p ] -
                pFunEffect[ yCati ] * pFunEffect[ p ] ) *
            allXVal[ -xPos ]
          derivCoef[ xPos, p ] <- derivCoef[ xPos, p ] +
            pFunRef[ yCati ] * pFunRef[ p ] * DRef[-nCat] -
            pFunEffect[ yCati ] * pFunEffect[ p ] * DEffect[-nCat]
        }
      }     
    }
    derivCoef <- c( derivCoef )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }

  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  
  # object to be returned
  result <- list()
  result$call <- match.call()
  result$allCoefVcov <- allCoefVcov
  result$derivCoef <- derivCoef
  result$effect <- effeG
  result$stdEr <- effeGSE
  class( result ) <- "urbin"
  return( result )
}
