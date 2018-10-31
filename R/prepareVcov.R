prepareVcov <- function( allCoefVcov, nCoef, xPos, xMeanSd = NULL, 
  nXVal = nCoef, iPos = 0, pCall = NULL ){

  if( !is.null( pCall ) ) {
    pCall <- paste( "In", deparse( pCall, width.cutoff = 500 ), ":\n  " )
  }

  # coefficients: number of 'replications' of explanatory variables 
  nCoefRep <- max( 1, round( nCoef / nXVal ) )
  if( nCoef != nCoefRep * nXVal && nCoef > nXVal ) {
    stop( "internal error: nCoef is not a multiple of nXVal" )
  }
  
  errMsgVcov <- paste( "argument 'allCoefVcov' must be a vector of length",
    nCoef, "or a symmetric matrix with dimension", nCoef )
  if( is.null( allCoefVcov ) ) {
    allCoefVcov <- matrix( NA, nrow = nCoef, ncol = nCoef )
    if( !is.null( xMeanSd ) ) {
      warning( pCall, "argument 'xMeanSd' is ignored,",
        " because argument 'allCoefVcov' has not been specified",
        call. = is.null( pCall ) )
    }
  } else if( is.matrix( allCoefVcov ) ) {
    if( nrow( allCoefVcov ) != nCoef || ncol( allCoefVcov ) != nCoef ) {
      stop( errMsgVcov )
    }
    if( !is.null( xMeanSd ) ) {
      warning( pCall, "argument 'xMeanSd' is ignored,",
        " the full variance-covariance matrix has been specified",
        " by argument 'allCoefVcov'",
        call. = is.null( pCall ) )
    }
  } else if( is.vector( allCoefVcov ) ) {
    if( length( allCoefVcov ) != nCoef ) {
      stop( errMsgVcov )
    } else {
      if( nCoef == 1 ) {
        allCoefVcov <- allCoefVcov^2
      } else {
        allCoefVcov <- diag( allCoefVcov^2 )
      }
    }
    # if( nCoefRep > 1 ) {
    #   for( xNo in 1:nXVal ) {
    #     for( i in 1:( nCoefRep - 1 ) ) {
    #       for( j in (i+1):nCoefRep ) {
    #         rowNo <- (i-1) * nXVal + xNo
    #         colNo <- (j-1) * nXVal + xNo
    #         allCoefVcov[ rowNo, colNo ] <- allCoefVcov[ colNo, rowNo ] <-
    #           0.35 * sqrt( allCoefVcov[ rowNo, rowNo ] ) * 
    #           sqrt( allCoefVcov[ colNo, colNo ] ) 
    #       }
    #     }
    #   }
    # }
    if( is.null( xMeanSd ) ) {
      if( length( xPos ) == 2 ) {
        warning( pCall, "the returned standard error is likely largely upward biased",
          " and, thus, in most cases meaningless;",
          " you can provide the full covariance matrix", 
          " via argument 'allCoefVcov' to avoid this bias",
          " or use argument 'xMeanSd' to substantially reduce this bias",
          call. = is.null( pCall ) )
      }
    } else {
      if( length( xPos ) != 2 ) {
        warning( pCall, "argument 'xMeanSd' is ignored,",
          " because the model does not include a quadratic term",
          " of the explanatory variable of interest",
          call. = is.null( pCall ) )
      } else {
        if( length( xMeanSd ) != 2 || !is.numeric( xMeanSd ) ) {
          stop( "argument 'xMeanSd' must be a vector of two numeric values")
        }
        set.seed( 123 )
        x <- rnorm( 1000, xMeanSd[ 1 ], xMeanSd[ 2 ] )
        X <- cbind( 1, x, x^2 )
        XXinv <- solve( t(X) %*% X )
        for( i in 1:nCoefRep ) {
          xPosRep <- (i-1) * nXVal + xPos
          sigmaSq <- sqrt( ( allCoefVcov[ xPosRep[1], xPosRep[1] ] / XXinv[2,2] ) * 
              ( allCoefVcov[ xPosRep[2], xPosRep[2] ] / XXinv[3,3] ) )
          allCoefVcov[ xPosRep[1], xPosRep[2] ] <- 
            allCoefVcov[ xPosRep[2], xPosRep[1] ] <-
            sigmaSq * XXinv[2,3]
          if( iPos != 0 ) {
            iPosRep <- (i-1) * nXVal + iPos
            allCoefVcov[ iPosRep, xPosRep[1] ] <- 
              allCoefVcov[ xPosRep[1], iPosRep ] <-
              sigmaSq * XXinv[1,2]
            allCoefVcov[ iPosRep, xPosRep[2] ] <- 
              allCoefVcov[ xPosRep[2], iPosRep ] <-
              sigmaSq * XXinv[1,3]
          }
        }
      }
    }
  } else {
    stop( errMsgVcov )
  }
  return( allCoefVcov )
} 
