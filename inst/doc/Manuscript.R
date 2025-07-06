## ----dataLoad,results='hide', echo=FALSE--------------------------------------
data( "Mroz87", package = "sampleSelection" )

## ----dataLfp3, echo=FALSE, results='hide'-------------------------------------
Mroz87$lfp3 <- factor( ifelse( Mroz87$hours == 0, "no",
  ifelse( Mroz87$hours <= 1300, "part", "full" ) ),
  ordered = TRUE, levels = c( "no", "part", "full" ) )

## ----dataKids, echo=FALSE, results='hide'-------------------------------------
Mroz87$kids <- Mroz87$kids5 + Mroz87$kids618

## ----dataAgeInt, echo=FALSE, results='hide'-----------------------------------
Mroz87$age30.37 <- Mroz87$age >= 30 & Mroz87$age <= 37
Mroz87$age38.44 <- Mroz87$age >= 38 & Mroz87$age <= 44
Mroz87$age45.52 <- Mroz87$age >= 45 & Mroz87$age <= 52
Mroz87$age53.60 <- Mroz87$age >= 53 & Mroz87$age <= 60
all.equal(
  Mroz87$age30.37 + Mroz87$age38.44 + Mroz87$age45.52 + Mroz87$age53.60,
  rep( 1, nrow( Mroz87 ) ) )

## ----tabDescStat,results='asis', echo=FALSE, message = FALSE------------------
library( "stargazer" )
stargazer( Mroz87, title = "Descriptive Statistics",
  label = "tab:descStats",
  float.env = "sidewaystable", header = FALSE, digits = 2, align = TRUE,
  omit.summary.stat = c( "p25", "p75" ) )

## ----estProbit, echo=FALSE, results='hide'------------------------------------
estProbit <- glm( lfp ~ kids + age + educ,
  family = binomial(link = "probit"), data = Mroz87 )
xMean <- c( 1, colMeans( Mroz87[ , c( "kids", "age", "educ" ) ] ) )

estProbitQ <- glm( lfp ~ kids + age + I(age^2) + educ,
  family = binomial(link = "probit"), data = Mroz87 )
xMeanQ <- c( xMean[ 1:3], xMean[3]^2, xMean[4] )

## ----tabEstProbit, echo=FALSE, results='asis'---------------------------------
stargazer( estProbit, estProbitQ,
  title = "Probit regression results with age as linear and quadratic covariate",
  label = "tab:probitSum", header = FALSE,
  digits = 2, intercept.bottom = FALSE )

## ----echo=FALSE, results='hide'-----------------------------------------------
library( "urbin" )

## -----------------------------------------------------------------------------
urbinEla( coef(estProbit), xMean, xPos = 3, model = "probit" )

## -----------------------------------------------------------------------------
urbinEla( coef(estProbitQ), xMeanQ, xPos = c( 3, 4 ),
  model = "probit" )

## -----------------------------------------------------------------------------
urbinEla( coef(estProbit), xMean, xPos = 3, model = "probit",
  allCoefVcov = vcov(estProbit)  )
urbinEla( coef(estProbitQ), xMeanQ, xPos = c( 3, 4 ),
  model = "probit", allCoefVcov = vcov(estProbitQ) )

## -----------------------------------------------------------------------------
urbinEla( coef(estProbit), xMean, xPos = 3, model = "probit",
  allCoefVcov = sqrt(diag(vcov(estProbit))),
  seSimplify = FALSE  )
urbinEla( coef(estProbitQ), xMeanQ, xPos = c( 3, 4 ),
  model = "probit", allCoefVcov = sqrt(diag(vcov(estProbitQ))),
  seSimplify = FALSE )

## -----------------------------------------------------------------------------
urbinEla( coef(estProbit), xMean, xPos = 3, model = "probit",
  allCoefVcov = sqrt(diag(vcov(estProbit))) )
urbinEla( coef(estProbitQ), xMeanQ, xPos = c( 3, 4 ),
  model = "probit", allCoefVcov = sqrt(diag(vcov(estProbitQ))),
  xMeanSd = c( mean(Mroz87$age), sd(Mroz87$age) ) )

## ----estLogitInt, echo=FALSE, results='hide'----------------------------------
estLogitInt <- glm( lfp ~ kids + age30.37 + age38.44 + age53.60 + educ,
  family = binomial(link = "logit"), data = Mroz87 )
xMeanInt <- c( xMean[1:2], mean( Mroz87$age30.37 ),
  mean( Mroz87$age38.44 ), mean( Mroz87$age53.60 ), xMean[4] )

## ----tabEstLogitInt, echo=FALSE, results='asis'-------------------------------
stargazer( estLogitInt,
  title = "Logistic regression results with age as interval-coded covariate",
  label = "tab:logitInt", header = FALSE,
  digits = 2, intercept.bottom = FALSE )

## -----------------------------------------------------------------------------
urbinElaInt( coef(estLogitInt), xMeanInt, xPos = c( 3, 4, 0, 5 ),
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), model = "logit" )

## -----------------------------------------------------------------------------
urbinElaInt( coef(estLogitInt), xMeanInt, xPos = c( 3, 4, 0, 5 ),
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), model = "logit",
  allCoefVcov = vcov(estLogitInt) )

## -----------------------------------------------------------------------------
urbinElaInt( coef(estLogitInt), xMeanInt, xPos = c( 3, 4, 0, 5 ),
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), model = "logit",
  allCoefVcov = sqrt(diag(vcov(estLogitInt))) )

## -----------------------------------------------------------------------------
urbinEffInt( coef(estProbit), xMean, xPos = 3,
  refBound = c( 30, 44 ), intBound = c( 53, 60 ),
  model = "probit" )

## -----------------------------------------------------------------------------
urbinEffInt( coef(estProbitQ), xMeanQ, xPos = c( 3, 4 ),
  refBound = c( 30, 44 ), intBound = c( 53, 60 ),
  model = "probit" )

## -----------------------------------------------------------------------------
urbinEffInt( coef(estProbit), xMean, xPos = 3,
  refBound = c( 30, 44 ), intBound = c( 53, 60 ),
  model = "probit", allCoefVcov = vcov(estProbit) )
urbinEffInt( coef(estProbitQ), xMeanQ, xPos = c( 3, 4 ),
  refBound = c( 30, 44 ), intBound = c( 53, 60 ),
  model = "probit", allCoefVcov = vcov(estProbitQ) )

## -----------------------------------------------------------------------------
urbinEffInt( coef(estProbit), xMean, xPos = 3,
  refBound = c( 30, 44 ), intBound = c( 53, 60 ),
  model = "probit", allCoefVcov = sqrt(diag(vcov(estProbit))) )
urbinEffInt( coef(estProbitQ), xMeanQ, xPos = c( 3, 4 ),
  refBound = c( 30, 44 ), intBound = c( 53, 60 ),
  model = "probit", allCoefVcov = sqrt(diag(vcov(estProbitQ))) )

## -----------------------------------------------------------------------------
urbinEffInt( coef(estProbitQ), xMeanQ, xPos = c( 3, 4 ),
  refBound = c( 30, 44 ), intBound = c( 53, 60 ),
  model = "probit", allCoefVcov = sqrt(diag(vcov(estProbitQ))),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )

## -----------------------------------------------------------------------------
urbinEffCat( coef(estLogitInt), xMeanInt, xPos = c( 3:5 ),
  xGroups = c( -1, -1, 1, 0 ), model = "logit" )

## -----------------------------------------------------------------------------
urbinEffCat( coef(estLogitInt), xMeanInt, c( 3:5 ),
  c( -1, -1, 1, 0 ), vcov(estLogitInt), model = "logit" )

## -----------------------------------------------------------------------------
urbinEffCat( coef(estLogitInt), xMeanInt, c( 3:5 ),
  c( -1, -1, 1, 0 ), sqrt(diag(vcov(estLogitInt))),
  model = "logit" )

## ----estOProbit,echo=FALSE, results='hide'------------------------------------
library( "MASS" )
estOProbitQ <- polr( lfp3 ~ kids + age + I(age^2) + educ,
  data = Mroz87, method = "probit", Hess = TRUE )
xMeanOProbit <- c( xMeanQ, -1 )

## ----tabEstOProbit,echo=FALSE, results='asis'---------------------------------
stargazer( estOProbitQ,
  title = "Ordered probit regression results with age as linear and quadratic covariate",
  label = "tab:estOProbit", header = FALSE,
  ord.intercepts = TRUE, digits = 2 )

## ----estMLogitInt,echo=FALSE, results='hide', message = FALSE-----------------
library( "mlogit" )
estMLogitInt <- mlogit(
  lfp3 ~ 0 | kids + age30.37 + age38.44 + age53.60 + educ,
  data = Mroz87, reflevel = "no", shape = "wide" )

## ----tabEstMLogit,echo=FALSE, results='asis'----------------------------------
stargazer( estMLogitInt,
  title = "Multinomial logistic regression results with age as interval-coded covariate",
  label = "tab:estMLogitInt", header = FALSE,
  digits = 2 )

## -----------------------------------------------------------------------------
urbinEla( coef(summary(estOProbitQ))[-6,1], c( xMeanQ[-1], -1 ),
  xPos = c( 2, 3 ), iPos = 5, model = "oprobit",
  vcov(estOProbitQ)[-6,-6] )

## -----------------------------------------------------------------------------
coefPermuteInt <- c( seq( 1, 11, 2 ), seq( 2, 12, 2 ) )
urbinElaInt( coef(estMLogitInt)[coefPermuteInt], xMeanInt,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "mlogit",
  vcov(estMLogitInt)[coefPermuteInt,coefPermuteInt],
  yCat = c( 1, 2 ) )

## ----eval=FALSE---------------------------------------------------------------
# data( "Mroz87", package = "sampleSelection" )

## ----eval=FALSE---------------------------------------------------------------
# Mroz87$kids <- Mroz87$kids5 + Mroz87$kids618

## ----eval=FALSE---------------------------------------------------------------
# Mroz87$age30.37 <- Mroz87$age >= 30 & Mroz87$age <= 37
# Mroz87$age38.44 <- Mroz87$age >= 38 & Mroz87$age <= 44
# Mroz87$age45.52 <- Mroz87$age >= 45 & Mroz87$age <= 52
# Mroz87$age53.60 <- Mroz87$age >= 53 & Mroz87$age <= 60
# all.equal(
#   Mroz87$age30.37 + Mroz87$age38.44 + Mroz87$age45.52 + Mroz87$age53.60,
#   rep( 1, nrow( Mroz87 ) ) )

## ----eval=FALSE---------------------------------------------------------------
# Mroz87$lfp3 <- factor( ifelse( Mroz87$hours == 0, "no",
#   ifelse( Mroz87$hours <= 1300, "part", "full" ) ),
#   ordered = TRUE, levels = c( "no", "part", "full" ) )

## ----eval=FALSE---------------------------------------------------------------
# estProbit <- glm( lfp ~ kids + age + educ,
#   family = binomial(link = "probit"), data = Mroz87 )
# xMean <- c( 1, colMeans( Mroz87[ , c( "kids", "age", "educ" ) ] ) )
# 
# estProbitQ <- glm( lfp ~ kids + age + I(age^2) + educ,
#   family = binomial(link = "probit"), data = Mroz87 )
# xMeanQ <- c( xMean[ 1:3], xMean[3]^2, xMean[4] )

## ----eval=FALSE---------------------------------------------------------------
# estLogitInt <- glm( lfp ~ kids + age30.37 + age38.44 + age53.60 + educ,
#   family = binomial(link = "logit"), data = Mroz87 )
# xMeanInt <- c( xMean[1:2], mean( Mroz87$age30.37 ),
#   mean( Mroz87$age38.44 ), mean( Mroz87$age53.60 ), xMean[4] )

## ----eval=FALSE---------------------------------------------------------------
# library( "MASS" )
# estOProbitQ <- polr( lfp3 ~ kids + age + I(age^2) + educ,
#   data = Mroz87, method = "probit", Hess = TRUE )
# xMeanOProbit <- c( xMeanQ, -1 )

## ----eval=FALSE---------------------------------------------------------------
# library( "mlogit" )
# estMLogitInt <- mlogit(
#   lfp3 ~ 0 | kids + age30.37 + age38.44 + age53.60 + educ,
#   data = Mroz87, reflevel = "no", shape = "wide" )

