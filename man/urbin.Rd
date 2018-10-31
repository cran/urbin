\name{urbin}
\alias{urbinEla}
\alias{urbinElaInt}
\alias{urbinEffInt}
\alias{urbinEffCat}

\title{Unifying Estimation Results with Binary Dependent Variables}

\description{
  These four functions calculate 
  the semi-elasticities and effects of explanatory variables
  in linear probability models, binary probit models,
  ordered probit models, multivariate probit models, binary logit models,
  and multinomial logit models.

  \code{urbinEla()} calculates 
  the semi-elasticity of a continuous variable
  that is used as a linear explanatory variable
  or as a linear and quadratic explanatory variable.

  \code{urbinElaInt()} calculates 
  the semi-elasticity of an interval-coded explanatory variable.

  \code{urbinEffInt()} calculates 
  the effect of a continuous variable
  that is used as a linear explanatory variable
  or as a linear and quadratic explanatory variable
  if this variable changes between discrete intervals.
  
  \code{urbinEffCat()} calculates 
  the effect of a categorical variable
  that is used as an explanatory variable,
  particularly if one wants to change the reference category and/or
  wants to group some of the categories together to a new category.

  The semi-elasticities calculated 
  by \code{urbinEla()} and \code{urbinElaInt()} indicate
  by how many percentage points the probability
  that the dependent variable has a value of one
  increases
  if the explanatory variable of interest
  increases by one percent.

  The effects calculated
  by \code{urbinEffInt()} and \code{urbinEffCat()} indicate
  by how much the probability
  that the dependent variable has a value of one
  increases
  if the explanatory variable of interest
  changes from the 'reference' interval/category 
  to a selected interval/category of interest
  (this effect multiplied by 100
  indicates the increase in percentage points).

  The four functions apply the Delta-method
  to calculate the approximate standard errors 
  of the calculated semi-elasticities and effects.
}

\usage{
urbinEla( allCoef, allXVal, xPos, model, allCoefVcov = NULL,
  seSimplify = !is.matrix( allCoefVcov ), xMeanSd = NULL, 
  iPos = 1, yCat = NULL )
  
urbinElaInt( allCoef, allXVal, xPos, xBound, model, 
  allCoefVcov = NULL, iPos = 1, yCat = NULL )
  
urbinEffInt( allCoef, allXVal = NULL, xPos, refBound, intBound, model,
  allCoefVcov = NULL, xMeanSd = NULL, iPos = 1, yCat = NULL )
  
urbinEffCat( allCoef, allXVal, xPos, xGroups, model, 
  allCoefVcov = NULL, iPos = 1, yCat = NULL )
}

\arguments{
  \item{allCoef}{a vector of all estimated coefficients
    (including intercept).}
  \item{allXVal}{a vector of the values of the explanatory variables,
    at which the semi-elasticity or the effect should be calculated
    (including intercept; the order of its elements must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{xPos}{a vector of non-negative integers 
    indicating the position(s) of the coefficient(s) and value(s)
    of the explanatory variable of interest
    in arguments \code{allCoef}, \code{allXVal}, 
    and eventually \code{allCoefVcov}
    (for \code{urbinElaInt()} and \code{urbinEffCat()}
    see \bold{Details} section below).}
  \item{xBound}{a numeric vector
    indicating the boundaries of the intervals.}
  \item{refBound}{a numeric vector of two elements
    that indicate the lower boundary and the upper boundary
    of the 'reference' interval, respectively.}
  \item{intBound}{a numeric vector of two elements
    that indicate the lower boundary and the upper boundary
    of the interval of interest, respectively.}
  \item{xGroups}{a vector consisting of \eqn{-1}s, \eqn{0}s, and \eqn{1}s, 
    where a \eqn{-1} indicates 
    that the category should belong to the new reference category, 
    a \eqn{1} indicates 
    that the category should belong to the new category of interest,
    and a zero indicates 
    that the category belongs to neither;
    the length of this vector must be equal to the number of categories;
    the last element of this vector corresponds to the reference category,
    while all other elements correspond to the categories 
    as indicated by argument \code{xPos}.}
  \item{model}{a character string that indicates the type of model 
    that was estimated to obtain the parameter estimates:
    \code{"lpm"} = linear probability model (see \bold{Details} section below), 
    \code{"probit"} = binary probit model or multivariate probit model
      (for multivariate probit models see \bold{Details} section below),
    \code{"oprobit"} = ordered probit model (see \bold{Details} section below), 
    \code{"logit"} = binary logit model, or
    \code{"mlogit"} = multinomial logit model (see \bold{Details} section below).}
  \item{allCoefVcov}{an optional argument
    that can be the variance-covariance matrix of all estimated coefficients
    OR the vector of the standard errors of all estimated coefficients
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{seSimplify}{logical value indicating whether the standard errors
    of the semi-elasticities should be calculated by a simplified equation
    (i.e., assuming that \eqn{\phi( \beta' x )} does not depend on \eqn{\beta}),
    which may be more appropriate 
    if the off-diagonal elements of the variance-covariance matrix are unknown.}
  \item{xMeanSd}{a vector with two elements: 
    the mean value and the standard deviation of the explanatory variable 
    of interest.
    If argument \code{allCoefVcov} is a vector of standard errors
    and the explanatory variables include a linear and a quadratic term
    of the variable of interest,
    the information provided in argument \code{xMeanSd} 
    is used to approximate the covariance 
    between the coefficient of the linear term 
    and the coefficient of the quadratic term.}
  \item{iPos}{a posivite integer 
    indicating the position of the intercept
    in arguments \code{allCoef}, \code{allXVal}, 
    and eventually \code{allCoefVcov}.
    A value of zero indicates that the model does not have an intercept.}
  \item{yCat}{a non-negative integer or a vector of non-negative integers
    that indicate(s) which category/categories of the dependent variable
    of a multinomial logit model should indicate a binary outcome of one.
    A zero indicates the base category of the dependent variable
    (see \bold{Details} section below).}
}

\details{
\bold{Argument \code{xPos} of \code{urbinElaInt()}}
must be a vector of non-negative integers
with length equal to the number of intervals.
Each element of this vector must refer to one interval,
whereas these intervals must be in ascending order.
Each element of this vector indicates 
the position of the coefficient of the respective dummy variable 
(i.e., the dummy variable that indicates 
whether the values of the explanatory variable of interest lie
in the corresponding interval or not)
in arguments \code{allCoef} 
and the position of the value of this dummy variable
(i.e., share of observations that lie in the respective interval
in the sample of interest)
in argument \code{allXVal}.
The position of the reference interval
(i.e., the interval without a corresponding dummy variable
as explanatory variable)
must be indicated by a zero.

\bold{Argument \code{xPos} of \code{urbinEffCat()}}
must be a vector of posivite integers 
that indicates the positions of the coefficients of the categorical variable
of interest in argument \code{allCoef}
and, equally, the positions of the shares of the non-reference categories 
in argument \code{allXVal}.
This vector must have one element less than the number of categories,
because argument \code{allCoef} does not include a coefficient 
of the reference category.

\bold{\code{urbinEffInt()} ignores argument \code{allXVal}}
in case of linear probability models.
In case of all other types of models, 
\code{urbinEffInt()} ignores the element(s) of argument \code{allXVal} 
that correspond to the variable of interest;
these values can be set to \code{NA}.

\bold{Linear Probability Model}:
If a user wants to calculate the semi-elasticity or effect
of a variable of interest
based on a linear probability model
using \code{urbin*( ..., model = "lpm")},
he/she can include only the coefficient(s)
of the explanatory variable of interest in argument \code{allCoef}
and omit all other coefficients.
In this case, arguments \code{allXVal} and \code{allCoefVcov}
must be adjusted accordingly.
When using \code{urbinEla()}, argument \code{allXVal} can be a scalar 
with the value of the variable of interest
even if the model includes both a linear and a quadratic term.
This simplified use of the \code{urbin*()} functions
is only unsuitable for linear probability models
and unsuitable for all other types of models.

\bold{Ordered Probit Model}:
The dependent variable has \eqn{P>2} distinct and strictly ordered categories.
In order to aggregate these \eqn{P} cateories into a binary outcome,
a user needs to choose a number \eqn{p \in \{1, \ldots, P-1\}}
that separates the \eqn{P} ordered categories 
into two mutually exclusive combined categories,
where the lower \eqn{p} categories indicate a binary outcome of zero
and the upper \eqn{P-p} categories indicate a binary outcome of one.
An ordered probit model does not have an intercept
but it has \eqn{P-1} thresholds between the \eqn{P} ordered categories.
If a user wants to calculate the semi-elasticity 
based on an ordered probit model 
using \code{urbinEla( ..., model = "oprobit")},
argument \code{allCoef} must include 
both the coefficients and the thresholds.
Consequently, argument \code{allCoefVcov} 
must also include the variences and covariances
or the standard arrors of both the coefficients and the thresholds.
The vector specified by argument \code{allXVal} 
must include elements that correspond to the thresholds,
where the element that corresponds to the \eqn{p}'s threshold
must be equal to \eqn{-1}, 
while the elements that correspond to the other thresholds
must be equal to \eqn{0}.
Argument \code{iPos} must specify the position of the \eqn{p}'s threshold
in arguments \code{allCoef}, \code{allXVal}, 
and eventually \code{allCoefVcov}.
While the \eqn{p}'s threshold must be specified in argument \code{allCoef}, 
all other thresholds can be omitted.
In this case, arguments \code{allXVal} and \code{allCoefVcov}
must be adjusted accordingly.

\bold{Multivariate Probit Model}:
This type of model has \eqn{P \geq 2}{P>=2} dependent variables.
The \code{urbin*()} functions
can calculate unconditional (but not conditional) semi-elasticities and effects
based on the estimation results from multivariate probit models.
Argument \code{allCoef} must include all coefficients
of the regression equation of the dependent variable of interest,
while the coefficients of the regression equations 
of the other dependent variables as well as the correlation parameters
must be omitted.
Argument \code{allCoefVcov} must be specified accordingly.

\bold{Multinomial Logit Model}:
The dependent variable has \eqn{P>2} distinct categories
(without any logical order).
For each category of the dependent variable
(except for the base category),
the multinomial logit model has a different set of coefficients.
Argument \code{allCoef} must be a vector of all these coefficients
with the following ordering:
\code{c(} coefficients of all explanatory variables 
for the first category of the dependent variable,
coefficients of all explanatory variables 
for the second category of the dependent variable, ...,
coefficients of all explanatory variables 
for the last category of the dependent variable \code{)},
where the coefficients of the base category
(that are all normalized to zero)
should not be included in argument \code{allCoef}.
The order of the coefficients for each and every category
must be the same as the order of the corresponding explanatory variables
in argument \code{allXVal}.
The order of the coefficients in argument \code{allCoefCov}
must be the same as the order of the coefficients in argument \code{allCoef}.
In order to aggregate the \eqn{P} cateories of the dependent variable 
into a binary outcome,
a user needs to use argument \code{yCat} to select 
\eqn{p \in \{1, \ldots, P-1\}} categories
that should indicate a binary outcome of one,
where a zero indicates the base category.
All other categories,
i.e., all categories that are not specified by argument \code{yCat}, 
indicate a binary outcome of zero.

\bold{Marginal Effects}:
If you know the marginal effects of the explanatory variable(s) of interest
on the dependent variable,
you can use these marginal effects 
as if they were coefficients of a linear probability model 
and calculate approximate semi-elasticities and effects
with \code{urbin*( ..., model = "lpm" )}
with argument \code{allCoef} equal to the marginal effect(s)
and eventually argument \code{allCoefVcov} 
equal to the variance covariance matrix or the standard errors
of the marginal effects.
}

\value{
The four \code{urbin*()} functions return a list of class \code{"urbin"} 
that contains the following components:
\item{call}{the matched call.}
\item{semEla}{the calculated semi-elasticity 
  (only \code{urbinEla()} and \code{urbinElaInt()}).}
\item{effect}{the calculated effect 
  (only \code{urbinEffInt()} and \code{urbinEffCat()}).}
\item{stdEr}{the approximate standard error 
  of the calculated semi-elasticity or effect.}
\item{allCoefVcov}{the variance covariance matrix 
  that was used to calculate the approximate standard error 
  of the calculated semi-elasticity or effect.}
\item{derivCoef}{the partial derivatives of the semi-elasticity or effect 
  with respect to the coefficients 
  that was used to calculate the approximate standard error 
  of the calculated semi-elasticity or effect, respectively.}
}

\examples{
# load data set
data( "Mroz87", package = "sampleSelection" )

# create dummy variable for kids
Mroz87$kids <- as.numeric( Mroz87$kids5 > 0 | Mroz87$kids618 > 0 )

# estimate probit model with linear and quadratic age terms
estProbitQuad <- glm( lfp ~ kids + age + I(age^2) + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitQuad )

# mean values of the explanatory variables
xMeanQuad <- c( 1, mean( Mroz87$kids ), mean( Mroz87$age ), 
  mean( Mroz87$age )^2, mean( Mroz87$educ ) )

# create dummy variables for age intervals
Mroz87$age30.37 <- Mroz87$age >= 30 & Mroz87$age <= 37
Mroz87$age38.44 <- Mroz87$age >= 38 & Mroz87$age <= 44
Mroz87$age45.52 <- Mroz87$age >= 45 & Mroz87$age <= 52
Mroz87$age53.60 <- Mroz87$age >= 53 & Mroz87$age <= 60

# probit estimation with age as interval-coded explanatory variable
estProbitInt <- glm( lfp ~ kids + age30.37 + age38.44 + age53.60 + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitInt )

# mean values of the explanatory variables
xMeanInt <- c( 1, colMeans( Mroz87[ , 
  c( "kids", "age30.37", "age38.44", "age53.60", "educ" ) ] ) )


##################################################################  
#######################  urbinEla()  #############################
##################################################################  
# semi-elasticity of age (full covariance matrix of coefficients)
urbinEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), model = "probit",
  vcov( estProbitQuad ) )

# semi-elasticity of age (standard errors of coefficients,
# mean and standard deviation of age provided, 
# simplified calculation of standard error)
urbinEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), model = "probit", 
  sqrt( diag( vcov( estProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
  
  
##################################################################  
#####################  urbinElaInt()  ############################
##################################################################  
# semi-elasticity of age (full covariance matrix of coefficients)
urbinElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), 
  model = "probit", vcov( estProbitInt ) )

# semi-elasticity of age (only standard errors of coefficients)
urbinElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), 
  model = "probit", sqrt( diag( vcov( estProbitInt ) ) ) )
 
 
##################################################################  
#####################  urbinEffInt()  ############################
##################################################################  
# effect of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix of coefficients)
urbinEffInt( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "probit",
  vcov( estProbitQuad ) )
  
# effect of age changing from the 30-40 interval to the 50-60 interval
# (with standard errors of coefficients as well as 
# mean and standard deviation of age) 
urbinEffInt( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "probit", 
  sqrt( diag( vcov( estProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
  
  
##################################################################  
#####################  urbinEffCat()  ############################
##################################################################  
# effect of age changing from 38-52 years (2nd category + reference category) 
# to 53-60 years (3rd category) (with full covariance matrix)
urbinEffCat( coef( estProbitInt ), xMeanInt, c( 3:5 ), c( 0, -1, 1, -1 ),
  model = "probit", vcov( estProbitInt ) )

# effect of age changing from 38-52 years (2nd category + reference category) 
# to 53-60 years (3rd category) (with standard errors only)
urbinEffCat( coef( estProbitInt ), xMeanInt, c( 3:5 ), c( 0, -1, 1, -1 ),
  model = "probit", sqrt( diag( vcov( estProbitInt ) ) ) )
}

\keyword{models}
\keyword{regression}
