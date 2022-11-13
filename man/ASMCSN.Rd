\name{ASMCSN}
\alias{ASMCSN}
\title{
Fit Additive Partially Linear Scale Mixtures of Centered Skew-Normal models (ASMCSN)
}
\description{
The ASMCSN function fits additive partially linear models using the scale mixtures of centered skew-normal distributions and natural cubic splines.
}
\usage{
ASMCSN(formula.pp, formula.npp, data, mixed = "1", k = 1/3,
                  iter.max = 300, alpha.fix = F, alpha0, nknot, a0, a1)
}
\arguments{
\item{formula.pp}{
an object of class "formula" (or one that can be coerced to that class) with only parametric predictos: a expression of the form "response ~ parametric predictors". See the documentation of \code{\link{glm}} and \code{\link{formula}} for details.
}
\item{formula.npp}{
an object of class "formula" (or one that can be coerced to that class) with only nonparametric predictors: a expression of the form "~ nonparametric predictors". See the documentation of \code{\link{glm}} and \code{\link{formula}} for details.
}
\item{data}{
data frame or environment containing the variables in the formula.pp and formula.npp.
}
\item{mixed}{
a character string specifying the mixing measure related to the scale mixtures. The options are:
\itemize{
\item "1": which generates the centered skew-normal distribution (CSN).
\item "gamma1p": which generates the centered skew t distribution (CST).
\item "gamma2p": which generates the centered skew generalized t distribution (CSGT).
\item "beta": which generates the centered skew slash distribution (CSS).
\item "binary": which generates the centered skew contaminated normal distribution (CSCN).
\item "betaprime": which generates the centered beta prime normal distribution (CSBPN).
\item "birnbaum-saunders": which generates the centered skew Birnbaum-Saunders normal distribution (CSBSN).
\item "generalizedgamma": which generates the centered skew generalized gamma normal distribution (CSGGN).
}
}
\item{k}{
A constant related to the stochastic approximation that determines the percentage of initial iterations with no memory.
}
\item{iter.max}{
The number of iterations.
}
\item{alpha.fix}{
Should the smooth parameter (alpha) fixed? True for yes.
}
\item{alpha0}{
If alpha.fix==T, alpha0 is the value of alpha for all iterations. If alpha.fix==F, alpha0 is the value of alpha in the first iteration. If none alpha0 is given, so alpha0 = 1.
}
\item{a0,a1}{
Bounds on the variables for the "L-BFGS-B" method on BIC minimization.
}
}
\details{
It is important that the input data be sorted by the argument called "id" so that the repeated measures structure is correctly specified.
}
\value{
An object of class "ASMCSN" representing the fit. The 'coefficients' list contains the estimated parameters. The 'comp' list contains the specification matrix (X), incidence matrix (N) and the penalty matrix (K). The 'trajectory' list contains the estimated parameters in each iteration.
}
\references{
}
\author{
João de Freitas, jvbfreitas@ime.unicamp.br.
}

\examples{
library(SemiPar)
library(ASMCSN)
data("ragweed")
data = ragweed
data$y = sqrt(ragweed$ragweed)
auxm = ASMCSN(y ~ rain + wind.speed - 1, ~ day.in.seas + temperature, data = data, mixed = "1", nknot = 30, a0 = c(100,100), a1 = c(5000,10000))
auxm = ASMCSN(y ~ rain + wind.speed+temperature - 1, ~ day.in.seas, data = data, mixed = "beta", nknot = 30, a0 = c(100), a1 = c(5000))
}
