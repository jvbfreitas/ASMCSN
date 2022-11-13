\name{qq.ASMCSN}
\alias{qq.ASMCSN}
\title{
  Compute QQ plot for standardized residuals of ASMCSN models.
}
\description{
  This function can be used to plot the QQ of stadardized residuals for ASMCSN models.
}
\usage{
  qq.ASMCSN(obj, confi = 0.95, repl = 100)
}
\arguments{
  \item{obj}{
    An R object returned by ASMCSN.
  }
  \item{confi}{
    Confidence level.
}
\item{repl}{
    Number of replicas to construct the simulated envelopes.
}
}
\details{
}
\value{
}
\references{
}
\author{
João de Freitas, jvbfreitas@ime.unicamp.br.
}
\note{
}

\seealso{
}
\examples{
library(SemiPar)
library(ASMCSN)
data("ragweed")
data = ragweed
data$y = sqrt(ragweed$ragweed)
auxm = ASMCSN(y ~ rain + wind.speed - 1, ~ day.in.seas + temperature, data = data, mixed = "1", nknot = 30, a0 = c(100,100), a1 = c(5000,10000))
qq.ASMCSN(auxm)
}
