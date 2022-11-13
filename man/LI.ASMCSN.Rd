\name{LI.ASMCSN}
\alias{LI.ASMCSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute local influence measures for ASMCSN model
}
\description{
  This function can be used to compute the local influence measures under certain perturbations for ASMCSN models.
}
\usage{
LI.ASMCSN(obj, type = c("CW","SP,RP","SKP","CP"), xc = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    an R object returned by ASMCSN.
  }
  \item{type}{
    a character string specifying the perturbation scheme. The following are permitted: "CW":"case-weight", "SP":"scale", "RP":"response", "SKP":"skewness", "CP":"covariate".
  }
  \item{xc}{
    position of the covariate of interest in the specification matrix if the covariate perturbation scheme is chosen.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
\references{
}
\author{
  João de Freitas, jvbfreitas@ime.unicamp.br.
}
\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
  }
\examples{
library(SemiPar)
library(ASMCSN)
data("ragweed")
data = ragweed
data$y = sqrt(ragweed$ragweed)
auxm = ASMCSN(y ~ rain + wind.speed - 1, ~ day.in.seas + temperature, data = data, mixed = "1", nknot = 30, a0 = c(100,100), a1 = c(5000,10000))
##### CW perturbation
lic = LI.ASMCSN(auxm, type = "CW")
ggplot() + geom_point(aes(x = 1:n, y = lic)) +
 xlab("Index") + ylab(expression(D[max])) +
 theme(
   panel.border = element_blank(),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   panel.background = element_blank(),
   axis.line = element_line(colour = "grey"),
   text=element_text(size=15)
 )
}
