\name{diag.ASMCSN}
\alias{diag.ASMCSN}
\title{
  Compute diagnostics measures for ASMCSN model
}
\description{
  This function can be used to compute some of the regression diagnostics for ASMCSN models.
}
\usage{
  diag.ASMCSN(obj, type = c("Gcooks", "Ghat"))
}
\arguments{
  \item{obj}{
    An R object returned by ASMCSN.
  }
  \item{type}{
    Type of measure: generalized Cook's distance or the generalized hat matrix.
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
GCD = diag.ASMCSN(auxm, type = "Gcooks")

n = length(GCD)
dadaux30 = data.frame(id = 1:n, gcd = GCD/sum(GCD))
ggplot(dadaux30,aes(x = id, y = gcd)) + geom_point() +
 xlab("Index") + ylab("Normalized generalized Cook's distance") +
 geom_text_repel(data = dadaux30[dadaux30$gcd>0.002987,], aes(label = id),
                 box.padding = 1, show.legend = FALSE) +
 theme(
   panel.border = element_blank(),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   panel.background = element_blank(),
   axis.line = element_line(colour = "grey"),
   text=element_text(size=15)
 )
  }
