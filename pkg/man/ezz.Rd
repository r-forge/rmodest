\name{ezz}
\Rdversion{1.1}
\alias{ezz}
\title{
Score test for survival data.
}
\description{
This function takes survival data from two groups and a list of quantiles and returns a table indicating whether one group or the other is significantly overrepresented in each quantile.
}
\usage{
ezz(d, d2 = NULL, g, quant = 0.9, step = 0.01, thresh = 0.05, gnames = NULL, noob = F, thetas = c(0.05, 0.95), debug = F)
}
\arguments{
  \item{d}{
An integer vector of survival times (usually days).
}
  \item{d2}{
An optional second integer vector of survival times (usually days). If this argument is used, 'd' is assumed to be the control group and 'd2' is assumed to be the experimental group. If this argument is left at its default of NULL, both groups are assumed to be in 'd' and the 'g' argument is used to separate them out.
}
  \item{g}{
A charcter vector containing two distinct group names, and equal in length to the combined length of 'd' and 'd2'. However, if both 'd' and 'd2' are used, this argument may also just have a length of 2.
}
  \item{quant}{
The quantiles to compare. The default is 0.9, but this argument can take a vector of numbers between 0 and 1.
}
  \item{step}{
Internal parameter used by score algorithm.
}
  \item{thresh}{
Significance threshold.
}
  \item{gnames}{
A character vector of length two indicating the group names as they should appear in the output. If not used, the group names are obtained from the 'g' argument instead.
}
  \item{noob}{
If this argument is set to TRUE then a simplified version of the output is printed.
}
  \item{thetas}{
Internal parameter used by score algorithm.
}
  \item{debug}{
Debug flag. Not necessary for normal use.
}
}
\details{
This test divides the observed difference in survival fractions at each quantile (of the combined dataset) by the estimated standard error under the null hyopthesis. The resulting table indicates whether one group or the other is overrepresented in each quantile beyond a level that would be expected by random chance.
}
\value{
Returns a data.frame with quantiles as the row names and the following columns:
\item{age}{
Age at the given quantile.
}
\item{srv.A}{
Number of surviving animals in group 'A' at the given age (A is substituted with the control group's name).
}
\item{srv.A.frc}{
Fraction of surviving animals in group 'A' at the given age (A is substituted with the control group's name).
}
\item{srv.B}{
Number of surviving animals in group 'B' at the given age (B is substituted with the experimental group's name).
}
\item{srv.B.frc}{
Fraction of surviving animals in group 'B' at the given age (B is substituted with the experimental group's name).
}
\item{survtotal.frc}{
Fraction of total surviving animals.
}
\item{zp}{
Score test statistic.
}
\item{p}{
Significance level.
}
\item{sig}{
Indicator of significance.
}
}
\references{
Wang,C., Li,Q., Redden,D.T., Weindruch,R., Allison,D.B. (2004). Statistical methods for testing effects on "maximal lifespan". Mechanisms of Ageing and Development 125, 629-632

Mehotra,D.V., Chan,I.S.F., Berger,R.L. (2003). A cautionary note on exact unconditional inference for a difference between two independent binomial proportions. Biometrics 59, 441-450
}
\author{
Alex F. Bokov
}

\seealso{
\code{\link{survgui}}
}
\keyword{ survival }
