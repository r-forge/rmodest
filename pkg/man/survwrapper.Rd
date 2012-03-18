\name{survwrapper}
\alias{survwrapper}
\title{
Function for fitting and selecting hazard functions to survival data.
}
\description{
Takes one or two vectors of event times (numeric format) and optionally
corresponding vectors of indicator variables to designate right-censored
events. Fits several mortality models, selects the best fitting one/s,
and if two vectors were given, tests hypotheses about the model parameters.
}
\usage{
survwrapper(x, y = NULL, models = c("g", "gm", "l", "lm"), cx = rep(1, length(x)), cy = rep(1, length(y)), ext = F, n = length(c(x, y)), AIC = F, BIC = F, breakties = "AIC", compare.matrix = NULL, constraint.matrix = NULL, thresh = 0.05, smooth = 7)
}
\arguments{
  \item{x}{
A numeric vector of event times. For example, number of days an
individual has survived. 
}
  \item{y}{
An optional second numeric vector of event times, in the same units as
\code{x}. Specifying this vector enables fitting joint models and
hypothesis testing.
}
  \item{models}{
A character vector of model names: \code{g} = Gompertz, \code{gm} =
Gompertz-Makeham, \code{l} = logistic, \code{lm} = logistic-Makeham, and
\code{w} = Weibull. The default is to fit all the models except Weibull.
}
  \item{cx}{
A vector of 0 and 1 the same length as \code{x}. Ordinary events
(e.g. death) are indicated by 1 and censored events (e.g. withdrawal
from the study) are indicated by 0.
}
  \item{cy}{
A vector of 0 and 1 the same length as \code{y}. Ordinary events
(e.g. death) are indicated by 1 and censored events (e.g. withdrawal
from the study) are indicated by 0.
}
  \item{ext}{
Not implemented.
}
  \item{n}{
Total sample size. Should normally be left for the script to
automatically calculate, but can be specified when survwrapper is called
from another script repeatedly in order to speed up runtimes.
}
  \item{AIC}{
Whether to calculate the AIC (Akaike Information Criterion) for each
candidate model.
}
  \item{BIC}{
Whether to calculate the BIC (Bayes Information Criterion) for each
candidate model.
}
  \item{breakties}{
What criterion to use for choosing a model if more than one is justified
by the comparisons.
}
  \item{compare.matrix}{
A matrix for specifying a customized comparison algorithm.
}
  \item{constraint.matrix}{
A matrix of 1's and 0's for specifying a customized set of parameter
constraints to test.
}
  \item{thresh}{
Significance cutoff.
}
  \item{smooth}{
Not yet supported.
}
}
\details{
In progress.
}
\value{
\item{x.m}{
}
\item{y.m}{
}
\item{xy.sm}{
}
\item{par.differences}{
}
\item{x}{
}
\item{y}{
}
\item{cx}{
}
\item{cy}{
}
\item{x.d}{
}
\item{y.d}{
}
\item{suggested.models}{
}
\item{nx}{
}
\item{ny}{
}
}
\references{
Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389
}
\author{
Alex F. Bokov (bokov@uthscsa.edu), Jon A. Gelfond
}
\note{
Uses Nelder-Mead algorithm to find maximum likelihood estimates of model
parameters.
}

\seealso{
  \code{\link{optim}}
}
\examples{
## Generate two sets of survival times.
population1 <- simsurv(629,type='g',p=c(7.33e-4,0.1227,0,0));
population2 <-
  simsurv(574,type='lm',p=c(5.4818e-5,0.1543,0.0023,0.6018));

## Fit models to the populations and compare the parameters.
models1vs2 <- survwrapper(population1,population2);
}
\keyword{ survival }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
