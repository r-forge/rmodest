\name{go}
\alias{go}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  The main function of Survomatic, which takes the raw survival data and performs analysis on it.
}
\usage{
go(x,y,xynames=c(),units="d",save=T,path,prompt=2,xlim=1350,slwd=4,
scol=c('darkred',1),spch=24:25,spcex=0,splwd=slwd,spbg=scol,sleg=T,sxcex=1.5,sxlwd=3,
qstart=1,demint=30, qrse='boot')}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of ages at death for control group.}
  \item{y}{Vector of ages at death for experimental group.}
  \item{xynames}{Character vector with labels for the two groups.}
  \item{units}{Not currently used.}
  \item{save}{Boolean value of whether to save the various output files or not.}
  \item{path}{Path for where to save the output files.}
  \item{prompt}{How aggressively to prompt you to make choices, versus deciding on your behalf.}
  \item{xlim}{The length of the x axis for the Kaplan-Meier curves.}
  \item{slwd}{Line width for the Kaplan-Meier survival plot.}
  \item{scol}{Colors for the plot lines.} 
  \item{spch}{Shapes for the Kaplan-Meier survival plot} 
  \item{spcex}{Size of the shapes for the Kaplan-Meier survival plot.}
  \item{splwd}{Line width for the shapes on the Kaplan-Meier survival plot.}
  \item{spbg}{Fill color for the shapes on the Kaplan-Meier survival plot.}
  \item{sleg}{Whether or not to plot a legend for the Kaplan-Meier survival plot.}
  \item{sxcex}{Font size for all the plot axes.}
  \item{sxlwd}{Line width for all the plot axes.}
  \item{qstart}{Which set of quantiles to try first.}
  \item{demint}{Default intervals to break up the hazard function into for plotting purposes.}
  \item{qrse}{Name of algorithm to use for calculating confidence intervals for quantile regression plot.}
}
\details{
  This is intended to be a wrapper that automagically calls various other 
  routines on your survival data. At the moment these
}
\value{
  A list object with a ton of potentially useful output, in the process
  of being documented.
  If it is a LIST, use
  \item{report}{A table listing some summary statistics and the quantiles flagged by quantile regression and the score test as potentially significant.}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ Alex F. Bokov }
\note{ 
This package does NOT yet support censoring. It assumes that all the ages
entered correspond to actual death events. Censoring is important and on the
short list of features that will soon be implemented. Also, be aware that
currently you are unlikely to observe significant differences in mortality
parameters for sample sizes smaller than 100. The mortality modeling portion
of this package is under active development and is likely to change often.
Please make sure you are using the latest revision of Survomatic.
}
\seealso{ \code{\link{survfit}}, \code{\link{rq}}, \code{\link{opsurv}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ survival }