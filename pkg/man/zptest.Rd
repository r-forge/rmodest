\name{zptest}
\alias{zptest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zptest(x1, x2, n1, n2, n, zponly = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{ ~~Describe \code{x1} here~~ }
  \item{x2}{ ~~Describe \code{x2} here~~ }
  \item{n1}{ ~~Describe \code{n1} here~~ }
  \item{n2}{ ~~Describe \code{n2} here~~ }
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{zponly}{ ~~Describe \code{zponly} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x1,x2,n1,n2,n,zponly=F){
	that1<-x1/n1; that2<-x2/n2;
	ttil<-(x1+x2)/n;
	denom<-sqrt(ttil*(1 - ttil)*(1/n1 + 1/n2));
	zp<-(that1 - that2)/denom;
	if(zponly){return(zp);}else{
		output<-cbind(x1,that1,x2,that2,ttil,zp);
		return(output);
	}
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
