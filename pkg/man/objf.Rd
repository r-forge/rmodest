\name{objf}
\alias{objf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
objf(par, lb = 0, ub = 1, cons = rep.int(1, 4), x, y = NULL, keep = 1:2, np = 2, ex = exg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{ ~~Describe \code{par} here~~ }
  \item{lb}{ ~~Describe \code{lb} here~~ }
  \item{ub}{ ~~Describe \code{ub} here~~ }
  \item{cons}{ ~~Describe \code{cons} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{keep}{ ~~Describe \code{keep} here~~ }
  \item{np}{ ~~Describe \code{np} here~~ }
  \item{ex}{ ~~Describe \code{ex} here~~ }
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
function(par,lb=0,ub=1,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,ex=exg){
	if(max(par<=lb)!=0|max(par>=ub)!=0){return(NA);}
	par1<-par2<-rep.int(NA,4); par1[keep]<-par2[keep]<-par[1:np];
	par2[keep][cons[keep]>0]<-par[(np+1):length(par)];
	env1<-list(a=par1[1],b=par1[2],c=par1[3],s=par1[4],x=x); 
	env2<-list(a=par2[1],b=par2[2],c=par2[3],s=par2[4],x=y);
	out<-sum(eval(ex,envir=env1))+sum(eval(ex,envir=env2));
	if(!is.na(out) & is.finite(out) & out < -50){return(out);}else{return(NA);}
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
