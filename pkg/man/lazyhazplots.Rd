\name{lazyhazplots}
\alias{lazyhazplots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lazyhazplots(ts, pars, ints, cols = rainbow(length(ts)), xlim = NULL, add = F, main = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ts}{ ~~Describe \code{ts} here~~ }
  \item{pars}{ ~~Describe \code{pars} here~~ }
  \item{ints}{ ~~Describe \code{ints} here~~ }
  \item{cols}{ ~~Describe \code{cols} here~~ }
  \item{xlim}{ ~~Describe \code{xlim} here~~ }
  \item{add}{ ~~Describe \code{add} here~~ }
  \item{main}{ ~~Describe \code{main} here~~ }
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
function(ts,pars,ints,cols=rainbow(length(ts)),xlim=NULL,add=F,main=""){
	mnt<-mnx<-mx<-c(); 
	for(j in ts){
		mx<-c(mx,max(j$ux,na.rm=T)); mnx<-c(mnx,min(j[j$ux>0,]$ux,na.rm=T));}
	mnt<-match(T,ts[[which.max(ints)]]$ux>0)
	mx<-max(mx); mnx<-min(mnx); if(is.null(xlim)){xrng<-c(mnt,max(ts[[which.min(ints)]]$time));}else{
		xrng<-c(10,xlim);}
	#print(xrng);
	if(add){points(1+ts[[1]]$time,ts[[1]]$ux,pch=3,col=cols[1]);} else {
		plot(1+ts[[1]]$time,ts[[1]]$ux,log='xy',type='p', pch=3,
			xlim=xrng,ylim=c(mnx,mx),xlab="Time",ylab="Hazard Rate",col=cols[1],main=main);
	} 
	lines(1:xrng[2],srvhaz(1:xrng[2],pars[1,1],pars[1,2],i=ints[1]),col=cols[1],lty=2);
	lines(1:xrng[2],srvhaz(1:xrng[2],pars[2,1],pars[2,2],pars[2,3],i=ints[1]),col=cols[1],lty=1);
	xrng<-xrng[1]:xrng[2];
	if(length(ts)>1){
		for(j in 2:length(ts)){
		points(1+ts[[j]]$time,ts[[j]]$ux,col=cols[j],pch=3); 
			lines(1:xrng[2],ints[j]*srvhaz(1:xrng[2],pars[1,1],pars[1,2],i=ints[j]),col=cols[j],lty=2);
			lines(1:xrng[2],ints[j]*srvhaz(1:xrng[2],pars[2,1],pars[2,2],pars[2,3],i=ints[j]),col=cols[j],lty=1);
		}
	}
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
