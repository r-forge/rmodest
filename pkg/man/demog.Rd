\name{demog}
\alias{demog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
demog(t, int = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ ~~Describe \code{t} here~~ }
  \item{int}{ ~~Describe \code{int} here~~ }
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
function(t,int=1){
	if(int>1){t<-round(t/int);}
	t<-sort(t); t.rle<-rle(t); d<-px<-lx<-c();
	ts<-1:max(t.rle$values);
	for(i in ts){
		
		j<-match(i,t.rle$values);
		if(is.na(j)){d<-c(d,0);}else{d<-c(d,t.rle$lengths[j]);}
	}
	for(i in ts) {
		lx<-c(lx,(sum(d)-sum(d[1:i]))/sum(d));
	}
	for(i in ts) {px<-c(px,lx[i+1]/lx[i]);}
	ux<- -log(px); ux[is.infinite(ux)]<-NA; lnux<-log(ux); lnux[is.infinite(lnux)]<-NA;
	out<-cbind(ts,d,lx,px,ux,lnux);
	colnames(out)<-c("time","deaths","lx","px","ux","lnux");
	data.frame(out);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
