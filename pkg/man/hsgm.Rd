\name{hsgm}
\alias{hsgm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hsgm(par, cons = NULL, x = NULL, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{ ~~Describe \code{par} here~~ }
  \item{cons}{ ~~Describe \code{cons} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
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
function(par,cons=NULL,x=NULL,y=NULL){
	if(is.null(cons)) cons<-rm.temp$cons;
	if(is.null(x)) { x<-rm.temp$x; y<-rm.temp$y; }
	par2<-par[1:3]; par2[cons>0]<-par[4:length(par)];
	env1<-list(a=par[1],b=par[2],c=par[3],x=rm.temp$x);
	env2<-list(a=par2[1],b=par2[2],c=par2[3],x=rm.temp$y);
	out1a<-c(sum(eval(D(D(exgm,'a'),'a'),envir=env1)),
		sum(eval(D(D(exgm,'b'),'a'),envir=env1)),
		sum(eval(D(D(exgm,'c'),'a'),envir=env1)));
	out1b<-c(sum(eval(D(D(exgm,'a'),'b'),envir=env1)),
		sum(eval(D(D(exgm,'b'),'b'),envir=env1)),
		sum(eval(D(D(exgm,'c'),'b'),envir=env1)));
	out1c<-c(sum(eval(D(D(exgm,'a'),'c'),envir=env1)),
		sum(eval(D(D(exgm,'b'),'c'),envir=env1)),
		sum(eval(D(D(exgm,'c'),'c'),envir=env1)));
	out1<-rbind(out1a,out1b,out1c); 
	out1<-rbind(cbind(out1,out1),cbind(out1,out1));
	out2a<-c(sum(eval(D(D(exgm,'a'),'a'),envir=env2)),
		sum(eval(D(D(exgm,'b'),'a'),envir=env2)),
		sum(eval(D(D(exgm,'c'),'a'),envir=env2)));
	out2b<-c(sum(eval(D(D(exgm,'a'),'b'),envir=env2)),
		sum(eval(D(D(exgm,'b'),'b'),envir=env2)),
		sum(eval(D(D(exgm,'c'),'b'),envir=env2)));
	out2c<-c(sum(eval(D(D(exgm,'a'),'c'),envir=env2)),
		sum(eval(D(D(exgm,'b'),'c'),envir=env2)),
		sum(eval(D(D(exgm,'c'),'c'),envir=env2)));
	out2<-rbind(out2a,out2b,out2c); 
	out2<-rbind(cbind(out2,out2),cbind(out2,out2));
	id1<-c(1,1,1,1-cons); id2<-c(1-cons,1,1,1);
	id1<-id1\%o\%id1; id2<-id2\%o\%id2;
	out<-out1*id1+out2*id2;
	if(rm.temp$debug){print("hessian:");print(out);}
	out[c(1:3,4:6*cons),c(1:3,4:6*cons)];
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
