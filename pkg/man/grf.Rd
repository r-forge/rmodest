\name{grf}
\alias{grf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grf(par, cons = rep.int(1, 4), x, y = NULL, keep = 1:2, np = 2, model = "g")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{ ~~Describe \code{par} here~~ }
  \item{cons}{ ~~Describe \code{cons} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{keep}{ ~~Describe \code{keep} here~~ }
  \item{np}{ ~~Describe \code{np} here~~ }
  \item{model}{ ~~Describe \code{model} here~~ }
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
function(par,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,model='g'){
 par1<-par2<-rep.int(NA,4); par1[keep]<-par2[keep]<-par[1:np]; cons<-cons[keep];
 par2[keep][cons>0]<-par[(np+1):length(par)];
 env1<-list(a=par1[1],b=par1[2],c=par1[3],s=par1[4],x=x); 
 env2<-list(a=par2[1],b=par2[2],c=par2[3],s=par2[4],x=y);
 out1<-out2<-c();
 for(i in dex[[model]]){
 	out1<-c(out1,sum(eval(i,envir=env1)));out2<-c(out2,sum(eval(i,envir=env2)));
 }
 out1<-c(out1,out1); out2<-c(out2,out2);
 id1<-c(rep.int(1,np),1-cons); id2<-c(1-cons,rep.int(1,np));
 out<-(id1*out1)+(id2*out2); 
 if(rm.temp$debug){cat("gradient:");print(out);}
 out[c(1:np,(np+1:np)[as.logical(cons)])];
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
