\name{ezz}
\alias{ezz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ezz(d, d2 = NULL, g, quant = 0.9, step = 0.01, thresh = 0.05, gnames = NULL, noob = F, thetas = c(0.05, 0.95), debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ ~~Describe \code{d} here~~ }
  \item{d2}{ ~~Describe \code{d2} here~~ }
  \item{g}{ ~~Describe \code{g} here~~ }
  \item{quant}{ ~~Describe \code{quant} here~~ }
  \item{step}{ ~~Describe \code{step} here~~ }
  \item{thresh}{ ~~Describe \code{thresh} here~~ }
  \item{gnames}{ ~~Describe \code{gnames} here~~ }
  \item{noob}{ ~~Describe \code{noob} here~~ }
  \item{thetas}{ ~~Describe \code{thetas} here~~ }
  \item{debug}{ ~~Describe \code{debug} here~~ }
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
function(d,d2=NULL,g,quant=.9,step=.01,thresh=.05,
	      gnames=NULL,noob=F,thetas=c(.05,.95),debug=F){
	tstart<-proc.time()[3];
	if(is.null(d2)){
		groups<-levels(as.factor(g));
		d1<-subset(d,g==groups[1]);
		d2<-subset(d,g==groups[2]);
	} else { d1<-d; d<-c(d1,d2); if(is.null(g)){groups<-c('a','b');}else{groups<-g;}}
	#mintheta<-min(1/d1,1/d2);
	n1<-length(d1); n2<-length(d2); n<-n1+n2;
	output<-c(); q<-quantile(d,quant);
	if(debug) cat("\nRunning zptest on x1 and x2.");
	x1<-x2<-c();
	for(i in q){x1<-c(x1,sum(d1>i)); x2<-c(x2,sum(d2>i));}
	output<-zptest(x1,x2,n1,n2,n);
	output<-cbind(q,output);
	if(is.null(gnames)){g1<-groups[1];g2<-groups[2];}
	else{g1<-gnames[1];g2<-gnames[2];}
	surv1<-paste("srv",g1,sep='.');
	surv2<-paste("srv",g2,sep='.');
	surv1n1<-paste("srv",g1,"frc",sep='.');
	surv2n2<-paste("srv",g2,"frc",sep='.');
	colnames(output)<-c('age',surv1,surv1n1,surv2,surv2n2,'survtotal.frc','zp');
	output<-data.frame(output);
	output$zp[is.nan(output$zp)]<-9999; #print(output$zp);

	p<-zptests<-c();
	if(debug) cat("\nRunning zptest on k n1\'s and n2\'s.");
	zptests<-outer(0:n1,0:n2,FUN="zptest",n1=n1,n2=n2,n=n,zponly=T);
	zptests[is.nan(zptests)]<-9999;

	for(i in 1:dim(output)[1]){
		sup<-c(); cat('.');
		keep<-abs(zptests)>abs(output$zp[i]);
		for(j in seq(thetas[1],thetas[2],step)){
			if(debug) cat("\nRunning zpprob on n1\'s and n2\'s.");
			zptemp<-outer(0:n1,0:n2,FUN="zpprob",n1=n1,n2=n2,th=j);
			sup<-c(sup,sum(zptemp[keep]));
		}
		p<-c(p,max(sup));
	}
	sig<-p<thresh; sig[sig==T]<-"*"; sig[sig==F]<-"";
	output<-cbind(output,p,sig);
	rownames(output)<-quant;
	if(noob){output<-output[,c(1:5,8,9)];}
	print(proc.time()[3]-tstart);
	output;
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
