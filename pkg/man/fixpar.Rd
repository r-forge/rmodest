\name{fixpar}
\alias{fixpar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fixpar(par, np, keep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{ ~~Describe \code{par} here~~ }
  \item{np}{ ~~Describe \code{np} here~~ }
  \item{keep}{ ~~Describe \code{keep} here~~ }
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
function(par,np,keep){
	lp<-length(par);cons<-rm.temp$cons;
	if(is.null(rm.temp$y)){
		if(lp>np){
			warning("You gave more parameters than necessary
for a single-model fit. Only the first ",np," parameters were used.");
			return(par[1:np]);
		}
		if(lp==np){return(par);}
		if(lp==1){return(rep(par,np));}
	}
	if(lp==1){return(rep(par,2*np)[c(rep.int(T,np),cons)]);}
	if(lp==np+sum(cons)){return(par);}
	if(lp==np){return(c(par,par)[c(rep.int(T,np),cons)]);}
	if(lp==2*np){return(par[c(rep.int(T,np),cons)]);}
	if(lp==4){return(fixpar(c(par,par)[keep],np,keep));}
	if(lp==8){return(fixpar(c(par)[keep],np,keep));}
	stop("If starting parameters (par) are specified,
they should be a numerical vector equal in length either
to 4, or 8, or the number of parameters in your model, or
twice the number of parameters in your model, or the
total number of *unique* parameters in your model (i.e.
two starting values for each unconstrained parameter and
one starting value for each constrained one).");
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
