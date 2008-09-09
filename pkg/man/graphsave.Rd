\name{graphsave}
\alias{graphsave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graphsave(path, xynames, suffix, prompt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{ ~~Describe \code{path} here~~ }
  \item{xynames}{ ~~Describe \code{xynames} here~~ }
  \item{suffix}{ ~~Describe \code{suffix} here~~ }
  \item{prompt}{ ~~Describe \code{prompt} here~~ }
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
function(path,xynames,suffix,prompt){
	if(prompt<2){return();}
	fname<-paste(path,paste(xynames,collapse=""),suffix,sep="");
	message<-paste("Do you wish to save this graph as",fname,"?\n");
	choices<-c("Save under suggested name.",
		   "Give a different name (point and click).",
		   "Type in a name from the console.",
		   "Don\'t save.")
	cat(message);
	input<-menu(choices,graphics=T,title=message);
	switch(input,
		dev.copy2eps(file=fname),
		{fname<-tclvalue(tkgetSaveFile(initialdir=path,initialfile=paste(xynames,collapse=""),
			defaultextension=suffix));dev.copy2eps(file=fname);},
		{res<-0;class(res)<-'try-error';
		 while(class(res)=='try-error'){
			fname<-readline("Please type the path and name to which you want to save the file: ");
			print(fname); res<-try(dev.copy2eps(file=fname));
			if(class(res)=='try-error'){cat('Oops. Bad file name or path, try again please.\n')}}
		},
		cat("\nContinuing without saving.")
	)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
