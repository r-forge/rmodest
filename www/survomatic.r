ltry<-try(library(Survomatic));
if(class(ltry)=='try-error'){
	install.packages('Survomatic',,c('http://www.rforge.net/','http://cran.r-project.org/');}
ltry<-try(library(Survomatic));
if(class(ltry)=='try-error'){
	cat('\nProblem installing library. Functions may not work but
your data is still safe. See if you have a problem with
your internet connection. If it\'s working properly and
you still get this message, contact
alex.bokov@gmail.com\n');} else {go();}
