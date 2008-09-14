ltry<-try(library(survomatic));
if(class(ltry)=='try-error'){
	install.packages('survomatic',repos='http://r-forge.r-project.org');}
ltry<-try(library(survomatic));
if(class(ltry)=='try-error'){
	cat('\nProblem installing library. Functions may not work but
your data is still safe. See if you have a problem with
your internet connection. If it\'s working properly and
you still get this message, contact
alex.bokov@gmail.com\n');} else {go();}
