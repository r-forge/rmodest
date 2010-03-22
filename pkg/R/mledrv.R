mledrv<-function(x,y=NULL,model='g',pars,what='hess'){
	mydrv<-get(paste('d',model,sep=''));
	if(!is.null(y)){return(pars[]<-NA);}
	out<-mydrv(pars,x);
	switch(what,
		hess = return(tryCatch(apply(attributes(out)$hessian,c(2,3),sum),
			      error=function(e){print(e);NA;})),
		sd = return(tryCatch(sqrt(abs(diag(solve(apply(attributes(out)$hessian,c(2,3),sum))))),
			      error=function(e){print(e);NA;})),
		grad = return(tryCatch(apply(attributes(out)$gradient,2,sum),error=function(e){print(e);NA;}))
	)
}