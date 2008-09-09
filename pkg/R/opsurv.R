`opsurv` <-
function(model,x,y=NULL,par=rep(1e-6,4),cons=1,old=F,usegr=T,usehs=F,debug=F,
		lb=c(1e-20,1e-10,1e-10,1e-8),ub=c(.5,1,.5,3),called=NULL,method='nm',...){
	if(!exists("rm.debug")){rm.debug<<-list(0);} 
	call<-sys.call(); tstart<-proc.time()[3];
	switch(model,
		g={np<-2;ex<-exg;hs<-hsg;keep<-c(1,2,5,6);},
		gm={np<-3;ex<-exgm;hs<-hsgm;keep<-c(1:3,5:7);},
		l={np<-3;ex<-exl;hs<-hsl;keep<-c(1,2,4,5,6,8);},
		lm={np<-4;ex<-exlm;hs<-hslm;keep<-1:8;},
		stop("The model you specified, ",model,", does not exist.
The only models supported in this version are: 
g, gm, l, and lm.")
	);
	if(!usegr) gr<-NULL; if(!usehs) hs<-NULL;
	rm.temp<<-list(cons=checkcons(cons,np,keep),cons2=rep(cons,le=4),x=x,y=y,
			keep=keep[1:np],old=old,debug=debug,called=called,np=np);
	par<-fixpar(par,np,keep);
	if(sum(is.na(par)>0)){
	stop("At least one of the starting parameters (par) you 
specified is an NA or NaN. Here are the parameters:\n",
paste(par,collapse=' '),
"\nPlease specify a different set of parameters or just omit them 
and use the defaults.");}
	lb<-fixpar(lb,np,keep);ub<-fixpar(ub,np,keep); rm.temp$lb<<-lb;rm.temp$ub<<-ub;
	fn<-function(par){objf(par,lb=lb,ub=ub,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],ex=ex);}
	gr<-function(par){grf(par,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],model=model);}
	change<-1; out<-list(par); totaliter<-0;
	#out.err<-list(estimate=rep.int(NA,np),maximum=NA,code=NA,message=NA);
        ctrl$parscale<-rep.int(1, length(par)); ctrl$ndeps<-rep.int(0.001, length(par));
	while(max(change)>0){
		out.old<-out;
		cat(".");
		options(show.error.messages=F);
		out<-try(.Internal(optim(out.old[[1]],fn,gr,"Nelder-Mead",ctrl,lb,ub)))
		options(show.error.messages=T);
		if(class(out)[1]!="try-error"){
			totaliter<-out[[3]][1]+totaliter;
			change<-abs(out.old[[1]]-out[[1]]);
		} else {
			p<-out.old[[1]]; lk<-NA;
			# Bad starting params, huh? We'll see about that.
			while(is.na(lk)|is.infinite(lk)){
			 p<-p/2; lbp<-p>lb;
			 p<-p*lbp+1.1*lb*(1-lbp);
			 # Damn starting params! You have bested me, I surrender.
			 if(max(p>1.1*lb)==0){cat('X\n');print(p);break();}
			 lk<-objf(p,lb,ub,rep(cons,le=4),
			          max(x),max(y),np=np,keep=keep[1:np],ex=ex);
			 cat('!');
			}
			if(!is.na(lk)){out<-list(p);change<-1;} else {break;}
		} 
	}
	cat(model," ");
	if(class(out)[1]!="try-error"){
		names(out)<-c('estimate','maximum','iterations','code','message');
		out$gradient<-gr(out$estimate);
		out$titer<-totaliter; 
		out$runtime<-proc.time()[3]-tstart;
	}
	if(debug){cat("Runtime:",proc.time()[3],"\n");}
	#out$totaliter<-totaliter;
	out;
}

