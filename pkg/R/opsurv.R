# This is the old version, which calls the interpreted versions of the functions and will be removed once the new compiled version is verified
`opsurv2` <-
function(model,x,y=NULL,par=rep(1e-6,4),cons=1,old=F,usegr=T,usehs=F,debug=F,
		lb=c(1e-20,1e-10,1e-10,1e-8),ub=c(.5,1,.5,3),called=NULL,method='nm',...){
	if(!exists("rm.debug")){rm.debug<<-list(0);} 
	call<-sys.call(); tstart<-proc.time()[3];
	# below needs to have the function names replaced with chars: "ofg","ofgm","ofl","oflm"
	switch(model,
		g={np<-2;ex<-exg;hs<-hsg;keep<-c(1,2,5,6);},
		gm={np<-3;ex<-exgm;hs<-hsgm;keep<-c(1:3,5:7);},
		l={np<-3;ex<-exl;hs<-hsl;keep<-c(1,2,4,5,6,8);},
		lm={np<-4;ex<-exlm;hs<-hslm;keep<-1:8;},
		stop("The model you specified, ",model,", does not exist.
The only models supported in this version are: 
g, gm, l, and lm.")
	);
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
	fn<-function(par){objf2(par,lb=lb,ub=ub,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],ex=ex);}
	gr<-function(par){grf2(par,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],model=model);}
	if(!usegr) gr<-NULL; if(!usehs) hs<-NULL;
	change<-1; out<-list(par); totaliter<-0;
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
			 if(max(p>1.1*lb)==0){cat('X\n');print(p);break();}
			 lk<-objf2(p,lb,ub,rep(cons,le=4),
			          max(x),max(y),np=np,keep=keep[1:np],ex=ex);
			 # Damn starting params! You have bested me, I surrender.
			 cat('!');
			}
			if(!is.na(lk)){out<-list(p);change<-1;} else {break;}
		} 
	}
	cat(model," ");
	if(class(out)[1]!="try-error"){
		names(out)<-c('estimate','maximum','iterations','code','message');
		if(usegr) { out$gradient<-gr(out$estimate); }
		out$call<-call;
		out$titer<-totaliter; 
		out$runtime<-proc.time()[3]-tstart;
	}
	if(debug){cat("Runtime:",proc.time()[3],"\n");}
	#out$totaliter<-totaliter;
	out;
}

`opsurv` <-
function(x,y=NULL,model='g',par=c(2.6e-6,.004,1e-10,0.1),cons=1,usegr=T,usehs=F,debug=F,
		lb=c(1e-20,1e-11,-4.940656e-324,-4.940656e-324),ub=c(.5,1,.5,3),called=NULL,cx=NULL,cy=NULL,
		mvers='',method='nm',debugcomment=NULL,...){
	callargs<-as.list(environment(),all.names=T); tstart<-proc.time()[3];
	if(!exists("rm.debug")){rm.debug<<-list(0);} 
	if(mean(cons)==1 & !is.null(y)){
		out<-list();
		options('warn'= -1);
		out0<-opsurv(x,model=model,par=par,cons=cons,usegr=usegr,usehs=usehs,debug=debug,
			     lb=lb,ub=ub,called=called,cx=cx,cy=cy,mvers=mvers,method=method,
			     debugcomment=debugcomment);
		out1<-opsurv(y,model=model,par=par,cons=cons,usegr=usegr,usehs=usehs,debug=debug,
			     lb=lb,ub=ub,called=called,cx=cx,cy=cy,mvers=mvers,method=method,
			     debugcomment=debugcomment);
		options('warn'=0);
		out$estimate<-c(out0$estimate,out1$estimate);
		out$maximum<-sum(out0$maximum,out1$maximum);
		out$iterations<-sum(out0$iterations,out1$iterations);
		out$code<-sum(out0$code,out1$code);
		out$message<-paste(out0$message,out1$message);
		out$gradient <- c(out0$gradient,out1$gradient);
		out$titer <- sum(out0$titer,out1$titer);
		out$runtime <- proc.time()[3] - tstart;
		out$group0<-out0; out$group1<-out1;
	} else {
		call<-sys.call();
		switch(model,
			g={np<-2;ex<-paste("ofg",mvers,sep='');hs<-hsg;keep<-c(1,2,5,6);},
			gm={np<-3;ex<-paste("ofgm",mvers,sep='');hs<-hsgm;keep<-c(1:3,5:7);},
			l={np<-3;ex<-paste("ofl",mvers,sep='');hs<-hsl;keep<-c(1,2,4,5,6,8);},
			lm={np<-4;ex<-paste("oflm",mvers,sep='');hs<-hslm;keep<-1:8;},
			stop("The model you specified, ",model,", does not exist.
The only models supported in this version are: 
g, gm, l, and lm.")
		);
		rm.temp<<-list(cons=checkcons(cons,np,keep),cons2=rep(cons,le=4),x=x,y=y,
				keep=keep[1:np],debug=debug,called=called,np=np);
		par<-fixpar(par,np,keep);
		if(sum(is.na(par)>0)){
		stop("At least one of the starting parameters (par) you 
specified is an NA or NaN. Here are the parameters:\n",
	paste(par,collapse=' '),
	"\nPlease specify a different set of parameters or just omit them 
and use the defaults.");}
		lb<-fixpar(lb,np,keep);ub<-fixpar(ub,np,keep); rm.temp$lb<<-lb;rm.temp$ub<<-ub;
		nx<-length(x);ny<-length(y);if(is.null(cx)){cx<-rep(1,nx);}
		if(is.null(cy)){if(is.null(y)){cy<-1;} else {cy<-rep(1,ny);}}
		fn<-function(par){objf(par,lb=lb,ub=ub,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],ex=ex,nx=nx,ny=ny,cx=cx,cy=cy);}
		gr<-function(par){grf(par,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],model=paste("g",model,mvers,sep=''),nx=nx,ny=ny,cx=cx,cy=cy);}
		if(!usegr) gr<-NULL; if(!usehs) hs<-NULL;
		change<-1; out<-list(par); totaliter<-0;
		ctrl$parscale<-rep.int(1, length(par)); ctrl$ndeps<-rep(0.001, length(par));
		# is there are reason not to jack up the iterations?
		ctrl$maxit<-5000;
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
				if(max(p>1.1*lb)==0){cat('X\n');print(p);break();}
				lk<-objf(p,lb=lb,ub=ub,cons=rep(cons,le=4),
					x=max(x),y=max(y),np=np,keep=keep[1:np],ex=ex,nx=nx,ny=ny,cx=cx,cy=cy);
				# Damn starting params! You have bested me, I surrender.
				cat('!');
				}
				if(!is.na(lk)){out<-list(p);change<-1;} else {break;}
			} 
		}
		cat(model," ");
		if(class(out)[1]!="try-error"){
			names(out)<-c('estimate','maximum','iterations','code','message');
			if(usegr) { out$gradient<-gr(out$estimate); }
			if(debug){out$call<-call; out$callargs<-callargs;}
			out$titer<-totaliter;
			out$runtime<-proc.time()[3]-tstart;
		}
		if(debug){cat("Runtime:",proc.time()[3],"\n");}
	}
	out;
}

