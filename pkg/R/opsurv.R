`opsurv` <-
function(x,y=NULL,model='g',par=c(2.6e-6,.004,1e-7,0.1),cons=1,usegr=F,usehs=F,debug=F,
	 lb=c(1e-14,1e-4,0,0),ub=c(.5,.5,.5,2),cx=NULL,cy=NULL,
 	 mvers='',method="Nelder-Mead",tlog=F){
	callargs<-as.list(environment(),all.names=T); tstart<-proc.time()[3];
	call<-sys.call();
	switch(model,
		w={np<-2;ex<-"ofw";hs<-NULL;keep<-c(1,2,5,6);ub=c(.5,10,0,0);},
		g={np<-2;ex<-paste("ofg",mvers,sep='');hs<-hsg;keep<-c(1,2,5,6);},
		gm={np<-3;ex<-paste("ofgm",mvers,sep='');hs<-hsgm;keep<-c(1:3,5:7);},
		l={np<-3;ex<-paste("ofl",mvers,sep='');hs<-hsl;keep<-c(1,2,4,5,6,8);},
		lm={np<-4;ex<-paste("oflm",mvers,sep='');hs<-hslm;keep<-1:8;},
		stop("The model you specified, ",model,", does not exist.
The only models supported in this version are: 
g, gm, l, and lm.")
		);
	if(mean(cons)==1 & !is.null(y)){
		out<-list();
		options('warn'= -1);
		out0<-opsurv(x,model=model,par=par,cons=cons,usegr=usegr,usehs=usehs,debug=debug,lb=lb,ub=ub,cx=cx,cy=cy,mvers=mvers,method=method,tlog=tlog);
			     if(length(par)==2*np){par1<-par[(np+1):(2*np)];
			     } else {if(length(par)==8){par1<-par[5:8];
			     } else par1<-par;}
		out1<-opsurv(y,model=model,par=par1,cons=cons,usegr=usegr,usehs=usehs,debug=debug,lb=lb,ub=ub,cx=cx,cy=cy,mvers=mvers,method=method,tlog=tlog);
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
		rm.temp<<-list(cons=checkcons(cons,np,keep),cons2=rep(cons,le=4),x=x,y=y,
				keep=keep[1:np],debug=debug,np=np);
		par<-fixpar(par,np,keep);
		if(sum(is.na(par)>0)){
		stop("At least one of the starting parameters (par) you 
specified is an NA or NaN. Here are the parameters:\n",
	paste(par,collapse=' '),
	"\nPlease specify a different set of parameters or just omit them 
and use the defaults.");}
# 		browser();
		lb<-lb[keep[1:np]]; ub<-ub[keep[1:np]] 
		lb<-fixpar(lb,np,keep);ub<-fixpar(ub,np,keep); 
		if(tlog){par<-log(par);}
		rm.temp$lb<<-lb;rm.temp$ub<<-ub;
		nx<-length(x);ny<-length(y);if(is.null(cx)){cx<-rep(1,nx);}
		if(is.null(cy)){if(is.null(y)){cy<-1;} else {cy<-rep(1,ny);}}
		fn<-function(par){
			objf(par,lb=lb,ub=ub,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],ex=ex,nx=nx,ny=ny,cx=cx,cy=cy,tlog=tlog);
		}
		gr<-function(par){grf(par,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],model=paste("g",model,mvers,sep=''),nx=nx,ny=ny,cx=cx,cy=cy,tlog=tlog);}
		if(!usegr) gr<-NULL; if(!usehs) hs<-NULL;
		change<-1; out<-list(par); totaliter<-0;
		ctrl$parscale<-rep.int(1, length(par)); ctrl$ndeps<-rep(0.001, length(par));
		# is there are reason not to jack up the iterations?
		ctrl$maxit<-5000;
		while(max(change)>0){
			out.old<-out;
			cat(".");
			options(show.error.messages=F);
			out<-try(.Internal(optim(out.old[[1]],fn,gr,method,ctrl,lb,ub)))
			options(show.error.messages=T);
			if(class(out)[1]!="try-error"){
				totaliter<-out[[3]][1]+totaliter;
				change<-abs(out.old[[1]]-out[[1]]);
			} else {
				p<-out.old[[1]]; lk<-NA;
				# Bad starting params, huh? We'll see about that.
				while(is.na(lk)|is.infinite(lk)){
				if(tlog){p[is.infinite(p)]<- 0; p<-(p-5)/2;
					if(max(p-log(lb))<1e-10){
						cat('X\n');print(p);print(exp(lb));
						print(lb);break();
					}
				}else{
					p<-p/2; lbp<-p>lb;
					p<-p*lbp+1.1*lb*(1-lbp);
					if(max(p>1.1*lb)==0){cat('X\n');print(p);break();}
				}
				if(is.null(y)){maxy<-NULL;} else {maxy<-max(y);}
				lk<-objf(p,lb=lb,ub=ub,cons=rep(cons,le=4),
					 x=max(x),y=maxy,np=np,keep=keep[1:np],ex=ex,nx=1,ny=1,cx=1,cy=1,
					 tlog=tlog);
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
			if(tlog){out$estimate<-exp(out$estimate);}
		}
		if(debug){cat("Runtime:",proc.time()[3],"\n");}
	}
	out;
}
