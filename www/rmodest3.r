# expressions for the density functions, evaluated by likelihood functions, gradients, and hessians
exg<-expression(log(a*exp(b*x-a*(exp(b*x)-1)/b)));
exgm<-expression(log((c+a*exp(b*x))*exp(-a*(exp(b*x)-1)/b-c*x)));
exl<-expression(log(a*exp(b*x)*(1+s*a*(exp(b*x)-1)/b)^(-(s+1)/s)));
exlm<-expression(log((c+a*exp(b*x)/(1+s*a*(exp(b*x)-1)/b))*(exp(-c*x)*(1+s*a*(exp(b*x)-1)/b)^(-1/s))));

# a nested list of derivatives, for constructing gradients
dex<-list();
dex$g<-list(a=D(exg,'a'),b=D(exg,'b'));
dex$gm<-list(a=D(exgm,'a'),b=D(exgm,'b'),c=D(exgm,'c'));
dex$l<-list(a=D(exl,'a'),b=D(exl,'b'),s=D(exl,'s'));
dex$lm<-list(a=D(exlm,'a'),b=D(exlm,'b'),c=D(exlm,'c'),s=D(exlm,'s'));

# generic gradient function
grf<-function(par,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,model='g'){
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

# generic likelihood function
objf<-function(par,lb=0,ub=1,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,ex=exg){
	if(max(par<=lb)!=0|max(par>=ub)!=0){return(NA);}
	par1<-par2<-rep.int(NA,4); par1[keep]<-par2[keep]<-par[1:np];
	par2[keep][cons[keep]>0]<-par[(np+1):length(par)];
	env1<-list(a=par1[1],b=par1[2],c=par1[3],s=par1[4],x=x); 
	env2<-list(a=par2[1],b=par2[2],c=par2[3],s=par2[4],x=y);
	out<-sum(eval(ex,envir=env1))+sum(eval(ex,envir=env2));
	if(!is.na(out) & is.finite(out) & out < -50){return(out);}else{return(NA);}
}

# control parameters for the optim internal library
ctrl<-list(trace=0,fnscale=-1,maxit = 100000, 
           abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
           beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
           factr = 1e+07, pgtol = 0, tmax = 10, temp = 10);

# main optimization function
opsurv<-function(model,x,y=NULL,par=rep(1e-6,4),cons=1,old=F,usegr=T,usehs=F,debug=F,
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

# checks constraints to make sure they have legal values
checkcons<-function(cons,np,keep){
	mykeep<-keep[1:np];
	cons<-as.logical(cons);
	if(length(cons)==1){return(rep.int(T,np));}
	if(length(cons)==4){return(cons[mykeep]);}
	if(length(cons)==np){return(cons);}
	stop("If a constraint (cons) is specified, is should be
  a vector of logical values equal in length either
  to 4 or to the number of parameters in your model.");
}

# checks and formats starting parameters
fixpar<-function(par,np,keep){
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

# hessians

hsg<-function(par,cons=NULL,x=NULL,y=NULL){
	if(is.null(cons)) cons<-rm.temp$cons;
	if(is.null(x)) { x<-rm.temp$x; y<-rm.temp$y; }
	par2<-par[1:2]; par2[cons>0]<-par[3:length(par)];
	env1<-list(a=par[1],b=par[2],x=rm.temp$x);
	env2<-list(a=par2[1],b=par2[2],x=rm.temp$y);
	out1a<-c(sum(eval(D(D(exg,'a'),'a'),envir=env1)),
		sum(eval(D(D(exg,'b'),'a'),envir=env1)));
	out1b<-c(sum(eval(D(D(exg,'a'),'b'),envir=env1)),
		sum(eval(D(D(exg,'b'),'b'),envir=env1)));
	out1<-rbind(out1a,out1b); 
	out1<-rbind(cbind(out1,out1),cbind(out1,out1));
	out2a<-c(sum(eval(D(D(exg,'a'),'a'),envir=env2)),
		sum(eval(D(D(exg,'b'),'a'),envir=env2)));
	out2b<-c(sum(eval(D(D(exg,'a'),'b'),envir=env2)),
		sum(eval(D(D(exg,'b'),'b'),envir=env2)));
	out2<-rbind(out2a,out2b); 
	out2<-rbind(cbind(out2,out2),cbind(out2,out2));
	id1<-c(1,1,1-cons); id2<-c(1-cons,1,1);
	id1<-id1%o%id1; id2<-id2%o%id2;
	out<-out1*id1+out2*id2;
	if(rm.temp$debug){print("hessian:");print(out);}
	out[c(1:2,3:4*cons),c(1:2,3:4*cons)];
}

hsgm<-function(par,cons=NULL,x=NULL,y=NULL){
	if(is.null(cons)) cons<-rm.temp$cons;
	if(is.null(x)) { x<-rm.temp$x; y<-rm.temp$y; }
	par2<-par[1:3]; par2[cons>0]<-par[4:length(par)];
	env1<-list(a=par[1],b=par[2],c=par[3],x=rm.temp$x);
	env2<-list(a=par2[1],b=par2[2],c=par2[3],x=rm.temp$y);
	out1a<-c(sum(eval(D(D(exgm,'a'),'a'),envir=env1)),
		sum(eval(D(D(exgm,'b'),'a'),envir=env1)),
		sum(eval(D(D(exgm,'c'),'a'),envir=env1)));
	out1b<-c(sum(eval(D(D(exgm,'a'),'b'),envir=env1)),
		sum(eval(D(D(exgm,'b'),'b'),envir=env1)),
		sum(eval(D(D(exgm,'c'),'b'),envir=env1)));
	out1c<-c(sum(eval(D(D(exgm,'a'),'c'),envir=env1)),
		sum(eval(D(D(exgm,'b'),'c'),envir=env1)),
		sum(eval(D(D(exgm,'c'),'c'),envir=env1)));
	out1<-rbind(out1a,out1b,out1c); 
	out1<-rbind(cbind(out1,out1),cbind(out1,out1));
	out2a<-c(sum(eval(D(D(exgm,'a'),'a'),envir=env2)),
		sum(eval(D(D(exgm,'b'),'a'),envir=env2)),
		sum(eval(D(D(exgm,'c'),'a'),envir=env2)));
	out2b<-c(sum(eval(D(D(exgm,'a'),'b'),envir=env2)),
		sum(eval(D(D(exgm,'b'),'b'),envir=env2)),
		sum(eval(D(D(exgm,'c'),'b'),envir=env2)));
	out2c<-c(sum(eval(D(D(exgm,'a'),'c'),envir=env2)),
		sum(eval(D(D(exgm,'b'),'c'),envir=env2)),
		sum(eval(D(D(exgm,'c'),'c'),envir=env2)));
	out2<-rbind(out2a,out2b,out2c); 
	out2<-rbind(cbind(out2,out2),cbind(out2,out2));
	id1<-c(1,1,1,1-cons); id2<-c(1-cons,1,1,1);
	id1<-id1%o%id1; id2<-id2%o%id2;
	out<-out1*id1+out2*id2;
	if(rm.temp$debug){print("hessian:");print(out);}
	out[c(1:3,4:6*cons),c(1:3,4:6*cons)];
}

hsl<-function(par,cons=NULL,x=NULL,y=NULL){
	if(is.null(cons)) cons<-rm.temp$cons;
	if(is.null(x)) { x<-rm.temp$x; y<-rm.temp$y; }
	par2<-par[1:3]; par2[cons>0]<-par[4:length(par)];
	env1<-list(a=par[1],b=par[2],s=par[3],x=rm.temp$x);
	env2<-list(a=par2[1],b=par2[2],s=par2[3],x=rm.temp$y);
	out1a<-c(sum(eval(D(D(exl,'a'),'a'),envir=env1)),
		sum(eval(D(D(exl,'b'),'a'),envir=env1)),
		sum(eval(D(D(exl,'s'),'a'),envir=env1)));
	out1b<-c(sum(eval(D(D(exl,'a'),'b'),envir=env1)),
		sum(eval(D(D(exl,'b'),'b'),envir=env1)),
		sum(eval(D(D(exl,'s'),'b'),envir=env1)));
	out1s<-c(sum(eval(D(D(exl,'a'),'s'),envir=env1)),
		sum(eval(D(D(exl,'b'),'s'),envir=env1)),
		sum(eval(D(D(exl,'s'),'s'),envir=env1)));
	out1<-rbind(out1a,out1b,out1s); 
	out1<-rbind(cbind(out1,out1),cbind(out1,out1));
	out2a<-c(sum(eval(D(D(exl,'a'),'a'),envir=env2)),
		sum(eval(D(D(exl,'b'),'a'),envir=env2)),
		sum(eval(D(D(exl,'s'),'a'),envir=env2)));
	out2b<-c(sum(eval(D(D(exl,'a'),'b'),envir=env2)),
		sum(eval(D(D(exl,'b'),'b'),envir=env2)),
		sum(eval(D(D(exl,'s'),'b'),envir=env2)));
	out2s<-c(sum(eval(D(D(exl,'a'),'s'),envir=env2)),
		sum(eval(D(D(exl,'b'),'s'),envir=env2)),
		sum(eval(D(D(exl,'s'),'s'),envir=env2)));
	out2<-rbind(out2a,out2b,out2s); 
	out2<-rbind(cbind(out2,out2),cbind(out2,out2));
	id1<-c(1,1,1,1-cons); id2<-c(1-cons,1,1,1);
	id1<-id1%o%id1; id2<-id2%o%id2;
	out<-out1*id1+out2*id2;
	if(rm.temp$debug){print("hessian:");print(out);}
	out[c(1:3,4:6*cons),c(1:3,4:6*cons)];
}

hslm<-function(par,cons=NULL,x=NULL,y=NULL){
	if(is.null(cons)) cons<-rm.temp$cons;
	if(is.null(x)) { x<-rm.temp$x; y<-rm.temp$y; }
	par2<-par[1:4]; par2[cons>0]<-par[5:length(par)];
	env1<-list(a=par[1],b=par[2],c=par[3],s=par[4],x=rm.temp$x);
	env2<-list(a=par2[1],b=par2[2],c=par2[3],s=par2[4],x=rm.temp$y);
	out1a<-c(sum(eval(D(D(exlm,'a'),'a'),envir=env1)),
		sum(eval(D(D(exlm,'b'),'a'),envir=env1)),
		sum(eval(D(D(exlm,'c'),'a'),envir=env1)),
		sum(eval(D(D(exlm,'s'),'a'),envir=env1)));
	out1b<-c(sum(eval(D(D(exlm,'a'),'b'),envir=env1)),
		sum(eval(D(D(exlm,'b'),'b'),envir=env1)),
		sum(eval(D(D(exlm,'c'),'b'),envir=env1)),
		sum(eval(D(D(exlm,'s'),'b'),envir=env1)));
	out1c<-c(sum(eval(D(D(exlm,'a'),'c'),envir=env1)),
		sum(eval(D(D(exlm,'b'),'c'),envir=env1)),
		sum(eval(D(D(exlm,'c'),'c'),envir=env1)),
		sum(eval(D(D(exlm,'s'),'c'),envir=env1)));
	out1s<-c(sum(eval(D(D(exlm,'a'),'s'),envir=env1)),
		sum(eval(D(D(exlm,'b'),'s'),envir=env1)),
		sum(eval(D(D(exlm,'c'),'s'),envir=env1)),
		sum(eval(D(D(exlm,'s'),'s'),envir=env1)));
	out1<-rbind(out1a,out1b,out1c,out1s); 
	out1<-rbind(cbind(out1,out1),cbind(out1,out1));
	out2a<-c(sum(eval(D(D(exlm,'a'),'a'),envir=env2)),
		sum(eval(D(D(exlm,'b'),'a'),envir=env2)),
		sum(eval(D(D(exlm,'c'),'a'),envir=env2)),
		sum(eval(D(D(exlm,'s'),'a'),envir=env2)));
	out2b<-c(sum(eval(D(D(exlm,'a'),'b'),envir=env2)),
		sum(eval(D(D(exlm,'b'),'b'),envir=env2)),
		sum(eval(D(D(exlm,'c'),'b'),envir=env2)),
		sum(eval(D(D(exlm,'s'),'b'),envir=env2)));
	out2c<-c(sum(eval(D(D(exlm,'a'),'c'),envir=env2)),
		sum(eval(D(D(exlm,'b'),'c'),envir=env2)),
		sum(eval(D(D(exlm,'c'),'c'),envir=env2)),
		sum(eval(D(D(exlm,'s'),'c'),envir=env2)));
	out2s<-c(sum(eval(D(D(exlm,'a'),'s'),envir=env2)),
		sum(eval(D(D(exlm,'b'),'s'),envir=env2)),
		sum(eval(D(D(exlm,'c'),'s'),envir=env2)),
		sum(eval(D(D(exlm,'s'),'s'),envir=env2)));
	out2<-rbind(out2a,out2b,out2c,out2s); 
	out2<-rbind(cbind(out2,out2),cbind(out2,out2));
	id1<-c(1,1,1,1,1-cons); id2<-c(1-cons,1,1,1,1);
	id1<-id1%o%id1; id2<-id2%o%id2;
	out<-out1*id1+out2*id2;
	if(rm.temp$debug){print("hessian:");print(out);}
	out[c(1:4,5:8*cons),c(1:4,5:8*cons)];
}