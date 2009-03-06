`scramble` <-
function(x,y,n,srt=T){
	lx<-length(x);ly<-length(y);
	xy <- c(x,y); lxy<-length(xy);
	outx<-c(); outy<-c();
	for(i in 1:n){
		ixy<-sample(xy,lxy)
		ix<-ixy[1:lx]; iy<-ixy[(lx+1):(lx+ly)];
		if(srt){ix<-sort(ix); iy<-sort(iy);}
		outx<-cbind(outx,ix); outy<-cbind(outy,iy);
	}
	return(list(x=outx,y=outy));
}

`permstat` <-
function(x,y,n,fn,srt=F,statname='statistic',...){
	fnout<-fn(x,y,...);
	perms<-scramble(x,y,n,srt);
	sout<-list(); length(sout)<-length(statname);
	names(sout)<-statname;
	for(i in 1:n){
		ifnout<-fn(perms$x[,i],perms$y[,i],...);
		for(j in statname){ sout[[j]]<-c(sout[[j]],ifnout[[j]]);}
	}
	for(j in statname){
		jecdf<-ecdf(sout[[j]]);
		cat(j,':',fnout[[j]],'permuted p:',1-jecdf(fnout[[j]]),'\n');
		cat('empirical quantiles:',quantile(sout[[j]],c(.05,.1,.5,.9,.95)),'\n');
	}
	invisible(list(f=fnout,perms=sout));
}

`twrapper` <-
function(data,foo,...){
	foox<-data[,1][foo[data[,2]==0]]; fooy=data[,1][foo[data[,2]==1]];
	out<-t.test(x=foox,y=fooy,alternative='two.sided',mu=0);
	out$statistic;
}

`startpars`<-
function(x,y,label,perms=10){
	d<-int2tab(x,y); lxy<-dim(d)[1]; pr<-matrix(nrow=lxy,ncol=perms);
	for(i in 1:perms){pr[,i]<-sample(d[,2],lxy);}
	fits<-list(lm=list(u=list(cons=c(1L,1L,1L,1L),out=c()), a=list(cons=c(0L,1L,1L,1L),out=c()),
			   b=list(cons=c(1L,0L,1L,1L),out=c()), c=list(cons=c(1L,1L,0L,1L),out=c()),
			   s=list(cons=c(1L,1L,1L,0L),out=c())),
		   l=list(u=list(cons=c(1L,1L,1L,1L),out=c()), a=list(cons=c(0L,1L,1L,1L),out=c()),
			  b=list(cons=c(1L,0L,1L,1L),out=c()), s=list(cons=c(1L,1L,1L,0L),out=c())),
		   gm=list(u=list(cons=c(1L,1L,1L,1L),out=c()), a=list(cons=c(0L,1L,1L,1L),out=c()),
			   b=list(cons=c(1L,0L,1L,1L),out=c()), c=list(cons=c(1L,1L,0L,1L),out=c())),
		   g=list(u=list(cons=c(1L,1L,1L,1L),out=c()), a=list(cons=c(0L,1L,1L,1L),out=c()),
			  b=list(cons=c(1L,0L,1L,1L),out=c())));	

	cat('\nCalculating unconstrained fits from default starting parameters.\n');
	for(i in 1:perms){
		for(j in 1:length(fits)){
			jout<-opsurv(d[pr[,i]==0,1],d[pr[,i]==1,1],names(fits)[j]);
			fits[[j]]$u$out<-rbind(fits[[j]]$u$out,
					       c(jout$maximum,jout$titer,jout$runtime,jout$estimate));
		}
	}

	cat('\nCalculating unconstrained fits from candidate optimal parameters.\n');
	for(j in 1:length(fits)){
		#browser();
		jcols<-4:dim(fits[[j]]$u$out)[2];
		jin<-fits[[j]]$u$out[-1,jcols];
		jmaxs<-numeric(4); jpars<-list();
		jmaxs[1]<-fits[[j]]$u$out[1,1];

		jpars[[1]]<-fits[[j]]$u$out[1,jcols];
		jpars[[2]]<-apply(jin,2,mean);
		jpars[[3]]<-apply(jin,2,median);
		jpars[[4]]<-apply(jin,2,gmean);

		for(l in 2:4){
			jmaxs[l]<-opsurv(x,y,names(fits)[j],par=jpars[[l]])$maximum;
		}
		fits[[j]]$u$label<-paste(label,'_',names(fits)[j],'u',sep='');
		jmaxmax<-which.max(jmaxs); if(jmaxmax>1){
			cat('\nBetter pars found for',fits[[j]]$u$label,':',jmaxmax,'\n');
		}
		fits[[j]]$u$par<-jpars[[jmaxmax]];
		#cat('\nmaxima:',jmaxs,'\n');
		write.table(cbind(jmaxs,t(as.data.frame(jpars))),
			    file=paste(fits[[j]]$u$label,'_opars.txt',sep=''),
			    col.names=F,row.names=F,sep='\t');
	}
	
	cat('\nCalculating constrained fits from unconstrained optimal parameters.\n');
	for(i in 1:perms){
		for(j in 1:length(fits)){
			for(k in 2:length(fits[[j]])){
				jkout<-opsurv(d[pr[,i]==0,1],d[pr[,i]==1,1],
					      names(fits)[j],
					      cons=fits[[j]][[k]]$cons,
					      par=fits[[j]]$u$par);
				fits[[j]][[k]]$out<-rbind(fits[[j]][[k]]$out,
							  c(jkout$maximum,jkout$titer,jkout$runtime,	
							    jkout$estimate));
			}
		}
	}

	cat('\nCalculating constrained from candidate optimal parameters.\n');
	for(j in 1:length(fits)){
		for(k in 2:length(fits[[j]])){
			jkcols<-4:dim(fits[[j]][[k]]$out)[2];
			jkin<-fits[[j]][[k]]$out[-1,jkcols];
			jkmaxs<-numeric(4); jkpars<-list();
			jkmaxs[1]<-fits[[j]][[k]]$out[1,1];
			jkpars[[1]]<-fits[[j]]$u$out[1,jkcols];
			jkpars[[2]]<-apply(jkin,2,mean);
			jkpars[[3]]<-apply(jkin,2,median);
			jkpars[[4]]<-apply(jkin,2,gmean);
	
			for(l in 2:4){
				jkmaxs[l]<-opsurv(x,y,names(fits)[j],cons=fits[[j]][[k]]$cons,
				                  par=jkpars[[l]])$maximum;
			}
			fits[[j]][[k]]$label<-paste(label,'_',names(fits)[j],names(fits[[j]])[k],sep='');
			jkmaxmax<-which.max(jkmaxs); if(jkmaxmax>1){
				cat('\nBetter pars found for',fits[[j]][[k]]$label,':',jkmaxmax,'\n');
			}
			fits[[j]][[k]]$par<-jkpars[[jkmaxmax]];
			write.table(cbind(jkmaxs,t(as.data.frame(jkpars))),
				    file=paste(fits[[j]][[k]]$label,'_opars.txt', sep=''), 
				    col.names=F,row.names=F,sep='\t');
		}
	}
	outname<-paste(label,'pars',sep='');
	assign(outname,fits,envir=.GlobalEnv);
	save(list=outname,file=paste(outname,'.rdata',sep=''));
	invisible(fits);
}

`ressurv` <-
function(x,y,label,resam=1000,pars=NULL,...){
	lx<-length(x);ly<-length(y);
	rsx<-matrix(nrow=lx,ncol=resam);
	rsy<-matrix(nrow=ly,ncol=resam);
	for(i in 1:resam){
		rsx[,i]<-sample(x,lx,replace=T);
		rsy[,i]<-sample(y,ly,replace=T);
	}
	rsx<-cbind(x,rsx); rsy<-cbind(y,rsy);
	if(is.character(pars)){parslabel<-pars;pars<-NULL;}
	else{parslabel<-paste(label,'pars',sep='');}
	# need to write some kind of validation function for fits
	if(is.null(pars)|length(pars)!=4){
		if(exists(parslabel)){pars<-get(parslabel);}
		else{pars<-startpars(x,y,label,perms=10);}
	}
	for(i in 1:(resam+1)){
		for(j in names(pars)){
			for(k in names(pars[[j]])){
				out<-opsurv(rsx[,i],rsy[,i],model=j,cons=pars[[j]][[k]]$cons,
					    par=pars[[j]][[k]]$par);
				outtab<-c(out$maximum,out$titer,out$runtime,out$estimate,
					  out$group0$estimate,out$group1$estimate);
				write.table(t(outtab),paste(label,'_resam_',j,k,'.txt',sep=''),sep='\t',
					    append=T,col.names=F,row.names=F);
			}
		}
	}
}

`permsurv` <-
function(x,y,label,perms=5000,pars=NULL){
	d<-int2tab(x,y); lxy<-dim(d)[1]; pr<-matrix(nrow=lxy,ncol=perms);
	for(i in 1:perms){pr[,i]<-sample(d[,2],lxy);}
	pr<-cbind(d[,2],pr);
	if(is.character(pars)){parslabel<-pars;pars<-NULL;}
	else{parslabel<-paste(label,'pars',sep='');cat('\nLooking for',parslabel,'\n');}
	# need to write some kind of validation function for fits
	if(is.null(pars)|length(pars)!=4){
		if(exists(parslabel)){cat('\nFound',parslabel,'\n');pars<-get(parslabel);}
		else{pars<-startpars(x,y,label,perms=10);}
	}
	cat('\nCommencing permutations.\n');
	for(i in 1:(perms+1)){
		for(j in 1:length(pars)){
			for(k in 1:length(pars[[j]])){
				jkts<-as.numeric(Sys.time());
				jklabel<-paste(pars[[j]][[k]]$label,jkts,sep='');
				jkfit<-opsurv(d[pr[,i]==0,1],d[pr[,i]==1,1],names(pars)[j],
						cons=pars[[j]][[k]]$cons, par=pars[[j]][[k]]$par);
				jktab<-c(jkfit$maximum,jkfit$titer,jkfit$runtime,jkfit$estimate,
						jkfit$group0$maximum,jkfit$group1$maximum);
				assign(jklabel,jkfit);
				dump(jklabel, paste(pars[[j]][[k]]$label,'.R',sep=''),
				     append=T,control='S_compatible');
				#print('');print(jktab);
				write.table(t(jktab),paste(pars[[j]][[k]]$label,'.txt',sep=''),sep='\t',
					    append=T,col.names=F,row.names=F);
			}
		}
	}
	save(list=ls(),file=paste(label,'debug.rdata',sep=''));
}


`hazwrapper` <-
function(data,foo,outlog='hazlog.R',debuglog='hazdebug.R',debugvars=c('estimate','maximum','titers','runtime'),...){
	ts<-as.numeric(Sys.time()); data[,1]<-data[,1][foo]; 
	outname<-paste('out',ts,sep='.');
	assign(outname,opsurv('lm',data[data[,2]==0,1],data[data[,2]==1,1],...));
	if(!is.null(outlog)){
		dump(outname,outlog,append=T,control='S_compatible')
	}
	get(outname)$maximum;
}

`gmean`<-
function(x){exp(sum(log(x))/length(x));}

`parnames`<-
function(i,j){
	model<-constraint<-parameters<-c();
	switch(i,
		{model<-'LM'; switch(j,
				{constraint<-'u'; parameters<-c('a','b','c','s','a','b','c','s');},
				{constraint<-'a'; parameters<-c('a','b','c','s','b','c','s');},
				{constraint<-'b'; parameters<-c('a','b','c','s','a','c','s');},
				{constraint<-'c'; parameters<-c('a','b','c','s','a','b','s');},
				{constraint<-'s'; parameters<-c('a','b','c','s','a','b','c');}
			);
		},
		{model<-'L'; switch(j,
				{constraint<-'u'; parameters<-c('a','b','s','a','b','s');},
				{constraint<-'a'; parameters<-c('a','b','s','b','s');},
				{constraint<-'b'; parameters<-c('a','b','s','a','s');},
				{constraint<-'s'; parameters<-c('a','b','s','a','b');}
			);
		},
		{model<-'GM'; switch(j,
				{constraint<-'u'; parameters<-c('a','b','c','a','b','c');},
				{constraint<-'a'; parameters<-c('a','b','c','b','c');},
				{constraint<-'b'; parameters<-c('a','b','c','a','c');},
				{constraint<-'c'; parameters<-c('a','b','c','a','b');}
			);
		},
		{model<-'G'; switch(j,
				{constraint<-'u'; parameters<-c('a','b','a','b');},
				{constraint<-'a'; parameters<-c('a','b','b');},
				{constraint<-'b'; parameters<-c('a','b','a');}
			);
		}
	)
	list(model=model,constraint=constraint,parameters=parameters);
}

`spyhaz`<-
function(label,show=c('models','lm','l','gm','g','pars')){
	fit<-input<-list(lm=list(u=list(),a=list(),b=list(),c=list(),s=list()), 
			 l=list(u=list(),a=list(),b=list(),s=list()),
			 gm=list(u=list(),a=list(),b=list(),c=list()),
			 g=list(u=list(),a=list(),b=list()));
	models<-list(g.l=c(),g.gm=c(),l.lm=c(),gm.lm=c());
	class(fit)<-'spyfit';class(input)<-'spyinput';class(models)<-'spymodels';
	for(i in 1:length(input) ){
		for(j in 1:length(input[[i]]) ){
			input[[i]][[j]]<-read.table(paste(label,'_',names(input)[i],names(input[[i]])[j],
							  '.txt',sep=''));
			ijcifile<-paste(label,'_resam_',names(input)[i],names(input[[i]])[j],
					'.txt',sep='');
			if(file.exists(ijcifile)){
				ijci<-read.table(ijcifile);
				#attr(input[[i]][[j]],'ci')<-apply(ijci,2,quantile,c(.05,.95));
				attr(input[[i]][[j]],'ci')<-ijci;
			} else ijci<-NULL;
			if(j>1){ 
				fit[[i]][[j]]<-2*(input[[i]]$u[,1] - input[[i]][[j]][,1]);
				if(!is.null(ijci) & exists('iuci')){
# 					attr(fit[[i]][[j]],'ci')<-quantile(2*(iuci[,1]-ijci[,1]),
# 									   c(.05,.95));
					attr(fit[[i]][[j]],'ci')<-2*(iuci[,1]-ijci[,1]);
				}
			} else { iuci<-ijci; }
		}
	}
	models$g.l<-2*(input$l$u[,1] - input$g$u[,1]);
	attr(models$g.l,'ci')<-2*(attr(input$l$u,'ci')[,1]-attr(input$g$u,'ci')[,1]);
	models$g.gm<-2*(input$gm$u[,1] - input$g$u[,1]);
	attr(models$g.gm,'ci')<-2*(attr(input$gm$u,'ci')[,1]-attr(input$g$u,'ci')[,1]);
	models$l.lm<-2*(input$lm$u[,1] - input$l$u[,1]);
	attr(models$l.lm,'ci')<-2*(attr(input$lm$u,'ci')[,1]-attr(input$l$u,'ci')[,1]);
	models$gm.lm<-2*(input$lm$u[,1] - input$gm$u[,1]);
	attr(models$gm.lm,'ci')<-2*(attr(input$lm$u,'ci')[,1]-attr(input$gm$u,'ci')[,1]);
	models$g.lm<-2*(input$lm$u[,1] - input$g$u[,1]);
	attr(models$g.lm,'ci')<-2*(attr(input$lm$u,'ci')[,1]-attr(input$g$u,'ci')[,1]);

	assign(paste(label,'fit',sep=''),fit,envir=.GlobalEnv);
	assign(paste(label,'input',sep=''),input,envir=.GlobalEnv);
	assign(paste(label,'models',sep=''),models,envir=.GlobalEnv);
}

`print.spyfit`<-
function(x,type=c('plot','table'),...){
	if(match('plot',type,nomatch=0)){
		layout(matrix(1:12,4,3,byrow=T));
		for(i in 4:1){
				switch(i,{model='Logistic-Makeham';pars=1:4;},
					 {model='Logistic';pars=c(1,2,4);},
					 {model='Gompertz-Makeham';pars=1:3;},
					 {model='Gompertz';pars=1:2;},);
				for(j in 2:length(x[[i]])){
					param<-switch(pars[j-1],'a','b','c','s');
					ijecdf<-ecdf(x[[i]][[j]]);
					ijp<-1-ijecdf(x[[i]][[j]][1]);
					ijxlab<-paste('parameter',param);
					ijsub<-paste('p =',ijp);
					hist(x[[i]][[j]],breaks=100,xlab=ijxlab,ylab='',
					     freq=F,main=model,sub=ijsub,...);
					abline(v=x[[i]][[j]][1],col='red');
					abline(v=attr(x[[i]][[j]],'ci')[1],col='green',lty=2);
					hist(attr(x[[i]][[j]],'ci'),breaks=100,freq=F,border='pink',
					     lty=2,add=T);
				}
			}
		layout(matrix(1,1,1));
	}
	if(match('table',type,nomatch=0)){
		rnames<-out<-c();
		for(i in 1:length(x)){
			for(j in 2:length(x[[i]])){
				out<-rbind(out,c(quantile(x[[i]][[j]],.05),median(x[[i]][[j]]),
						 x[[i]][[j]][1],mean(x[[i]][[j]]), 
						 quantile(x[[i]][[j]],.95),
						 1-ecdf(x[[i]][[j]])(x[[i]][[j]][1]),
						 pchisq(x[[i]][[j]][i],1,lower=F)));
				out<-rbind(out,c(quantile(attr(x[[i]][[j]],'ci'),.05),
						 median(attr(x[[i]][[j]],'ci')),
						 attr(x[[i]][[j]],'ci')[1],mean(attr(x[[i]][[j]],'ci')),
						 quantile(attr(x[[i]][[j]],'ci'),.95),NA,NA));
				ijlabel<-paste(toupper(names(x))[i],',',names(x[[i]])[j],'constrained');
				rnames<-c(rnames,paste('Permuted',ijlabel),paste('Resampled',ijlabel));
			}
		}
		colnames(out)<-c('5%','median','actual','mean','95%','emp p','chisq p');
		rownames(out)<-rnames;
		oldwidth<-options('width')[[1]];options('width'=130);
		print(out);options('width'=oldwidth);invisible(out);
	}
}

`print.spyinput`<-
function(x,type=c('plot','table'),...){
	if(match('plot',type,nomatch=0)){
		layout(matrix(1:12,4,3,byrow=T));
		for(i in 4:1){
				switch(i,{model='Logistic-Makeham';pars=1:4;},
					 {model='Gompertz-Makeham';pars=1:3;},
					 {model='Logistic';pars=c(1,2,4);},
					 {model='Gompertz';pars=1:2;},);
				for(j in 1:length(pars)){
					param<-switch(pars[j],'a','b','c','s');
					ijxlab<-paste('parameter',param);
					ijx0<-x[[i]]$u[,3+j];
					ijx1<-x[[i]]$u[,3+j+length(pars)];
					hist(ijx0[ijx0>0],breaks=100,xlab=ijxlab,freq=F,
					     ylab='',main=model,...);
					hist(ijx1[ijx1>0],breaks=100,lty=2,add=T,freq=F);
					hist(attr(x[[i]]$u,'ci')[,3+j],breaks=100,border='red',freq=F,
					     lty=2,add=T);
					hist(attr(x[[i]]$u,'ci')[,3+j+length(pars)],breaks=100,
					     border='red',freq=F,lty=3,add=T);
				}
			}
		layout(matrix(1,1,1));
	}
	if(match('table',type,nomatch=0)){
		rnames<-out<-c();
		for(i in 1:length(x)){
			for(j in 1:length(x[[i]])){
				ijnames<-parnames(i,j); 
				ijnames$parameters<-c('MLE','iter','rt',ijnames$parameters);
				for(k in 1:length(x[[i]][[j]])){
					out<-rbind(out,c(quantile(x[[i]][[j]][[k]],.05),
						   median(x[[i]][[j]][[k]]),x[[i]][[j]][[k]][1],
						   mean(x[[i]][[j]][[k]]),
						   quantile(x[[i]][[j]][[k]],.95)));
					out<-rbind(out,c(quantile(attr(x[[i]][[j]],'ci')[[k]],.05),
						   median(attr(x[[i]][[j]],'ci')[[k]]),
						   attr(x[[i]][[j]],'ci')[[k]][1],
						   mean(attr(x[[i]][[j]],'ci')[[k]]),
						   quantile(attr(x[[i]][[j]],'ci')[[k]],.95)));
					ijkname<-paste(ijnames$model,ijnames$constraint,
						       ijnames$parameters[k]);
					rnames<-c(rnames,paste(c('Permuted','Resampled'),ijkname))
				}
			}
		}
		colnames(out)<-c('5%','median','actual','mean','95%');
		rownames(out)<-rnames;
		oldwidth<-options('width')[[1]];options('width'=130);
		print(out);options('width'=oldwidth);invisible(out);
	}
}

`print.spymodels`<-
function(x,type=c('plot','table'),...){
	if(match('plot',type,nomatch=0)){
		layout(matrix(c(1:5,0),2,3,byrow=T));
		for(i in 1:5){
			model<-switch(i,'Gompertz vs. Logistic',
					'Gompertz vs. Gompertz-Makeham',
					'Logistic vs. Logistic-Makeham',
					'Gompertz-Makeham vs. Logistic-Makeham',
					'Gompertz vs. Logistic-Makeham');
			iecdf<-ecdf(x[[i]]);
			ip<-1-iecdf(x[[i]][1]);
			isub<-paste('p =', ip);
			hist(x[[i]],breaks=100,xlab='',ylab='',main=model,freq=F,sub=isub,...)
			abline(v=x[[i]][1],col='red');
			abline(v=attr(x[[i]],'ci')[1],col='green',lty=2);
			hist(attr(x[[i]],'ci'),breaks=100,freq=F,border='pink',lty=2,add=T);
		}
		layout(matrix(1,1,1));
	}
	if(match('table',type,nomatch=0)){
		out<-c();
		compars<-c('G vs. L', 'G vs. GM', 'L vs. LM', 'GM vs. LM', 'G vs. LM');
		compars<-expand.grid(c('Permuted','Resampled'),compars);
		compars<-paste(compars[,1],compars[,2]);
		for(i in 1:length(x)){
			if(i!=5){idf<-1;}else{idf<-2;}
			out<-rbind(out,c(quantile(x[[i]],.05),median(x[[i]]),x[[i]][1],mean(x[[i]]),
					 quantile(x[[i]],.95),1-ecdf(x[[i]])(x[[i]][1]),
					 pchisq(x[[i]][i],idf,lower=F)));
			out<-rbind(out,c(quantile(attr(x[[i]],'ci'),.05),median(attr(x[[i]],'ci')),
					 attr(x[[i]],'ci')[1],mean(attr(x[[i]],'ci')),
					 quantile(attr(x[[i]],'ci'),.95),NA,NA));
		}
		colnames(out)<-c('5%','median','actual','mean','95%','emp p','chisq p');
		rownames(out)<-compars;
		oldwidth<-options('width')[[1]];options('width'=130);
		print(out);options('width'=oldwidth);invisible(out);
	}
}

`spyhaz2`<-
function(logfile,realval=NULL){
	env<-new.env();sys.source(logfile,env);logtemp<-as.list(env);
	out<-list(); out$log<-logtemp;
	out$loglen<-loglen<-length(logtemp);
	if(is.null(realval)){out$realval<-realval<-logtemp[[loglen]];}
	out$realmax<-realmax<-realval$maximum;
	ests<-maxs<-iters<-rts<-c();
	for(i in logtemp){
		ests<-rbind(ests,i$estimate);
		maxs<-c(maxs,i$maximum);iters<-c(iters,i$titer);rts<-c(rts,i$runtime);
	};
	out$ests<-ests; out$maxs<-maxs; out$iters<-iters; out$rts<-rts;
	lests<-dim(ests)[2];
	hist(maxs,breaks=1000,xlim=c(min(maxs),max(maxs)),freq=F);
	abline(v=realmax,col='red');myecdf<-ecdf(maxs);
	out$permp<-1-myecdf(realmax); out$avgrt<-sum(rts)/loglen; out$loglen<-loglen;
	maxests<-minests<-means<-gmeans<-medians<-numeric(lests);
	for(i in 1:lests){
		means[i]<-mean(ests[,i]);medians[i]<-median(ests[,i]);gmeans[i]<-gmean(ests[,i]);
		maxests[i]<-max(ests[,i]);minests[i]<-min(ests[,i]);
	}
	avpars<-rbind(means,medians,gmeans,maxests,minests);
	out$avpars<-avpars;
	class(out)<-'spyhaz';
	print(out);
	invisible(out);
}

`print.spyhaz2`<-
function(x,...){
	cat('permuted p:',x$permp,'\n');
	cat('avg runtime:',x$avgrt,'\n');
	cat('# perms:',x$loglen,'\n');
	print(x$avpars);
}

`raw2wm`<-
function(t){
	if(!is.list(t)){t<-list(t);}
	out<-c();
	for(i in 1:length(t)){
		irle<-rle(sort(t[[i]]));
		iwm<-cbind(irle$values,irle$lengths,i,1);
		out<-rbind(out,iwm);
	}
	out<-rbind(out,c(-1,NA,NA,NA));
	colnames(out)<-NULL;
	out;
}

`raw2tab`<-
function(x,y,labels=c(0,1)){cbind(c(x,y),c(rep(labels[1],length(x)),rep(labels[2],length(y))));}

`int2tab`<-
function(x,y,labels=c(0L,1L))
	{cbind(as.integer(c(x,y)),c(rep(labels[1],length(x)),rep(labels[2],length(y))));}
