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

`startpars.old`<-
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

`ressurv.old` <-
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
		else{pars<-startpars.old(x,y,label,perms=10);}
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

`findpars`<-
function(x,y,nil=1e-7,bnil=0.005,pf=mean,labels=NULL,summary=F,id=0){
	lx<-length(x);ly<-length(y);sumx<-sum(x);sumy<-sum(y);
	xemax<-lx*log(lx/sumx)-lx;yemax<-ly*log(ly/sumy)-ly;
	es.e<-list(maximum=xemax+yemax,
		   estimate=c(lx/sumx,ly/sumy),group0=list(maximum=xemax),
		   group1=list(maximum=yemax),titer=0,runtime=0);
	partemp<-modpars(es.e$estimate,'e','g',nil=nil,trim=T);
	partemp[c(2,4)]<-bnil;
	es.g<-opsurv(x,y,'g',par=partemp,minout=es.e$maximum);
	es.ga<-opsurv(x,y,'g',cons=c(F,T,T,T),
		      par=modpars(es.g$estimate,'g',cno=c(F,T,T,T),pf=pf,trim=T));
	if(es.ga$maximum>es.g$maximum){
		es.g<-opsurv(x,y,'g',
			     par=modpars(es.ga$estimate,'g',cni=c(F,T,T,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.gb<-opsurv(x,y,'g',cons=c(T,F,T,T),
		      par=modpars(es.g$estimate,'g',cno=c(T,F,T,T),pf=pf,trim=T));
	if(es.gb$maximum>es.g$maximum){
		es.g<-opsurv(x,y,'g',
			     par=modpars(es.gb$estimate,'g',cni=c(T,F,T,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	partemp<-modpars(es.g$estimate,'g','gm',nil=nil,trim=T);
	es.gm<-opsurv(x,y,'gm',par=partemp);
	es.gma<-opsurv(x,y,'gm',cons=c(F,T,T,T),
		       par=modpars(es.gm$estimate,'gm',cno=c(F,T,T,T),pf=pf,trim=T));
	if(es.gma$maximum>es.gm$maximum){
		es.gm<-opsurv(x,y,'gm',
			     par=modpars(es.gma$estimate,'gm',cni=c(F,T,T,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.gmb<-opsurv(x,y,'gm',cons=c(T,F,T,T),
		       par=modpars(es.gm$estimate,'gm',cno=c(T,F,T,T),pf=pf,trim=T));
	if(es.gmb$maximum>es.gm$maximum){
		es.gm<-opsurv(x,y,'gm',
			     par=modpars(es.gmb$estimate,'gm',cni=c(T,F,T,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.gmc<-opsurv(x,y,'gm',cons=c(T,T,F,T),
		       par=modpars(es.gm$estimate,'gm',cno=c(T,T,F,T),pf=pf,trim=T));
	if(es.gmc$maximum>es.gm$maximum){
		es.gm<-opsurv(x,y,'gm',
			     par=modpars(es.gmc$estimate,'gm',cni=c(T,T,F,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.l<-opsurv(x,y,'l',par=partemp);
	es.la<-opsurv(x,y,'l',cons=c(F,T,T,T),
		      par=modpars(es.l$estimate,'l',cno=c(F,T,T,T),pf=pf,trim=T));
	if(es.la$maximum>es.l$maximum){
		es.l<-opsurv(x,y,'l',
			     par=modpars(es.la$estimate,'l',cni=c(F,T,T,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.lb<-opsurv(x,y,'l',cons=c(T,F,T,T),
		      par=modpars(es.l$estimate,'l',cno=c(T,F,T,T),pf=pf,trim=T));
	if(es.lb$maximum>es.l$maximum){
		es.l<-opsurv(x,y,'l',
			     par=modpars(es.lb$estimate,'l',cni=c(T,F,T,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.ls<-opsurv(x,y,'l',cons=c(T,T,T,F),
		      par=modpars(es.l$estimate,'l',cno=c(T,T,T,F),pf=pf,trim=T));
	if(es.ls$maximum>es.l$maximum){
		es.l<-opsurv(x,y,'l',
			     par=modpars(es.ls$estimate,'l',cni=c(T,T,T,F),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.gmlm<-opsurv(x,y,'lm',
		      par=modpars(es.gm$estimate,'gm','lm',nil=nil,trim=T));
	es.llm<-opsurv(x,y,'lm',
		      par=modpars(es.l$estimate,'l','lm',nil=nil,trim=T));
	# Might wish to consider getting rid of es.gmlm and have the gm output go straight to es.lm
	if(es.gmlm$maximum<es.llm$maximum){es.lm<-es.llm;} else {es.lm<-es.gmlm;}
	es.lma<-opsurv(x,y,'lm',cons=c(F,T,T,T),minout=es.lm$maximum,
		       par=modpars(es.lm$estimate,'lm',cno=c(F,T,T,T),pf=pf,trim=T));
	if(es.lma$maximum>es.lm$maximum){
		es.lm<-opsurv(x,y,'lm',
			     par=modpars(es.lma$estimate,'lm',cni=c(F,T,T,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.lmb<-opsurv(x,y,'lm',cons=c(T,F,T,T),minout=es.lm$maximum,
		       par=modpars(es.lm$estimate,'lm',cno=c(T,F,T,T),pf=pf,trim=T));
	if(es.lmb$maximum>es.lm$maximum){
		es.lm<-opsurv(x,y,'lm',
			     par=modpars(es.lmb$estimate,'lm',cni=c(T,F,T,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.lmc<-opsurv(x,y,'lm',cons=c(T,T,F,T),minout=es.lm$maximum,
		       par=modpars(es.lm$estimate,'lm',cno=c(T,T,F,T),pf=pf,trim=T));
	if(es.lmc$maximum>es.lm$maximum){
		es.lm<-opsurv(x,y,'lm',
			     par=modpars(es.lmc$estimate,'lm',cni=c(T,T,F,T),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	es.lms<-opsurv(x,y,'lm',cons=c(T,T,T,F),minout=es.lm$maximum,
		       par=modpars(es.lm$estimate,'lm',cno=c(T,T,T,F),pf=pf,trim=T));
	if(es.lms$maximum>es.lm$maximum){
		es.lm<-opsurv(x,y,'lm',
			     par=modpars(es.lms$estimate,'lm',cni=c(T,T,T,F),
					 cno=c(T,T,T,T),pf=pf,trim=T));}
	if(!is.null(labels)){
		write.table(t(c(es.e$max,NA,NA,modpars(es.e$estimate,'e',nil=nil,trim=F),id,'e',nil,NA)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.g$max,es.g$titer,es.g$runtime,
				modpars(es.g$estimate,'g',nil=nil,trim=F),id,'g',nil,es.g$max-es.e$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.ga$max,es.ga$titer,es.ga$runtime,
				modpars(es.ga$estimate,'g',cni=c(F,T,T,T),nil=nil,trim=F),id,'ga',nil,
			        es.g$max-es.ga$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.gb$max,es.gb$titer,es.gb$runtime,
				modpars(es.gb$estimate,'g',cni=c(T,F,T,T),nil=nil,trim=F),id,'gb',nil,
			        es.g$max-es.gb$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.gm$max,es.gm$titer,es.gm$runtime,
				modpars(es.gm$estimate,'gm',nil=nil,trim=F),id,'gm',nil,
			        es.gm$max-es.g$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.gma$max,es.gma$titer,es.gma$runtime,
				modpars(es.gma$estimate,'gm',cni=c(F,T,T,T),nil=nil,trim=F),
				id,'gma',nil,es.gm$max-es.gma$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.gmb$max,es.gmb$titer,es.gmb$runtime,
				modpars(es.gmb$estimate,'gm',cni=c(T,F,T,T),nil=nil,trim=F),
				id,'gmb',nil,es.gm$max-es.gmb$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.gmc$max,es.gmc$titer,es.gmc$runtime,
				modpars(es.gmc$estimate,'gm',cni=c(T,T,F,T),nil=nil,trim=F),
				id,'gmc',nil,es.gm$max-es.gmc$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.l$max,es.l$titer,es.l$runtime,
				modpars(es.l$estimate,'l',nil=nil,trim=F),id,'l',nil,es.l$max-es.g$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.la$max,es.la$titer,es.la$runtime,
				modpars(es.la$estimate,'l',cni=c(F,T,T,T),nil=nil,trim=F),id,'la',nil,
			        es.l$max-es.la$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.lb$max,es.lb$titer,es.lb$runtime,
				modpars(es.lb$estimate,'l',cni=c(T,F,T,T),nil=nil,trim=F),id,'lb',nil,
			        es.l$max-es.lb$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.ls$max,es.ls$titer,es.ls$runtime,
				modpars(es.ls$estimate,'l',cni=c(T,T,T,F),nil=nil,trim=F),id,'ls',nil,
			        es.l$max-es.ls$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.lm$max,es.lm$titer,es.lm$runtime,
				modpars(es.lm$estimate,'lm',nil=nil,trim=F),id,'lm',nil,
			        es.lm$max-max(es.l$max,es.gm$max))),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.lma$max,es.lma$titer,es.lma$runtime,
				modpars(es.lma$estimate,'lm',cni=c(F,T,T,T),nil=nil,trim=F),
				id,'lma',nil,es.lm$max-es.lma$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.lmb$max,es.lmb$titer,es.lmb$runtime,
				modpars(es.lmb$estimate,'lm',cni=c(T,F,T,T),nil=nil,trim=F),
				id,'lmb',nil,es.lm$max-es.lmb$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.lmc$max,es.lmc$titer,es.lmc$runtime,
				modpars(es.lmc$estimate,'lm',cni=c(T,T,F,T),nil=nil,trim=F),
				id,'lmc',nil,es.lm$max-es.lmc$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.lms$max,es.lms$titer,es.lms$runtime,
				modpars(es.lms$estimate,'lm',cni=c(T,T,T,F),nil=nil,trim=F),
				id,'lms',nil,es.lm$max-es.lms$max)),
			    file=labels$main,sep='\t',append=T,col.names=F,row.names=F);

		write.table(t(c(es.e$group0$maximum,es.e$group1$maximum,id,'e',nil)),
			    file=labels$unc,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.g$group0$maximum,es.g$group1$maximum,id,'g',nil)),
			    file=labels$unc,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.gm$group0$maximum,es.gm$group1$maximum,id,'gm',nil)),
			    file=labels$unc,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.l$group0$maximum,es.l$group1$maximum,id,'l',nil)),
			    file=labels$unc,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.lm$group0$maximum,es.lm$group1$maximum,id,'lm',nil)),
			    file=labels$unc,sep='\t',append=T,col.names=F,row.names=F);

		write.table(t(c(es.gmlm$max,es.gmlm$titer,es.gmlm$runtime,
				modpars(es.gmlm$estimate,'lm',trim=F),id,'gmlm',nil)),
			    file=labels$gml,sep='\t',append=T,col.names=F,row.names=F);
		write.table(t(c(es.llm$max,es.llm$titer,es.llm$runtime,
				modpars(es.llm$estimate,'lm',trim=F),id,'llm',nil)),
			    file=labels$gml,sep='\t',append=T,col.names=F,row.names=F);
	}
	if(summary){
		maxs<-c(es.e$max,es.g$max,es.ga$max,es.gb$max,es.gm$max,es.gma$max,es.gmb$max,es.gmc$max,
			es.l$max,es.la$max,es.lb$max,es.ls$max,es.gmlm$max,es.llm$max,es.lma$max,
			es.lmb$max,es.lmc$max,es.lms$max);
		iters<-c(0,es.g$titer,es.ga$titer,es.gb$titer,es.gm$titer,es.gma$titer,es.gmb$titer,
			es.gmc$titer,es.l$titer,es.la$titer,es.lb$titer,es.ls$titer,es.gmlm$titer,
			es.llm$titer,es.lma$titer,es.lmb$titer,es.lmc$titer,es.lms$titer);
		difs<-sign(c(0,es.g$max-es.e$max,es.g$max-es.ga$max,es.g$max-es.gb$max,es.gm$max-es.g$max,
			es.gm$max-es.gma$max,es.gm$max-es.gmb$max,es.gm$max-es.gmc$max,es.l$max-es.g$max,
			es.l$max-es.la$max,es.l$max-es.lb$max,es.l$max-es.ls$max,es.gmlm$max-es.gm$max,
			es.llm$max-es.l$max,es.lm$max-es.lma$max,es.lm$max-es.lmb$max,
			es.lm$max-es.lmc$max,es.lm$max-es.lms$max));
		out<-cbind(maxs,iters,difs);
		#rownames(ests)<-
		rownames(out)<-c('E','G','Ga','Gb','GM','GMa','GMb','GMc','L','La','Lb','Ls',
						'GMLM','LLM','LMa','LMb','LMc','LMs');
		print(out,digits=22);#print(ests);
		print(sum(out[,2]));invisible(list(results=out));
	}
}

`modpars`<-
function(x,modeli,modelo=NULL,cni=rep(T,4),cno=NULL,nil=1e-7,pf=mean,trim=F){
	if(sum(is.null(x))>0){print('input params contain nulls!');traceback();browser();};
	if(is.null(modelo)){modelo<-modeli;}; if(is.null(cno)){cno<-cni;};
	keepi<-switch(modeli,e=c(1,5),g=c(1,2,5,6),gm=c(1:3,5:7),l=c(1,2,4:6,8),lm=1:8);
	keepo<-switch(modelo,e=c(1,5),g=c(1,2,5,6),gm=c(1:3,5:7),l=c(1,2,4:6,8),lm=1:8);
	# Will eventually do more input validation, but for now let's agree to pass 
	# constraints as vectors of either 4 or 1 logical values
	cni <- !cni; if(length(cni)==1){cni<-rep(cni,4);};
	cno <- !cno; if(length(cno)==1){cno<-rep(cno,4);};
	if(length(cni)+length(cno)!=8){stop('In this version constraints must be specified
either as a single boolean value or a vector
of four boolean values. Please fix your input
and try again.');}
	keepi<-keepi[c(rep(T,4),!cni)[keepi]];
	if(length(x)!=length(keepi)){
		print('Mismatch between model and parameter length.');traceback();browser();
	}
	# Here we convert the unique-values-only parameter vector into standard form
	out<-rep.int(NA,8); out[keepi]<-x; 
	if(sum(cni)>0){out[5:8][cni]<-out[1:4][cni];}
	# We insure that any parameters used by the target model are populated
	out[keepo][is.na(out[keepo])]<-nil;
	# Now we apply the target model's constraints
	if(sum(cno)>0){
		# If the target constraint is different from the input constraint,
		# average the two parameters using the function specified by pf
		# No point in doing this is the constraints are the same
		if(sum(cni==cno)<4&sum(out[cno]>0)){
			out[cno]<-apply(out*as.matrix(diag(cno)[c(1:4,1:4),cno]),2,
				       function(x,f){f(na.omit(x[x!=0]));},pf);
		}
		if(trim){out[5:8][cno]<-NA;}
	}
	# Gotta remember to trash the unused parameters for the target model.
	out[-keepo]<-NA;
	if(trim){out<-as.numeric(na.omit(out));} 
	return(out);
}

`compsurvs`<-
function(a,b){
	out1<-cbind(c(a$maximum,a$titer),c(b$maximum,b$titer));
	rownames(out1)<-c('maximum','titer');
	out2<-rbind(a$estimate,b$estimate);
	print(out1,digits=22);print(out2);
	invisible(list(out1,out2));
}


`empdist` <-
function(x,y,label,nil=2e-7,pf=mean,rounds=5000,type='p'){
	d<-int2tab(x,y); lx<-length(x); ly<-length(y); lxy<-lx+ly;
	labels<-list(main=paste(label,'.main.txt',sep=''),
		     unc=paste(label,'.unc.txt',sep=''),
		     gml=paste(label,'.gml.txt',sep=''));
	if(type=='p'){
		m<-matrix(nrow=lxy,ncol=rounds);
		for(i in 1:rounds){m[,i]<-sample(d[,2],lxy);}
		m<-cbind(d[,2],m);
		write.table(m,paste(label,'_perms.txt',sep=''),col.names=F,row.names=F,sep='\t');
		cat('\nCommencing permutations.\n');
		for(i in 1:(rounds+1)){
			findpars(d[m[,i]==0,1],d[m[,i]==1,1],nil=nil,pf=pf,labels=labels,id=i);
		}
	}
	if(type=='b'){
		rsx<-matrix(nrow=lx,ncol=rounds);
		rsy<-matrix(nrow=ly,ncol=rounds);
		for(i in 1:rounds){
			rsx[,i]<-sample(x,lx,replace=T);
			rsy[,i]<-sample(y,ly,replace=T);
		}
		rsx<-cbind(x,rsx);rsy<-cbind(y,rsy);
		write.table(rsx,paste(label,'_xboot.txt',sep=''),col.names=F,row.names=F,sep='\t');
		write.table(rsy,paste(label,'_yboot.txt',sep=''),col.names=F,row.names=F,sep='\t');
		cat('\nCommencing bootstrap.\n');
		for(i in 1:(rounds+1)){
			findpars(rsx[,i],rsy[,i],nil=nil,pf=pf,labels=labels,id=i);
		}
	}
}


`permsurv.old` <-
function(x,y,label,perms=5000,pars=NULL){
	d<-int2tab(x,y); lxy<-dim(d)[1]; pr<-matrix(nrow=lxy,ncol=perms);
	for(i in 1:perms){pr[,i]<-sample(d[,2],lxy);}
	pr<-cbind(d[,2],pr);
	if(is.character(pars)){parslabel<-pars;pars<-NULL;}
	else{parslabel<-paste(label,'pars',sep='');cat('\nLooking for',parslabel,'\n');}
	# need to write some kind of validation function for fits
	if(is.null(pars)|length(pars)!=4){
		if(exists(parslabel)){cat('\nFound',parslabel,'\n');pars<-get(parslabel);}
		else{pars<-startpars.old(x,y,label,perms=10);}
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


`hazwrapper.old` <-
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

`dropn`<-
function(x,n,which='smallest'){
	if(match('largest',which,nomatch=0)){
		for(i in 1:n){
			x<-x[-which.max(x)];
		}
	}
	if(match('largest',which,nomatch=0)){
		for(i in 1:n){
			x<-x[-which.min(x)];
		}
	}
}

`spyhaz`<-
function(label,path='./',extramain=0){
	main<-paste(label,'main',sep='.');
	gml<-paste(label,'gml',sep='.');
	unc<-paste(label,'unc',sep='.');
	gmlcols<-maincols<-c("max", "iters", "rt", "a1", "b1", "c1", "s1", "a2", "b2", "c2",
			      "s2", "id", "fit", "nilval");
	unccols<-c('max1','max2','id','fit','nilval');
	length(maincols)<-length(maincols)+extramain;
	assign(main,read.table(paste(path,main,'.txt',sep=''),col.names=maincols),envir=.GlobalEnv);
	assign(gml,read.table(paste(path,gml,'.txt',sep=''),col.names=gmlcols),envir=.GlobalEnv);
	assign(unc,read.table(paste(path,unc,'.txt',sep=''),col.names=unccols),envir=.GlobalEnv);
	cat('Output collected from',max(get(gml)$id),'iterations.\n');
}

`plot.mainhaz`<-
function(perm,breaks=100,add=F,col=1,abcol='red',onescreen=T,what='cons'){
	if(what=='cons'){
		if(onescreen){layout(matrix(1:12,4,3,byrow=T));}
		plotlr(perm,'g','ga',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'g','gb',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'gm','gma',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'gm','gmb',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'gm','gmc',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'l','la',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'l','lb',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'l','ls',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','lma',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','lmb',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','lmc',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','lms',breaks=breaks,col=col,abcol=abcol,add=add);
	}
	if(what=='model'){
		if(onescreen){layout(matrix(1:6,3,2,byrow=T));}
		plotlr(perm,'g','e',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'gm','g',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'l','g',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','gm',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','l',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','g',breaks=breaks,col=col,abcol=abcol,add=add,df=2);
	}
	if(onescreen){layout(matrix(1,1,1));}
}

`plotlr`<-
function(perm,fit1='gm',fit2='g',breaks=100,df=1,col=1,abcol='red',add=F){
	lrs<-2*(subset(perm,fit==fit1)$max-subset(perm,fit==fit2)$max);
	badlrs<-lrs<0; sumbad<-sum(badlrs);
	cat(fit1,'vs',fit2,'fraction <0:',sumbad/length(lrs),'\n');
	print(summary(lrs[badlrs]));
	if(sumbad<10){cat(lrs[badlrs],'\n');} else {cat(lrs[badlrs][1:10],'...\n');}
	empp<-1-ecdf(lrs)(lrs[1]);
	chip<-pchisq(lrs[1],df,lower=F);
	hist(lrs,breaks=breaks,border=col,main=paste(fit1,'vs',fit2),sub='',add=add,
	     xlab=paste('Emp. p:',format(empp,digits=3),'  Chi-sq p:',format(chip,digits=3)));
	abline(v=lrs[1],col=abcol);
}


`plotest`<-
function(perm,val='gm',par='a1',breaks=100,add=F,col=1){
	ests<-subset(perm,fit==val)[[par]];
	hist(ests,breaks=breaks,add=add,main=paste(par,'parameter for the',val,'fit'),border=col);
	abline(v=ests[1],col='red');
}

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

`spyhaz.old`<-
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
					     freq=T,main=model,sub=ijsub,...);
					abline(v=x[[i]][[j]][1],col='red');
					abline(v=attr(x[[i]][[j]],'ci')[1],col='green',lty=2);
					hist(attr(x[[i]][[j]],'ci'),breaks=100,freq=T,border='red',
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
			hist(x[[i]],breaks=100,xlab='',ylab='',main=model,freq=T,sub=isub,...)
			abline(v=x[[i]][1],col='red');
			abline(v=attr(x[[i]],'ci')[1],col='green',lty=2);
			hist(attr(x[[i]],'ci'),breaks=100,freq=T,border='red',lty=2,add=T);
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
