.First<-function(x=NULL,skip=F){
	if(skip){invisible();}
	updpath='http://biochem2.uthscsa.edu/~bokov/survomatic.';
	welcome<-paste("Welcome ",if(!is.null(x)){'back '},"to Survomatic, the survival analysis package\nthat at least tries to be user-friendly.\nHow can I help you today?",sep="");
	choices<-"Analyze brand new data (you can do this manually by typing go()).";
	if(!is.null(x)){
		choices<-c(choices,
		paste("Print out the results of the previous analysis (you can do this manually by typing ",
		paste(x,"$go()",sep=""),").",sep=""));
	}
	choices<-c(choices,"Update Survomatic","Let me do my own thing, I know what I want.");
	input<-menu(choices,graphics=T,title=welcome);
	if(is.null(x)&input>1){input<-input+1;}
	switch(input,
		go(),
		{d<-get(x);d$go();},
		{
			options(warn=-1);options(show.error.messages=F);
			upd<-try(source(updpath),silent=T);
			if(class(upd)=='try-error'){
				cat("\nUpdate failed. Please contact bokov@uthscsa.edu. But in the meantime you can still continue using the current version.\n");
				.First();
			};
			options(warn=0);options(show.error.messages=T);invisible();
		},
		{cat('\n');invisible();}
	);
}


instpx<-function(libs){
	missing<-c();
	for(i in libs){if(!require(i,character.only=T,warn.conflicts=F)){missing<-c(missing,i)}}
	return(missing);
}

checkpx<-function(libs){
	missing<-instpx(libs);
	if(length(missing)>0){
		message<-"You are missing one or more needed R packages. If you are connected to the internet, I can automatically install them for you. Do you wish to do that?";
		choices<-c('Yes','No (quit program)');
		input<-menu(choices,graphics=T,title=message);
		switch(input,
			install.packages(missinglib,repos='http://cran.fhcrc.org'),
			{cat('When you do get this computer connected to the internet, please try running this script again in order to install the missing packages:',paste(missing,collapse=""),'\n');return(-1);}
		)
		missing<-instpx();
		if(length(missing)>0){cat("For some reason, the following needed packages cannot be installed: ",paste(missing, collapse=","),"\nMake sure you can connect to the internet and try this again. If you know for certain that your internet connection works and the packages still don\'t install, you need to talk to the authors of Survomatic and/or your local R expert to get this problem fixed.\n");return(-1);} else {return(1);}
	} else {return(1);}
}

go<-function(x,y,names=c(),units="d",save=T,path,prompt=2,xlim=1350){
	options(warn=-1);options(show.error.messages=F);
	libs=c('SparseM','splines','survival','quantreg','surv2sample','tcltk');
	if(checkpx(libs)<0){return(-1);}
	quantiles<-list(seq(.01,.99,.01),seq(.05,.95,.05),seq(.1,.9,.1),seq(.2,1,.2));
	mycall<-strsplit(as.character(sys.call()[1]),"\\$")[[1]];
	if(length(mycall)>1){
		top=F;d<-get(mycall[1]);
		if(length(mycall)>2){
			for(i in 2:(length(mycall)-1)){
				d<-d[[mycall[i]]];
		}}
		x<-d$x;y<-d$y;names<-d$names;smry<-d$smry;zsc<-d$zsc;xy<-d$xy;
		group<-d$group;lr<-d$lr;qreg<-d$qreg;qreg.sum<-d$qreg.sum;tt<-d$tt;
		qreg.tab<-d$qreg.tab;xmod<-d$xmod;ymod<-d$ymod;xymod<-d$xymod;
		xd1<-d$xd1;xd7<-d$xd7;xd30<-d$xd30;yd1<-d$yd1;yd7<-d$yd7;yd30<-d$yd30;
		xhzpars<-d$xhzpars;yhzpars<-d$yhzpars;xymod.tab<-d$xymod.tab;
		modx<-d$modx;mody<-d$mody;modxy<-d$modxy;xysurvfit<-d$xysurvfit;
		sig.tests<-d$sig.tests;tests<-d$tests; report<-d$report;
	}else{top=T;}
	if(missing(path)&top&save){
		message<-"Where would you like to save your results?";
		choices<-c("Let me choose using the point and click method.",
			   "Let me type in the location from the console.",
			   "I don\'t want to save any files. Run analysis without saving.");
		input<-menu(choices,graphics=T,title=message);
		switch(input,
			{path<-paste(tclvalue(tkchooseDirectory()),"/",sep="");},
			{cat("\nPlease type in the location where you would like to save files:\n");
			   path<-paste(scan("",'character',nlines=1,quiet=T),"/",sep="");},
			{save<-F;cat('\nContinuing without saving.');}
		)
		if(save){cat('\nResults will be saved to:',path);}
	}
	if(missing(x)&top){
		cat("\nYou need a group to compare. Please copy and paste one column of survival times from a spreadsheet.
They don\'t have to be in order, but make sure there are no empty spaces. Hit the enter key when you 
are done.\n");
		x<-scan();	
	}
	if(length(names)<2&top){
		cat("Please give a brief name for the first group.\n");
		names[1]<-paste(scan("",'character',nlines=1,quiet=T),collapse=".");
	}
	if(missing(y)&top){
		cat("\nYou need a second group to compare. Please copy and paste one column of survival times from a 
spreadsheet. They don\'t have to be in order, but make sure there are no empty spaces. Hit 
the enter key when you are done.\n");
		y<-scan();
	}
	if(length(names)<2&top){
		cat("Please give a brief name for the second group.\n");
		names[2]<-paste(scan("",'character',nlines=1,quiet=T),collapse=".");
	}
	if(save&top){
		assign(names[1],x,envir=.GlobalEnv);
		assign(names[2],y,envir=.GlobalEnv);
		cat("The survival times for the two groups have been saved as",names[1],"and",names[2],"inside R so you won\'t have to paste them in next time. Instead you can type \'go(",names[1],",",names[2],")\'\n");
	}

	if(top){
		tests<-c("t test"="tt","z score"="zsc","log rank"="lr","quantile regression"="qreg.tab",
			"mortality parameters"="xymod.tab");
		sig.tests<-c();
	}

	cat("\nDoing summary statistics...");
	if(top){smry<-t(cbind(summary(x),summary(y))); rownames(smry)<-names;}
	cat('\n'); print(smry);
	cat("\nDoing t-test on log-transformed data...");
	if(top){tt<-t.test(log(x),log(y));if(tt$p.value<=.05) sig.tests<-c(sig.tests,1);}
	cat('\n'); print(tt);
	cat("\nDoing z-score test to see which quantiles differ..."); 
	if(top){
		zsc<-c();
		zsc<-ezz(x,y,names,quant=seq(0,1,.05));
		if(length(levels(zsc$sig))>1) sig.tests<-c(sig.tests,2);
	} 
	cat('\n'); print(zsc);
	cat("\nDoing log-rank test with permutations...");
	if(top){
		xy<-c(x,y); group<-c(rep(names[1],length(x)),rep(names[2],length(y)));
		print(group);
		lr<-surv2.logrank(Surv(xy),factor(group));
		if(lr$pval<=.05) sig.tests<-c(sig.tests,3);
	}
	print(lr);
	cat("\nDoing quantile regression...");
	if(top){
		i<-1;qreg.sum<-0;class(qreg.sum)<-"try-error";
		while(class(qreg.sum)=="try-error"){
			if(i>length(quantiles)){
				qreg.sum<-"Unable to perform quantile regression.";
				next;
			}
			#cat('\nTrying quantiles...',paste(quantiles[[i]],collapse=" "));
			qreg<-rq(log(xy)~group,tau=quantiles[[i]]);i<-i+1;
			qreg.sum<-try(summary(qreg,se='ker'));
		}
		if(class(qreg.sum)!="character"){
			qreg.tab<-c();
			for(i in qreg.sum){
				qreg.tab<-rbind(qreg.tab,c(quantile=i$tau,i$coefficients[2,]));
				if(min(qreg.tab[,5])<=.05) sig.tests<-c(sig.tests,4);
			}}else{qreg.tab<-qreg.sum;}
	}
	cat('\n'); print(qreg.tab);

	cat("\nFinding the best mortality models..."); 
	if(top){xmod<-ymod<-xymod<-list();
		xd1<-demog(x); xd7<-demog(x,7); xd30<-demog(x,30);
		yd1<-demog(y); yd7<-demog(y,7); yd30<-demog(y,30);
		xmod$g<-opsurv('g',x,par=rep(5e-7,4)); ymod$g<-opsurv('g',y,par=rep(5e-7,4));
		xmod$gm<-opsurv('gm',x,par=rep(5e-7,4)); ymod$gm<-opsurv('gm',y,par=rep(5e-7,4));
		xmod$l<-opsurv('l',x,par=rep(5e-7,4)); ymod$l<-opsurv('l',y,par=rep(5e-7,4));
		xmod$lm<-opsurv('lm',x,par=rep(5e-7,4)); ymod$lm<-opsurv('lm',y,par=rep(5e-7,4));
		likxu<-likyu<-dlikxu<-dlikyu<-c();
		for(i in 1:4){
			if(class(xmod[[i]])[1]=="try-error"){likxu[i]<- -Inf;}
			else{likxu[i]<-xmod[[i]]$maximum;}
			if(class(ymod[[i]])[1]=="try-error"){likyu[i]<- -Inf;}
			else{likyu[i]<-ymod[[i]]$maximum;}
		}
		for(i in 2:4){
			dlikxu[i-1]<- -2*(likxu[1]-likxu[i]);
			dlikyu[i-1]<- -2*(likyu[1]-likyu[i]);
		}
		dlikxu<-pchisq(dlikxu,c(1,1,2),lower=F);
		dlikyu<-pchisq(dlikyu,c(1,1,2),lower=F);
		if(min(dlikxu)<=.05){modx<-max(which(dlikxu<=.05))+1;} else {modx<-1;}
		if(min(dlikyu)<=.05){mody<-max(which(dlikyu<=.05))+1;} else {mody<-1;}
		modxy<-min(modx,mody);
	}
	cat("\nThe best model for",names[1],"is",names(xmod)[modx],".");
	cat("\nThe best model for",names[2],"is",names(ymod)[mody],".");
	cat("\nChecking for shared mortality parameters between groups...");
	if(top){
		idxcons<-switch(modxy,c(a=1,b=2),c(a=1,b=2,c=3),c(a=1,b=2,s=4),c(a=1,b=2,c=3,s=4));
		for(i in idxcons){
			cons<-rep(1,4); cons[i]<-0;
			xymod[[names(idxcons)[i]]]<-opsurv(names(xmod)[modxy],x,y,rep(5e-7,4),cons);
		}
		xymod.tab<-c();
		for(i in xymod){
			xymod.tab<-rbind(xymod.tab, c(i$estimate,i$maximum,
					pchisq(-2*(likxu[modxy]+likyu[modxy]-i$maximum),1,lower=F)));
		}
		rownames(xymod.tab)<-names(xymod);
		colnames(xymod.tab)<-c(rep("",dim(xymod.tab)[2]-2),"Likelihood","p");
		if(xymod.tab[,5]<=.05) sig.tests<-c(sig.tests,5);
	}
	cat('\n');print(xymod.tab); 

	if(top){xysurvfit<-survfit(Surv(xy)~group);}
	tmain<-paste(names,collapse=" vs. ");
	cat("\nPlotting survival curves.");
	plot(xysurvfit,lwd=1,bty="n",xlim=c(0,xlim));
	for(i in 1:2){points(xysurvfit[i],pch=22+i,cex=1.5,lwd=1.5,bg=rainbow(2)[i]);}
	title(main=paste(tmain,"Survival"));
	legend("topright",legend=names[2:1],pch=23:24,y.intersp=1.5,
		#inset=-0.006*psize,-+
		bty="n",pt.cex=1.5,cex=1.5,lwd=1.5,pt.bg=rainbow(2));
	
	graphsave(path,names,".surv.eps",prompt);
	if(class(qreg.sum)!="character"){
		cat("\n\nPlotting quantile regression curves.");
		plot(qreg.sum);
		graphsave(path,names,".qreg.eps",prompt);
	title(main=paste(tmain,"Quantile Regression"));
	}

	if(top){
		xhzpars<-rbind(c(xmod$g$estimate,NA),c(xmod$gm$estimate));
		yhzpars<-rbind(c(xmod$g$estimate,NA),c(xmod$gm$estimate));}
	cat("\n\nPlotting hazard curves for",names[1],".");
	lazyhazplots(list(xd1,xd7,xd30),xhzpars,c(1,7,30),xlim=xlim);
	title(main=paste(names[1],"Hazard"));
	graphsave(path,names[1],".haz.eps",prompt);
	cat("\n\nPlotting hazard curves for",names[2],".");
	lazyhazplots(list(yd1,yd7,yd30),yhzpars,c(1,7,30),xlim=xlim);
	title(main=paste(names[2],"Hazard"));
	graphsave(path,names[2],".haz.eps",prompt);

	if(top){
	myzsc<-zsc[match("*",zsc$sig),c(1,8)];
	if(class(qreg.tab)!="character"){
		myqreg<-matrix(qreg.tab[qreg.tab[,5]<.05,c(1,5)],ncol=2);
	} else {myqreg<-cbind(NA,NA);}
	lz<-dim(myzsc)[1];lqr<-dim(myqreg)[1];
	report<-c(mean(x),round(quantile(x,c(.5,.9)),0),max(x),lr$stat,NA,xmod$g$estimate,NA,myqreg[,1],myzsc[,1]);
	report<-cbind(report,c(mean(y),round(quantile(y,c(.5,.9)),0),max(y),NA,NA,ymod$g$estimate,NA,rep(NA,lqr+lz)));
	rownames(report)<-c('mean','median','90%','max','log-rank',"",'initial mortality','mortality acceleration',
			    "quantile regression",rep("",lqr),paste("z-score",rownames(myzsc)));
	colnames(report)<-names;
	report<-cbind(report,p=
		c(tt$p.value,zsc[match(c(.5,.9),rownames(zsc)),8],rep(NA,3),xymod.tab[,5],NA,myqreg[,2],myzsc[,2]));
	}

	if(top){
		sig.tests<-tests[sig.tests]; 
		if(length(sig.tests)==0){sig.tests<-0;names(sig.tests)<-"None";}
		rname<-paste(paste(names,collapse=""),"surv",sep=".");
		fname<-paste(path,rname,".rdata",sep="");
		out<-list(smry=smry,zsc=zsc,x=x,y=y,xy=xy,names=names,
			group=group,lr=lr,tt=tt,tests=tests,sig.tests=sig.tests,
			qreg=qreg,qreg.tab=qreg.tab,qreg.sum=qreg.sum,
			xd1=xd1,xd7=xd7,xd30=xd30,
			yd1=yd1,yd7=yd7,yd30=yd30,xhzpars=xhzpars,yhzpars=yhzpars,
			xmod=xmod,ymod=ymod,modx=modx,mody=mody,
			modxy=modxy,xymod=xymod,report=report,
			xymod.tab=xymod.tab,xysurvfit=xysurvfit)
		class(out)<-"survomatic";
		go<-get(as.character(sys.call()[1]));
		out$go<-go; out$sys<-function(q){print(group);print(sys.status());print(sys.call());}
		if(save){
			assign(rname,out);
			formals(.First)$x<-rname;
			save(list=c(rname,"checkcons","checkpx","ctrl","demog","dex",
				    "exg","exgm","exl","exlm","ezz","fixpar","go","graphsave",
				    "grf","hsg","hsgm","hsl","hslm","instpx","lazyhazplots",
				    "lazymltplot","objf","opsurv","srvhaz","zpprob","zptest",
			            ".First",names),file=fname);
			cat("\nYour results have been saved as ",fname);
		}
	}

	cat("\nThe following test/s may have produced significant results:\n");
	print(names(sig.tests));
	cat("\nExecutive Summary:\n");
	print(report[1:4,],na.print="");
	print('\n-----------------------\n');
	print(report[5:dim(report)[1],],na.print="");
	cat("\n\nMethods:");
	cat("\nMeans of log-transformed survival times were compared using the Student\'s t-test. Quantiles, including the median, were compared using a modified version of the score test described in Wang et. al. (2004). The log-rank test (Fleming and Harrington, 1991) was used to compare distributions of survival times with the p-value adjusted based on 2000 permutations of the data (Andersen et. al., 1993). The mortality parameters were fitted and compared using the maximum likelihood method adapted from Pletcher et. al. (2000). The R statistical language (R Development Core Team, 2008) was used for all calculations except where specified. The R scripts used to perform these calculations are available from the authors upon request.\n");
	cat("\nThank you for using Survomatic! "); 
	if(top&save) cat("All the analysis you have done will now be saved for later use as an R file entitled",fname,".")
	cat("\n\n\n");
	options('warn'= 0);options(show.error.messages=T);
	if(top) invisible(out);
}

attr(go,'v')<-0.1;

graphsave<-function(path,names,suffix,prompt){
	if(prompt<2){return();}
	fname<-paste(path,paste(names,collapse=""),suffix,sep="");
	message<-paste("Do you wish to save this graph as",fname,"?");
	choices<-c("Save under suggested name.",
		   "Give a different name (point and click).",
		   "Type in a name from the console.",
		   "Don\'t save.")
	input<-menu(choices,graphics=T,title=message);
	switch(input,
		dev.copy2eps(file=fname),
		{fname<-tclvalue(tkgetSaveFile(initialdir=path,initialfile=paste(names,collapse=""),
			defaultextension=suffix));},
		{cat("\nPlease type the name and path to which you want to save the file.");
			fname<-scan("","character",nlines=1,quiet=T);print(fname);
			dev.copy2eps(file=fname);},
		cat("\nContinuing without saving.")
	)

}

demog<-function(t,int=1){
	if(int>1){t<-round(t/int);}
	t<-sort(t); t.rle<-rle(t); d<-px<-lx<-c();
	ts<-1:max(t.rle$values);
	for(i in ts){
		
		j<-match(i,t.rle$values);
		if(is.na(j)){d<-c(d,0);}else{d<-c(d,t.rle$lengths[j]);}
	}
	for(i in ts) {
		lx<-c(lx,(sum(d)-sum(d[1:i]))/sum(d));
	}
	for(i in ts) {px<-c(px,lx[i+1]/lx[i]);}
	ux<- -log(px); ux[is.infinite(ux)]<-NA; lnux<-log(ux); lnux[is.infinite(lnux)]<-NA;
	out<-cbind(ts,d,lx,px,ux,lnux);
	colnames(out)<-c("time","deaths","lx","px","ux","lnux");
	data.frame(out);
}

lazyhazplots<-function(ts,pars,ints,cols=rainbow(length(ts)),xlim=NULL,add=F){
	mnt<-mnx<-mx<-c(); 
	for(j in ts){
		mx<-c(mx,max(j$ux,na.rm=T)); mnx<-c(mnx,min(j[j$ux>0,]$ux,na.rm=T));}
	mnt<-match(T,ts[[which.max(ints)]]$ux>0)
	mx<-max(mx); mnx<-min(mnx); print(mnx); if(is.null(xlim)){xrng<-c(mnt,max(ts[[which.min(ints)]]$time));}else{
		xrng<-c(10,xlim);}
	print(xrng);
	if(add){points(1+ts[[1]]$time,ts[[1]]$ux,pch=3,col=cols[1]);} else {
		plot(1+ts[[1]]$time,ts[[1]]$ux,log='xy',type='p', pch=3,
			xlim=xrng,ylim=c(mnx,mx),xlab="Time",ylab="Hazard Rate",col=cols[1]);
	} 
	xrng<-xrng[1]:xrng[2];
	for(j in 2:length(ts)){
	points(1+ts[[j]]$time,ts[[j]]$ux,col=cols[j],pch=3); 
		lines(1+xrng,ints[j]*srvhaz(1+xrng,pars[1,1],pars[1,2],i=ints[j]),col=cols[j],lty=2);
		lines(1+xrng,ints[j]*srvhaz(1+xrng,pars[2,1],pars[2,2],pars[2,3],i=ints[j]),col=cols[j],lty=1);
		lines(1+xrng,srvhaz(1+xrng,pars[1,1],pars[1,2]),col=cols[1],lty=2);
		lines(1+xrng,srvhaz(1+xrng,pars[2,1],pars[2,2],pars[2,3]),col=cols[1],lty=1);
	}
}

srvhaz<-function(x,a,b,c=0,s=0,i=1){x<-x*i;c+a*exp(b*x)/(1+s*a*(exp(b*x)-1)/b);}

lazymltplot<-function(d,subset,xmax=1350){
	names<-lwds<-ltys<-comps<-groups<-exps<-c(); exp<-0;
	for(i in subset){
		inames<-d[[i]]$names; 
		m1<-match(inames[1],names); m2<-match(inames[2],names);
		if(is.na(m1)){
			names<-c(names,inames[1]);
			lwds<-c(lwds,1);
			ltys<-c(ltys,2);
			comps<-c(comps,i);
			groups<-c(groups,2);
			exp<-exp+1; exps<-c(exps,exp);
		}
		if(is.na(m2)){
			names<-c(names,inames[2]);
			ltys<-c(ltys,1);
			if(min(d[[i]]$sig.tests)>0){
				if(is.na(match("lr",d[[i]]$sig.tests))){
					lwds<-c(lwds,2);
				} else { lwds<-c(lwds,4);}
			} else { lwds<-c(lwds,1);}
			comps<-c(comps,i);
			groups<-c(groups,1);
			exps<-c(exps,exp);
		}
	}
	plt<-data.frame(names,lwds,ltys,comps,groups,exps);
	nexp<-length(unique(plt$exps));
	cols<-c(); 
	for(i in 1:nexp){
		li<-sum(plt$exps==i);
		cols<-c(cols,hsv(i*seq(.8,1,len=li)/nexp,v=seq(1,.5,len=li)));
	}
	#rainbow(sum(plt$exps==i),start=(i-.001)/nexp,end=i/nexp));}
	plot(d[[plt$comps[1]]]$xysurvfit[plt$groups[1]],lwd=plt$lwds[1],lty=plt$ltys[1],col=cols[1],xlim=c(0,xmax),bg='gray75');
	for(i in 2:dim(plt)[1]){
		lines(d[[plt[i,]$comps]]$xysurvfit[plt[i,]$groups],lwd=plt[i,]$lwds,
		lty=plt[i,]$ltys,col=cols[i]);
	}
	legend("bottomleft",legend=plt$names,lwd=plt$lwds,lty=plt$ltys,col=cols)
	graphsave('~/temp/aging08/',c('female','survival'),'.eps',prompt=2);
	lazyhazplots(list(d[[1]]$xd1,d[[1]]$xd7,d[[1]]$xd30),d[[1]]$xhzpars,c(1,7,30),xlim=xmax,cols=rep(cols[1],3));
	for(i in 2:dim(plt)[1]){
		switch(plt[i,]$groups,
			{d1<-d[[i]]$yd1;d7<-d[[i]]$yd7;d30<-d[[i]]$yd30;pars<-d[[i]]$yhzpars;},
			{d1<-d[[i]]$xd1;d7<-d[[i]]$xd7;d30<-d[[i]]$xd30;pars<-d[[i]]$xhzpars;});
		lazyhazplots(list(d1,d7,d30),pars,c(1,7,30),xlim=xmax,cols=rep(cols[i],3),add=T);
	}
	legend("topleft",legend=plt$names,col=cols,pch=3)
	print(plt); print(cols);
}

########################################################################################

# Alex F. Bokov, GPL 2008
# Based on the work of 
# Chenxi Wang, Qing Li, David T. Redden, Richard Weindruch, and David B. Allison
# as published in Mechanisms of Ageing and Development 125 (2004) pp.629-632
# with some modifications by Alex F. Bokov
# But at present time not vetted or endorsed by the above, or by any professional
# statistician. Use at own risk!
# For technical questions or bug reports, please contact bokov@uthscsa.edu

# Usage:
# ezz() is the main function; the others are called by it
# d = vector of survival times (required)
# g = vector of group assignments (required)
# The following have default values and are optional:
# quant = the quantile to be tested, or a vector of quantiles
# step = governs granularity when finding the supermum
# thresh = significance cutoff
# gnames = names of the groups being compared
# noob = gives abbreviated output if set to true

# The script outputs a table (data frame) with the following data:
#  quantiles,
#  pooled age at death for that quantile, 
#  number of survivors from group 1 at that quantile
#  fraction of survivors from group 1 at that quantile
#  number of survivors from group 2 at that quantile
#  fraction of survivors from group 2 at that quantile
#  total fraction of survivors at that quantile
#  the Zp test statistic
#  the corresponding p value
#  an asterisk indicating significance relative to the chosen thresshold

# This script assumes that only two groups are being compared.
# Not designed/tested for multiple groups.

ezz<-function(d,d2=NULL,g,quant=.9,step=.01,thresh=.05,
	      gnames=NULL,noob=F,thetas=c(.05,.95),debug=F){
	if(is.null(d2)){
		groups<-levels(as.factor(g));
		d1<-subset(d,g==groups[1]);
		d2<-subset(d,g==groups[2]);
	} else { d1<-d; if(is.null(g)){groups<-c('a','b');}else{groups<-g;}}
	#mintheta<-min(1/d1,1/d2);
	n1<-length(d1); n2<-length(d2);
	n<-n1+n2;
	output<-c(); q<-quantile(d,quant);
	if(debug) cat("\nRunning zptest on x1 and x2.");
	for(i in q){
		x1<-sum(d1>i); x2<-sum(d2>i);
		output<-rbind(output,zptest(x1,x2,n1,n2,n));
		if(debug) cat('.');
	}
	output<-cbind(q,output);
	if(is.null(gnames)){g1<-groups[1];g2<-groups[2];}
	else{g1<-gnames[1];g2<-gnames[2];}
	surv1<-paste("srv",g1,sep='.');
	surv2<-paste("srv",g2,sep='.');
	surv1n1<-paste("srv",g1,"frc",sep='.');
	surv2n2<-paste("srv",g2,"frc",sep='.');
	colnames(output)<-c('age',surv1,surv1n1,surv2,surv2n2,'survtotal.frc','zp');
	output<-data.frame(output);
	output$zp[is.nan(output$zp)]<-9999; #print(output$zp);

	p<-zptests<-c();
	if(debug) cat("\nRunning zptest on k n1\'s and n2\'s.");
	for(k in 0:n1){zptests<-c(zptests,zptest(k,0:n2,n1,n2,n)[,6]);if(debug) cat('.');}
	zptests[is.nan(zptests)]<-9999;

	for(i in 1:dim(output)[1]){
		sup<-c();
		for(j in seq(thetas[1],thetas[2],step)){
			zptemp<-c();
			if(debug) cat("\nRunning zpprob on n1\'s and n2\'s.");
			for(k in 0:n1){zptemp<-c(zptemp,zpprob(k,0:n2,n1,n2,j));
				if(debug) cat('.');}
			keep<-abs(zptests)>abs(output$zp[i]);
			sup<-c(sup,sum(zptemp[keep]));
		}
		p<-c(p,max(sup));
	}
	sig<-p<thresh; sig[sig==T]<-"*"; sig[sig==F]<-"";
	output<-cbind(output,p,sig);
	rownames(output)<-quant;
	if(noob){output<-output[,c(1:5,8,9)];}
	output;
}

zptest<-function(x1,x2,n1,n2,n){
	that1<-x1/n1; that2<-x2/n2;
	ttil<-(x1+x2)/n;
	denom<-sqrt(ttil*(1 - ttil)*(1/n1 + 1/n2));
	zp<-(that1 - that2)/denom;
	output<-cbind(x1,that1,x2,that2,ttil,zp);
	output;
}

zpprob<-function(x1,x2,n1,n2,th){
	choose(n1,x1)*choose(n2,x2)*th^(x1+x2)*(1-th)^(n1+n2-x1-x2);
}

########################################################################################

# Alex F. Bokov, GPL 2008
# Based on the work of 
# Scott D. Pletcher, Aziz A. Khazaeli, and James W. Curtsinger
# as published in Journal of Gerontology: Biological Sciences 55A (2000) pp.B381-B389
# with some modifications by Alex F. Bokov
# But at present time not vetted or endorsed by the above, or by any professional
# statistician. Use at own risk!
# For technical questions or bug reports, please contact bokov@uthscsa.edu

# Usage:


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


.First();