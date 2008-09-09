`go` <-
function(x,y,xynames=c(),units="d",save=T,path,prompt=2,xlim=1350){
	options(warn=-1);options(show.error.messages=F);
	quantiles<-list(seq(.01,.99,.01),seq(.05,.95,.05),seq(.1,.9,.1),seq(.2,1,.2));
	mycall<-strsplit(as.character(sys.call()[1]),"\\$")[[1]];
	if(length(mycall)>1){
		top=F;d<-get(mycall[1]);
		if(length(mycall)>2){
			for(i in 2:(length(mycall)-1)){
				d<-d[[mycall[i]]];
		}}
		x<-d$x;y<-d$y;xynames<-d$xynames;smry<-d$smry;zsc<-d$zsc;xy<-d$xy;
		group<-d$group;lr<-d$lr;qreg<-d$qreg;qreg.sum<-d$qreg.sum;tt<-d$tt;
		qreg.tab<-d$qreg.tab;xmod<-d$xmod;ymod<-d$ymod;xymod<-d$xymod;
		xd1<-d$xd1;xd7<-d$xd7;xd30<-d$xd30;yd1<-d$yd1;yd7<-d$yd7;yd30<-d$yd30;
		xhzpars<-d$xhzpars;yhzpars<-d$yhzpars;xymod.tab<-d$xymod.tab;
		modx<-d$modx;mody<-d$mody;modxy<-d$modxy;xysurvfit<-d$xysurvfit;
		sig.tests<-d$sig.tests;tests<-d$tests; report<-d$report; path<-d$path;
		sigqreg<-d$sigqreg;sigzsc<-d$sigzsc;
	}else{top=T;}
	if(missing(path)&top&save){
		message<-"Where would you like to save your results?\n";
		choices<-c("Let me choose using the point and click method.",
			   "Let me type in the location from the console.",
			   "I don\'t want to save any files. Run analysis without saving.");
		cat(message);
		input<-menu(choices,graphics=T,title=message);
		switch(input,
			{path<-paste(tclvalue(tkchooseDirectory()),"/",sep="");},
			{path<-readline("Please type in the location where you would like to save files: ");
			   path<-paste(path,"/",sep="");},
			{save<-F;cat('\nContinuing without saving.');}
		)
		if(save){cat('\nResults will be saved to:',path,'\n');}
	} else {path=NULL;}
	if(missing(x)&top){
		cat("\nYou need a group to compare. Please copy and paste one column of survival times from a spreadsheet.
They don\'t have to be in order, but make sure there are no empty spaces. Hit the enter key when you 
are done.\n");
		x<-scan(); 
	}
	if(length(xynames)<2&top){xynames[1]<-readline("Please give a brief name for the first group: ");}
	if(missing(y)&top){
		cat("\nYou need a second group to compare. Please copy and paste one column of survival times from a 
spreadsheet. They don\'t have to be in order, but make sure there are no empty spaces. Hit 
the enter key when you are done.\n");
		y<-scan(); 
	}
	if(length(xynames)<2&top){xynames[2]<-readline("Please give a brief name for the second group: ");}
	if(save&top){
		# ought to globally replace spaces with dots here
		assign(xynames[1],x,envir=.GlobalEnv);
		assign(xynames[2],y,envir=.GlobalEnv);
		cat("The survival times for the two groups have been saved as",xynames[1],"and",xynames[2],"inside R 
so you won\'t have to paste them in next time. Instead you can type 
\'go(",xynames[1],",",xynames[2],")\'\n");
	}

	if(top){
		tests<-c("t test"="tt","z score"="zsc","log rank"="lr","quantile regression"="qreg.tab",
			"mortality parameters"="xymod.tab");
		sig.tests<-c();
	}

	cat("\nDoing summary statistics...");
	if(top){smry<-t(cbind(summary(x),summary(y))); rownames(smry)<-xynames;}
	cat('\n'); print(smry);
	cat("\nDoing t-test on log-transformed data...");
	if(top){xt<-x[which(x>0)]; yt<-y[which(y>0)];tt<-t.test(log(xt),log(yt));if(tt$p.value<=.05) sig.tests<-c(sig.tests,1);}
	cat('\n'); print(tt);
	cat("\nDoing z-score test to see which quantiles differ..."); 
	if(top){
		zsc<-c();
		zsc<-ezz(x,y,xynames,quant=seq(0,.95,.05));
		if(length(levels(zsc$sig))>1) sig.tests<-c(sig.tests,2);
	} 
	cat('\n'); print(zsc);
	cat("\nDoing log-rank test with permutations...");
	if(top){
		xy<-c(x,y); group<-c(rep(xynames[1],length(x)),rep(xynames[2],length(y)));
		#print(group);
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
			qreg<-rq(log(xy)~relevel(factor(group),2),tau=quantiles[[i]]);i<-i+1;
			qreg.sum<-try(summary(qreg,se='ker')); cat('.');
		}
		if(class(qreg.sum)!="character"){
			qreg.tab<-c(); qreg.sig=F;
			for(i in qreg.sum){
				qreg.tab<-rbind(qreg.tab,c(quantile=i$tau,i$coefficients[2,]));
				if(min(qreg.tab[,5])<=.05){qreg.sig=T;}
			}
			if(qreg.sig){sig.tests<-c(sig.tests,4);}
		}else{qreg.tab<-qreg.sum;}
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
	cat("\nThe best model for",xynames[1],"is",names(xmod)[modx],".");
	cat("\nThe best model for",xynames[2],"is",names(ymod)[mody],".");
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

	if(top){
	sigzsc<-zsc[which(zsc$sig=="*"),c(2,4,8)];
	if(class(qreg.tab)!="character"){
		sigqreg<-matrix(qreg.tab[qreg.tab[,5]<.05,c(1:2,5)],ncol=3);
	} else {sigqreg<-matrix(NA,ncol=3);}
	if(dim(sigqreg)[1]==0){sigqreg<-matrix(NA,ncol=3);}
	lz<-dim(sigzsc)[1];lqr<-dim(sigqreg)[1];
	report<-c(length(x),mean(x),round(quantile(x,c(.5,.9)),0),max(x),lr$stat,NA,xmod$g$estimate,
		NA,if(match(NA,sigqreg[,1],no=F)){NA}else{quantile(x,sigqreg[,1])},sigzsc[,1]);
	report<-cbind(report,c(length(y),mean(y),round(quantile(y,c(.5,.9)),0),max(y),NA,NA,ymod$g$estimate,
		NA,if(match(NA,sigqreg[,1],no=F)){NA}else{quantile(y,sigqreg[,1])},sigzsc[,2]));
	rownames(report)<-c('n','mean','median','90%','max','log-rank',"",'initial mortality','mortality acceleration',"",if(match(NA,sigqreg[,1],no=F)){NA}else{paste("(qr) longevity of quantile",sigqreg[,1])},paste("(z score) number alive at quantile",rownames(sigzsc)))[1:dim(report)[1]];
	colnames(report)<-xynames;
	report<-cbind(report,p=
	   c(NA,tt$p.value,zsc[match(c(.5,.9),rownames(zsc)),8],NA,lr$pval,NA,xymod.tab[,5],NA,sigqreg[,3],sigzsc[,3]));
	}

	if(top){xysurvfit<-survfit(Surv(xy)~group);}
	tmain<-paste(xynames,collapse=" vs. ");
	cat("\nPlotting survival curves.\n");
	plot(xysurvfit,lwd=1,bty="n",xlim=c(0,xlim));
	for(i in 1:2){points(xysurvfit[i],pch=23+i,cex=1.5,lwd=1.5,bg=c('gray20','white')[i]);}
	title(main=paste(tmain,"Survival"));
	legend("topright",legend=xynames[2:1],pch=24:25,y.intersp=1.5,
		#inset=-0.006*psize,-+
		bty="n",pt.cex=1.5,cex=1.5,lwd=1.5,pt.bg=c('gray20','white'));
	
	if(save){graphsave(path,xynames,".srv.eps",prompt);}else{readline("Press enter to continue.")}
	if(class(qreg.sum)!="character"){
		cat("\n\nPlotting quantile regression curves.\n");
		plot(qreg.sum,parm=2,ols=F,main=paste(tmain,"Quantile Regression"));
		points(sigqreg[,1],sigqreg[,2],col='red',lwd=3);
		if(save){graphsave(path,xynames,".qr.eps",prompt);}else{readline("Press enter to continue.")}
	}

	if(top){
		xhzpars<-rbind(c(xmod$g$estimate,NA),c(xmod$gm$estimate));
		yhzpars<-rbind(c(xmod$g$estimate,NA),c(xmod$gm$estimate));}
	cat("\n\nPlotting hazard curves.\n");
	#lazyhazplots(list(xd1,xd7,xd30),xhzpars,c(1,7,30),xlim=xlim);
	lazyhazplots(list(xd30),xhzpars,30,xlim=xlim/30,cols='black',main=paste(tmain,"Hazard"));
	lazyhazplots(list(yd30),yhzpars,30,xlim=xlim/30,cols='blue',add=T);
	if(save){graphsave(path,xynames[2],".hz.eps",prompt);}else{readline("Press enter to continue.")}

	if(top){
		sig.tests<-tests[sig.tests]; 
		if(length(sig.tests)==0){sig.tests<-0;names(sig.tests)<-"None";}
		out<-list(smry=smry,zsc=zsc,x=x,y=y,xy=xy,xynames=xynames,path=path,
			group=group,lr=lr,tt=tt,tests=tests,sig.tests=sig.tests,
			qreg=qreg,qreg.tab=qreg.tab,qreg.sum=qreg.sum,
			xd1=xd1,xd7=xd7,xd30=xd30,
			yd1=yd1,yd7=yd7,yd30=yd30,xhzpars=xhzpars,yhzpars=yhzpars,
			xmod=xmod,ymod=ymod,modx=modx,mody=mody,
			modxy=modxy,xymod=xymod,report=report,sigqreg=sigqreg,
			sigzsc=sigzsc,xymod.tab=xymod.tab,xysurvfit=xysurvfit)
		class(out)<-"survomatic";
		go<-get(as.character(sys.call()[1]));
		out$go<-go; out$sys<-function(q){print(group);print(sys.status());print(sys.call());}
		if(save){
			rname<-paste(paste(xynames,collapse=""),"surv",sep=".");
			fname<-paste(path,rname,".rdata",sep="");
			assign(rname,out);
			formals(.First)$x<-rname;
			save(list=c(rname),file=fname);
			cat("\nYour results have been saved as ",fname);
		}
	}

	cat("\nThe following test/s may have produced significant results:\n");
	print(names(sig.tests));
	cat("\nExecutive Summary:\n");
	print(report[1:5,],na.print="");
	cat('\n-----------------------\nLog-rank: p =',report[6,3],'\n');
	cat('\n-----------------------\nMortality\n');
	print(report[8:9,],na.print="");
	cat('\n-----------------------\nQuantile Regression & Z-score\n');
	print(report[11:dim(report)[1],],na.print="");
	cat("\n\nMethods:");
	cat("\nMeans of log-transformed survival times were compared using the Student\'s 
t-test. Quantiles, including the median, were compared using a modified 
version of the score test described in Wang et. al. (2004). Quantile 
regression was done as described by Koenker and Geling (2001). The 
log-rank test (Fleming and Harrington, 1991) was used to compare 
distributions of survival times with the p-value adjusted based on 2000 
permutations of the data (Andersen et. al., 1993). The mortality 
parameters were fitted and compared using the maximum likelihood method 
adapted from Pletcher et. al. (2000). The R statistical language (R 
Development Core Team, 2008) was used for all calculations except where 
specified. The R scripts used to perform these calculations are 
available from the authors upon request.\n");
	cat("\nThank you for using Survomatic! "); 
	if(top&save) cat("All the analysis you have done will now be 
saved for later use as an R file entitled
",fname,".")
	cat("\n\n\n");
	options('warn'= 0);options(show.error.messages=T);
	if(top) invisible(out);
}

