`lazymltplot` <-
function(d,subset,xmax=1350){
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

