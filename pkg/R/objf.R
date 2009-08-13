`objf` <-
function(par,lb=0,ub=10,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,ex="ofg",nx,ny,cx,cy,tlog=F){
	#cat('cx:',cx,'\n','cy:',cy,'\n');
	#if(ex=='ofw'&mean(cons)!=1){browser();}
	if(tlog){par<-exp(par);}
	if(max(par<lb)!=0|max(par>=ub)!=0){return(NA);}
	par1<-par2<-rep.int(0,4); par1[keep]<-par2[keep]<-par[1:np];
	sum1<-try(.C(ex,a=as.double(par1[1]),b=as.double(par1[2]),c=as.double(par1[3]),
		     s=as.double(par1[4]),x=as.integer(x),size=as.integer(nx),censor=as.integer(cx),
		     ans=double(1),PACKAGE="Survomatic"));
	if(class(sum1)[1]!="try-error"){sum1<-sum1$ans;}else{return(NA);}
	if(!is.null(y)){
		par2[keep][cons[keep]>0]<-par[(np+1):length(par)];
		sum2<-try(.C(ex,a=as.double(par2[1]),b=as.double(par2[2]),c=as.double(par2[3]),
			     s=as.double(par2[4]),x=as.integer(y),size=as.integer(ny),
			     censor=as.integer(cy),ans=double(1),PACKAGE="Survomatic"));
		if(class(sum2)[1]!="try-error"){sum2<-sum2$ans;}else{return(NA);}
	} else { sum2 <- 0; }
	out<-sum1+sum2;
	if(!is.na(out) & is.finite(out) & out < -50){return(out);}else{return(NA);}
}