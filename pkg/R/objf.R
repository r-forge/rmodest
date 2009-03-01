# This is the old, interpreted version, which will be removed once the compiled version is verified
`objf2` <-
function(par,lb=0,ub=1,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,ex=exg){
	if(max(par<=lb)!=0|max(par>=ub)!=0){return(NA);}
	par1<-par2<-rep.int(NA,4); par1[keep]<-par2[keep]<-par[1:np];
	par2[keep][cons[keep]>0]<-par[(np+1):length(par)];
	env1<-list(a=par1[1],b=par1[2],c=par1[3],s=par1[4],x=x); 
	env2<-list(a=par2[1],b=par2[2],c=par2[3],s=par2[4],x=y);
	out<-sum(eval(ex,envir=env1))+sum(eval(ex,envir=env2));
	if(!is.na(out) & is.finite(out) & out < -50){return(out);}else{return(NA);}
}

`objf` <-
function(par,lb=0,ub=1,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,ex="ofg",nx,ny,cx,cy){
	if(max(par<=lb)!=0|max(par>=ub)!=0){return(NA);}
	par1<-par2<-rep.int(0,4); par1[keep]<-par2[keep]<-par[1:np];
	sum1<-.C(ex,a=as.double(par1[1]),b=as.double(par1[2]),c=as.double(par1[3]),s=as.double(par1[4]),x=as.integer(x),size=as.integer(nx),censor=as.integer(cx),ans=double(1),PACKAGE="Survomatic")$ans;
	if(!is.null(y)){
		par2[keep][cons[keep]>0]<-par[(np+1):length(par)];
		sum2<-.C(ex,a=as.double(par2[1]),b=as.double(par2[2]),c=as.double(par2[3]),s=as.double(par2[4]),x=as.integer(y),size=as.integer(ny),censor=as.integer(cy),ans=double(1),PACKAGE="Survomatic")$ans;
	} else { sum2 <- 0; }
	out<-sum1+sum2;
	if(!is.na(out) & is.finite(out) & out < -50){return(out);}else{return(NA);}
}
