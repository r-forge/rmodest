`objf` <-
function(par,lb=0,ub=1,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,ex=exg){
	if(max(par<=lb)!=0|max(par>=ub)!=0){return(NA);}
	par1<-par2<-rep.int(NA,4); par1[keep]<-par2[keep]<-par[1:np];
	par2[keep][cons[keep]>0]<-par[(np+1):length(par)];
	env1<-list(a=par1[1],b=par1[2],c=par1[3],s=par1[4],x=x); 
	env2<-list(a=par2[1],b=par2[2],c=par2[3],s=par2[4],x=y);
	out<-sum(eval(ex,envir=env1))+sum(eval(ex,envir=env2));
	if(!is.na(out) & is.finite(out) & out < -50){return(out);}else{return(NA);}
}

