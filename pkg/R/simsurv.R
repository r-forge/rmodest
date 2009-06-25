`simsurv`<-
function(n,type='g',p=c(2.433083e-05,0.005,3e-11,0.0015)){
	model<-switch(type,w=1,l=2,lm=2,g=3,gm=3,l2=4,e=5);
	if((model==2|model==4)&p[4]<6e-10){model=3;}
	if((model==3|model==1)&p[2]<2e-19){model=5;}
	if(model==5){out<-as.integer(rexp(n,as.double(p[1])));} else {
		out<-
		.C('simsurv',a=as.double(p[1]),b=as.double(p[2]),c=as.double(p[3]),s=as.double(p[4]),
		   size=as.integer(n),model=as.integer(model),ans=double(n),PACKAGE="Survomatic")$ans
	}
	attr(out,'type')<-type;attr(out,'pars')<-p;
	return(out);
}