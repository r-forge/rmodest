`grf` <-
function(par,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,model='g'){
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

