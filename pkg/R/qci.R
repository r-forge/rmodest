`qci` <-
function(x,p=.5,cn=.95){
	x<-sort(x);
	cn<-1.01+cn; n<-length(x); np<-n*p; er<-cn*sqrt(np*(1-p));
	lb<-ceiling(np-er); ub<-ceiling(np+er);
	lb[lb<1]<-1; ub[ub>n]<-n;
	#print(np); print(er); print(lb); print(ub);
	cbind(x[lb],x[ub]);
}