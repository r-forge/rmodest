`srvshp`<-
function(x,a,b,c=0,s=0,i=1,model='g'){
	model<-switch(model,w=1,l=2,g=3,stop('Model must be one of \'w\', \'l\', or \'g\'.'));
	n<-length(x);
	out <- .C('srvshp',a=as.double(a),b=as.double(b),c=as.double(c),s=as.double(s),
	   size=as.integer(n),model=as.integer(model),x=as.integer(x),ans=double(n),
	   PACKAGE="Survomatic")$ans;
	return(out);
}
