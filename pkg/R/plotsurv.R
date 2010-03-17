`plotsurv`<-
function(x,y,cx=NULL,cy=NULL,col=c(1,'darkred'),lwd=2,legend=NULL,lloc='bottomleft',pcex=1,pch=24:25,ylim=0:1,bg=NULL,...){
	lx<-length(x);ly<-length(y);args<-list(...);
	if(is.null(cx)){cx<-rep(1,lx);}; if(is.null(cy)){cy<-rep(1,ly);};
	if(is.null(bg)){bg<-col;}
	sf<-survfit(Surv(c(x,y),event=c(cx,cy))~c(rep(0,lx),rep(1,ly)));
	plot(sf,col=col,lwd=lwd,ylim=ylim,...);
	# not a great way to do it, but will do for now
	if(!is.null(args$fun)){
		if(args$fun=='event')for(i in 1:2)
			{points(sf[i]$time,1-sf[i]$surv,pch=pch[i],cex=pcex,lwd=lwd,bg=bg[i],col=col[i])}
	} else {for(i in 1:2){points(sf[i],pch=pch[i],cex=pcex,lwd=lwd,bg=bg[i],col=col[i]);}}
	if(!is.null(legend)){legend(x=lloc,col=col,pt.bg=bg,lwd=lwd,pch=pch,legend=legend);}
	invisible(sf);
}
