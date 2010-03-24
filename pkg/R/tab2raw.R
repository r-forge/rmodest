tab2raw<-function(d,output='times',timecol=1,eventcol=2,groupcol=3,censcol=4,choosegroup=0){
	if(d[nrow(d),1]== -1) {d<-d[-nrow(d),];}
	if(choosegroup){d<-d[d[,groupcol]==choosegroup,];}
	cens=F;time=F;
	if(match('times',output,nomatch=F)){tout<-rep(d[,timecol],d[,eventcol]); time=T;}
	if(match('censor',output,nomatch=F)){cout<-rep(d[,censcol],d[,eventcol]); cens=T;}
	if(time & cens){return(data.frame(time=tout,censor=cout));}
	if(time){return(tout);}
	if(cens){return(cout);}
}