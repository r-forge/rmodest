dev.copy2emf<-function(file=Rplot,...){
	if(Sys.info()[['sysname']]=='Linux' & 
	   length(system('which fig2dev',intern=T))>0){
		dev.print(device=xfig,file=paste(file,'.fig',sep=''),...);
		system(paste('fig2dev -L emf ',file,'.fig ',
			     '> ',file,'.emf',sep=''));
	} else stop('This command only works for Linux with the transfig package installed. If you are running Linux, please make sure that the \'fig2dev\' command exists and is in your path, then try again.');
}
