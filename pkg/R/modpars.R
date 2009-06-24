`modpars`<-
function(x,modeli,modelo=NULL,cni=rep(T,4),cno=NULL,nil=1e-7,pf=mean,trim=F,onegrp=F){
	if(onegrp){x<-c(x,x); cni<-rep(T,4);}
	if(sum(is.null(x))>0){print('input params contain nulls!');traceback();browser();};
	if(is.null(modelo)){modelo<-modeli;}; if(is.null(cno)){cno<-cni;};
	keepi<-switch(modeli,e=c(1,5),w=c(1,2,5,6),g=c(1,2,5,6),gm=c(1:3,5:7),l=c(1,2,4:6,8),lm=1:8);
	keepo<-switch(modelo,e=c(1,5),w=c(1,2,5,6),g=c(1,2,5,6),gm=c(1:3,5:7),l=c(1,2,4:6,8),lm=1:8);
	# Will eventually do more input validation, but for now let's agree to pass 
	# constraints as vectors of either 4 or 1 logical values
	cni <- !cni; if(length(cni)==1){cni<-rep(cni,4);};
	cno <- !cno; if(length(cno)==1){cno<-rep(cno,4);};
	if(length(cni)+length(cno)!=8){stop('In this version constraints must be specified
either as a single boolean value or a vector
of four boolean values. Please fix your input
and try again.');}
	keepi<-keepi[c(rep(T,4),!cni)[keepi]];
	if(length(x)!=length(keepi)){
		print('Mismatch between model and parameter length.');traceback();browser();
	}
	# Here we convert the unique-values-only parameter vector into standard form
	out<-rep.int(NA,8); out[keepi]<-x; 
	if(sum(cni)>0){out[5:8][cni]<-out[1:4][cni];}
	# We insure that any parameters used by the target model are populated
	out[keepo][is.na(out[keepo])]<-nil;
	# Now we apply the target model's constraints
	if(sum(cno)>0){
		# If the target constraint is different from the input constraint,
		# average the two parameters using the function specified by pf
		# No point in doing this is the constraints are the same
		if(sum(cni==cno)<4&sum(out[cno]>0)){
			out[cno]<-apply(out*as.matrix(diag(cno)[c(1:4,1:4),cno]),2,
				       function(x,f){f(na.omit(x[x!=0]));},pf);
		}
		if(trim){out[5:8][cno]<-NA;}
	}
	# Gotta remember to trash the unused parameters for the target model.
	out[-keepo]<-NA;
	if(trim){out<-as.numeric(na.omit(out));} 
	return(out);
}
