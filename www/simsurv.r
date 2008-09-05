# Arguments to simsurv explained:
# n:		sample size; can be a vector; if smaller than g, the last 
#		size gets repeated multiple times
# g:		number of populations
# type:		model ('g' or 'l'); can be a vector; if smaller than g, the #		last type gets repeated multiple times
# seed:		if a duplicate set of simulations is desired for some reason, 
#		the seed from that simulation can be reused and the same 
#		exact values will be generated. Reproduceability not tested 
#		accross multiple hardware/software platforms, though, so use 
#		at own risk.
# names:	vector of group names; if it runs out of names it will just 
#		name successive groups NULL
# tests:	whether to do summary statistics, log-ranks, etc. on the 
#		groups; will NOT work in this version
# surv:		whether to output a surv-list object; will NOT work in this 
#		version
# fp:		path to save temporary output files, in case process gets 
#		interrupted
# minpars:	floor on parameters; can also be a singular value
# maxpars:	ceiling on parameters; can also be a singular value
# minmed:	populations with a hazard rate >= 0.5 occuring before this
#		timepoint will be rejected
# maxmed:	populations with a hazard rate >= 0.5 occuring after this
#		timepoint will be rejected
# pars:		if not NULL, the exact pars specified will be used instead
#		of randomly generated; can be a matrix, and if the number
#		of rows is less than g, the last row will be repeated; if
#		you use this parameter, the script assumes you know what
#		you are doing and disables the minmed/maxmed sanity 
#		checking
# prevsurv:	the output from a previous run of simsurv; this uses
#		some of the settings from that run, overriding the n, g, 
#		pars, type, and names parameters (thereby also disabling 
#		sanity checking).

simsurv<-function(n=100,g=2,type='g',seed=runif(1,0,1e9),names=letters[1:g],
	tests=F,surv=F,fp='~/temp/',minpars=c(-19,-12,-17,-3),
	maxpars=c(-3,-4,-3,-.1),minmed=20,maxmed=6e4,pars=NULL,
	prevsurv=NULL){
	# Start the timer, for measuring runtime.
	tstart<-as.numeric(Sys.time());
	# Override the parameters if prevsurv is present
	if(!is.null(prevsurv)){
		n<-prevsurv$n; g<-prevsurv$g; pars<-prevsurv$pars;
		type<-prevsurv$type; names<-names(prevsurv$r);
	}
	# Check the lengths of various arguments and fill them out 
	# as necessary.
	ln<-length(n); lt<-length(type);
	if(ln<g){n<-c(n,rep(n[ln],g-ln));}
	if(lt<g){type<-c(type,rep(type[lt],g-lt));}
	if(!is.null(pars)){
		pars<-rbind(pars); pt<-pars[dim(pars)[1],];
		while(dim(pars)[1]<g){pars<-rbind(pars,pt);}
	}
	# Initialize some variables
	out<-list(); out$r<-as.list(1:g); simlog<<-c();
	set.seed(seed); 
	# Here can be seen the structure of the output object.
	# The actual survival times will be in list $r
	out$call<-sys.call(); out$seed<-seed; out$n<-n; out$g<-g; 
	out$pars<-pars; out$type<-type; out$risk<-c();
	# Main loop.
	for(i in 1:g){
		cat("\n-------\nLine",i);
		# Initialize some loop variables.
		group<-c(); risk<-c(); imed<-0;
		# Pick the proper survivorship function to use 
		sf<-switch(type[i],g=survfng,l=survfnl);
		if(is.null(sf)){stop('For simsurv, type must be \'g\', \'l\', or left at default.');}
		# If pars are specified, use them.
		if(!is.null(pars)){
			ai<-pars[i,1]; bi<-pars[i,2]; 
			ci<-pars[i,3]; si<-pars[i,4];}
		cat("\nCalculating params & risks. Median:");
		# Keep generating parameters until we get a valid
		# 50% hazard rate between minmed and maxmed (unless 
		# pars have been specified).
		while(is.na(imed)|imed<minmed|imed>maxmed){
			if(is.null(pars)){
				# Set a and b
				ai<-exp(runif(1,minpars[1],maxpars[1]));
				bi<-exp(runif(1,minpars[2],maxpars[2]));
				# Half the time, c will be zero
				ci<-exp(runif(1,minpars[3],maxpars[3]))*round(runif(1,0,1));
				# s gets set only if the model for this 
				# population is a logistic one, otherwise
				# it's 0
				if(type[i]=='l'){si<-exp(runif(1,minpars[4],maxpars[4]));} else si<-0;
			# If parameters have been specified, the above does
			# not happen and instead this loop exits after
			# the first iteration
			} else imed<-minmed+1;
			# tt is an index variable, hardcoded to a million
			# places, which is probably overkill.
			tt<-1:1e5;
			# Initialize the hazard rate with the above 
			# variables and drop any NA values that might
			# be generated at very late timepoints
			risk<-1-sf(tt+1,ai,bi,ci,si)/sf(tt,ai,bi,ci,si);
			risk<-risk[!is.na(risk)];
			# Record the location of the 0.5 hazard rate to
			# determine whether it's necessary to continue
			# the while loop.
			imed<-match(T,0.5<=risk); cat(" ",imed);
		}
		# Now we have a useable time-ordered list of hazard rates
		# If the pars have been generated from scratch, save them
		# to the output object. 
		if(is.null(pars)){out$pars<-rbind(out$pars,c(ai,bi,ci,si));}
		# We also save the hazard rates; thanks for the idea, Jon.
		out$risk<-rbind(out$risk,risk);
		cat("\nModel:",type[i],"Pars:",ai,bi,ci,si,"\nProgress: ");
		# For each individual j out of the n individuals in 
		# population i, roll a vector of random numbers in (0,1)
		# equal in length to risk. Record the first such number 
		# that is larger than the corresponding element of the
		# risk vector and add it onto the growing list of ages
		# at death (named group).
		for(j in 1:n[i]){
			t<-match(T,runif(length(risk))<risk);
			cat(' ',j,":",t);
			group<-c(group,t);
		}
		# Assign group names
		out$r[[i]]<-group; names(out$r)[i]<-names[i];
		# Write a temp output object to the parent environment,
		# and then write it out to a file. To save soul-crushing
		# anguish in case the script or the session has to be
		# interrupted.
		simsurv.temp<<-out;
		save(simsurv.temp,file=paste(fp,'simsurv.',tstart,'.rdata',sep=''));
	}
	# Not implemented here, and very slow when it used to be.
	if(surv){out<-raw2surv(out$r,names,tests=tests);}
	# Delete the temp file
	simsurv.temp<<-NULL;
	# Calculate and report the runtime.
	cat("\n"); ndp<-sum(n*g); 
	trun<-as.numeric(Sys.time())-tstart;
	tdp<-trun/ndp;
	cat("Simulated a total of",ndp,"datapoints in",trun,"seconds,",tdp,"seconds per datapoint.\n");
	# Return output.
	out;
}

# Gompertz survivorship function (with c>0, it's Gompertz-Makeham
# An unused s parameter is specified, since the calling script doesn't
# know which function it's calling and give all four parameters to both
# Gompertz and Logistic functions.
survfng<-function(x,a,b,c,s=0){exp(-a*(exp(b*x)-1)/b -c*x);}
# Logistic survivorship function (with c>0, it's Logistic-Makeham
survfnl<-function(x,a,b,c,s){exp(-c*x)*(1+s*a*(exp(b*x)-1)/b)^(-1/s);}

# This is a quick little script that converts the simulated populations to
# a matrix that can be exported as a table and directly used by WinModest
# Usage: foo<-raw2wm(bar$r); write.table(foo,file='foo.txt',sep='\t')
raw2wm<-function(t){
	if(!is.list(t)){t<-list(t);}
	out<-c();
	for(i in 1:length(t)){
		irle<-rle(sort(t[[i]]));
		iwm<-cbind(irle$values,irle$lengths,i,1);
		out<-rbind(out,iwm);
	}
	out<-rbind(out,c(-1,NA,NA,NA));
	colnames(out)<-NULL;
	out;
}
