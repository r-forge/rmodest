#include <R.h>
#include <Rmath.h>

void zpprob(double *x1, double *x2, double *n1, double *n2, int *ln, double *th, double *ans){
	int i;
	for(i=0;i< *ln; i++){
		ans[i] = choose(*n1,x1[i]) * choose(*n2,x2[i]) * pow(*th,(x1[i]+x2[i])) * pow((1- *th),(*n1+ *n2 - x1[i]- x2[i]));
	}
}

void runzpprobs(double *n1, double *n2, double *th, double *myans, double *sup, double *keep){
	int k; int l; int m; double mysup;// initialize x1 and x2 of length n1+1 and n2+2
	m=0; mysup = 0;
	for(l=0;l<= *n2;l++){
		for(k=0;k<= *n1;k++){
			// x1[m] = k; x2[m] = l;
			m++;
		}
	}
	//zpprob(*x1,*x2,*n1,*n2,m,*th,*myans);
	for(l=0;l<= m;l++){
		if(keep[l]){mysup += myans[l];}
	}
	*sup = mysup;
}
/*	# below for statement is candidate for porting to C
	# args: float *zptests, float *zp, float *thetas1, float *thetas2, float *step, 
	#	int *n1, int *n2, and the zpprob function that should also be ported
	for(i in 1:dim(output)[1]){
		sup<-c(); cat('.');
		keep<-abs(zptests)>abs(output$zp[i]);
		for(j in seq(thetas[1],thetas[2],step)){
			if(debug) cat("\nRunning zpprob on n1\'s and n2\'s.");
			zptemp<-outer(0:n1,0:n2,FUN="zpprob2",n1=n1,n2=n2,th=j);
			sup<-c(sup,sum(zptemp[keep]));
		}
		p<-c(p,max(sup));
	}
*/