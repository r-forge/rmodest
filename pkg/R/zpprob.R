`zpprob` <-
function(x1,x2,n1,n2,th){
	#choose(n1,x1)*choose(n2,x2)*th^(x1+x2)*(1-th)^(n1+n2-x1-x2);
	exp(lchoose(n1,x1)+lchoose(n2,x2)+(x1+x2)*log(th)+(n1+n2-x1-x2)*log(1-th));
}

