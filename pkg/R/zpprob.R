`zpprob` <-
function(x1,x2,n1,n2,th){
	choose(n1,x1)*choose(n2,x2)*th^(x1+x2)*(1-th)^(n1+n2-x1-x2);
}

