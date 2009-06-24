`tcl2nArray`<-
function(d,rows,cols){
	lr<-length(rows);lc<-length(cols);
	out<-matrix(nrow=lr,ncol=lc);
	for(i in 1:lr){
		for(j in 1:lc){
			out[i,j]<-as.numeric(d[[rows[i],cols[j]]]);
		}
	}
	return(out);
}

`array2tcl`<-
function(d,out,rows,cols){
	for(i in 1:length(rows)){
		for(j in 1:length(cols)){
			out[[rows[i],cols[j]]]<-d[i,j];
		}
	}
}