`choosemodel`<-function(x){
    # This function takes a vector of concatenated model abbreviations where the first abbreviation 
    # of each pair is the simple model and the second is the complex one
    w <- 'we' %in% x; g <- 'ge' %in% x; gm <- 'gmg' %in% x; l <- 'lg' %in% x;
    lmg <- 'lmg' %in% x; lmgm <- 'lmgm' %in% x; lml <- 'lml' %in% x;
    out<-c();
    if(w) {out<-c(out,'w');}
    if(!g) {warning('Either there is an error in the data or this data does not follow the Gompertz distribution. It may instead follow an exponential distribution. Do not rely on Gompertz-based models for this data.'); return(out);}
    if(!(gm | l)){if(lmg){return(c(out,'lm'));} else {return(c(out,'g'));}} else {
	if(gm & l){ if(lmgm & lml) {return(c(out,'lm'));} else {return(c(out,'gm','l','lm'));} } else {
	    if(gm & !lmgm){return(c(out,'gm'));} else {if(l & !lml){return(c(out,'l'));}}
	    if((gm & lmgm)|(l & lml)){return(c(out,'gm','l','lm'));}
	}
    }
}