fitmods <- function(d,cx=rep(1,length(d)),models=c('g','gm','l','lm'),ext=F){
  models <- unique(c(modelinfo(models,'deps',extended=ext),models));
  # below, ext=F is intentionally there-- this is for param fitting, not model comparison
  t <- modelinfo(models,'tblcomp',extended=F)[,-(3:5)]; t <- rbind(data.frame(simple='e',complex='e'),t);
  nd <- length(d); lambda <- 1/mean(d);
  t <- cbind(t,LL=NA,a=NA,a.se=NA,b=NA,b.se=NA,c=NA,c.se=NA,s=NA,s.se=NA); 
  fc <- function(x){
    pf <- parent.frame(3); mo <- as.character(t$complex[x[1]]);
    if(all(t$complex[x]=='e')) {pf$t$LL[x] <- nd*(log(lambda)-1); pf$t$a[x] <- lambda; pf$t$a.se[x] <- lambda; return();}
    if(length(x)>1){
      # get the geometric mean of the nested model params, and turn the unused params into a small but non-zero number
      p <- apply(subset(t,complex %in% as.character(t$simple[x]))[,c('a','b','c','s')],2,function(y) gmean(ifelse(is.na(y),1e-14,y)));
      mi <- mo;
    } else {
      p <- subset(t,complex==as.character(t$simple[x]))[,c('a','b','c','s')];
      mi <- as.character(t$simple[x]);
    }
    o <- opsurv(d,model=mo,par=modpars(unlist(p),modeli=mi,modelo=mo,trim=T),lb=c(1e-14,0,0,0),cx=cx,getses=T);
    pf$t$LL[x] <- o$maximum; p <- modpars(o$estimate,modeli=mo,trim=F)[1:4];
    #browser();
    pse <- modpars(o$se,modeli=mo,trim=F)[1:4];
    pf$t$a[x] <- p[1]; pf$t$b[x] <- p[2]; pf$t$c[x] <- p[3]; pf$t$s[x] <- p[4];
    pf$t$a.se[x] <- pse[1]; pf$t$b.se[x] <- pse[2]; pf$t$c.se[x] <- pse[3]; pf$t$s.se[x] <- pse[4];
  }
  tapply(1:nrow(t),t$complex,fc); names(t)[2] <- 'model';
  return(unique(t[,-1]));
  # compare models using comptab-- if it's empty, defaults to left two columns of t; otherwise merge comptab with t
  # LR and pchi will be generated along with potentially AIC and BIC
  # if a choice table is specified, use it to decide which model is best
}
