`survgui`<-
function(dat=NULL){
	# make sure the environment can support this GUI
	require(tcltk);
	tclRequire('Tktable');
	# declare a bunch of tcl variables for later use
	a<-tktoplevel(); n<-1000;
	frame1<-tkframe(a,pady=5,padx=5);
	frame2<-tkframe(a);frame3<-tkframe(a);frame4<-tkframe(frame2);
	ar10<-tkfont.create(family='arial',size=10); ar8<-tkfont.create(family='arial',size=8);
	qincr<-c(.01,.05,.1,.2); tcl('set','qrse','boot'); penv<-environment(); cursor<-1;
	digits=4;
	oldwidth<-options('width'); options(width=250);
	loadFunc<-function(filename=NULL){ 
		tkconfigure(a,cursor='watch');
		if(is.null(filename)){
			filename<-tclvalue(tkgetOpenFile(initialdir=tclvalue('path'),
							 filetypes="{{R Data} {.rdata .Rdata}}"));
		}
		if(filename!=""){
			for(i in c('summ.out','lrnk.out','scor.out','qreg.out','qreg.sum','mort.out')){
				if(exists(i,envir=penv)){rm(list=i,envir=penv);}
			}
			ls(env=penv);
			for(i in c('plQregBut','plSurvBut','svSummBut','svScorBut','svQregBut',
				   'svMortBut','svplSurvBut','svplQregBut','svAllBut','prSummBut',
				   'prScorBut','prQregBut','prMortBut','prAllBut')){
				tkconfigure(get(i),state='disabled');
			}
			load(filename); objcount<-0;
			if(dev.cur()>1){dev.off(dev.cur());}
			if(exists('xtab')&exists('ytab')){
				array2tcl(xtab,vard1,1:dim(xtab)[1],0:1);
				array2tcl(ytab,vard2,1:dim(ytab)[1],0:1);
				tkconfigure(plSurvBut,state='active')
				objcount<-objcount+1;
			}
			if(exists('xynames')){
				tcl("set","name1",xynames[1]);tcl("set","name2",xynames[2]);
				objcount<-objcount+1;
			}
			if(exists('summ.out')){
				tkconfigure(svSummBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(svAllBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(prSummBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(prAllBut,state='active',foreground='black',
					    activeforeground='black');
				assign('summ.out',summ.out,env=penv)
				objcount<-objcount+1;
			}
			if(exists('lrnk.out')){assign('lrnk.out',lrnk.out,env=penv);objcount<-objcount+1;}
			if(exists('scor.out')){
				tkconfigure(svScorBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(svAllBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(prScorBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(prAllBut,state='active',foreground='black',
					    activeforeground='black');
				assign('scor.out',scor.out,env=penv)
				objcount<-objcount+1;
			}
			if(exists('qreg.out')){
				tkconfigure(svQregBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(svAllBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(prQregBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(prAllBut,state='active',foreground='black',
					    activeforeground='black');
				assign('qreg.out',qreg.out,env=penv)
				objcount<-objcount+1;
			}
			if(exists('qreg.sum')){
				tkconfigure(plQregBut,state='active',foreground='black',
					    activeforeground='black');
				assign('qreg.sum',qreg.sum,env=penv)
				objcount<-objcount+1;
			}
			if(exists('mort.out')){
				tkconfigure(svMortBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(svAllBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(prMortBut,state='active',foreground='black',
					    activeforeground='black');
				tkconfigure(prAllBut,state='active',foreground='black',
					    activeforeground='black');
				assign('mort.out',mort.out,env=penv)
				objcount<-objcount+1;
			}
		}
		tkfocus(a);
		tkconfigure(a,cursor='');
	}

	quitFunc<-function(){tkdestroy(a);}

	doGeneralFunc<-function(){
		# the [,1:2] subscript might be extraneous; remove it when nothing has broken in a while
		assign('xtab',na.exclude(tcl2nArray(vard1,1:n,0:1))[,1:2],env=penv);
		assign('ytab',na.exclude(tcl2nArray(vard2,1:n,0:1))[,1:2],env=penv);
		assign('x',xtab[xtab[,2]==1,1],env=penv);
		assign('y',ytab[ytab[,2]==1,1],env=penv);
		assign('xynames',c(tclvalue('name1'),tclvalue('name2')),env=penv);
		assign('tmain',paste(xynames,collapse=' vs. '),env=penv);
		assign('xy',c(x,y),env=penv);
		assign('lx',length(x),env=penv);assign('ly',length(y),env=penv);
		assign('group',factor(rep(xynames,c(lx,ly)),levels=xynames),env=penv);
	}

	doSummFunc<-function(){
		tkconfigure(a,cursor='watch');
		doGeneralFunc(); 
		sdx<-sd(x); sdy<-sd(y); mx<-mean(x); my<-mean(y); 
		qx<-round(quantile(xtab[,1],c(.5,.9)),0); qy<-round(quantile(ytab[,1],c(.5,.9)),0);
		qcix<-qci(xtab[,1],c(.5,.9)); qciy<-qci(ytab[,1],c(.5,.9));
		xyt<-t.test(log(x[x>0]),log(y[y>0]));
		xytab<-rbind(cbind(xtab,group=1),cbind(ytab,group=2));
		xylr<-surv2.logrank(Surv(xytab[,1],event=xytab[,2]),xytab[,3]);
		summ.out<-data.frame(NA,lx,NA,NA,ly,NA,NA);
		colnames(summ.out)<-c(paste(xynames[1],'lci'),xynames[1],paste(xynames[1],'uci'),
				      paste(xynames[2],'lci'),xynames[2],paste(xynames[2],'uci'),'p');
		rownames(summ.out)<-'n';
		summ.out<-rbind(summ.out,mean=c(sdx,mx,sdx,sdy,my,sdy,xyt$p.value));
		summ.out<-rbind(summ.out,median=c(qcix[1,1],qx[1],qcix[1,2],qciy[1,1],qy[1],qciy[1,2],
				NA)); #if qreg.out exists, get p from there
		summ.out<-rbind(summ.out,`90th percentile`=c(qcix[1,1],qx[1],qcix[1,2],
				qciy[1,1],qy[1],qciy[1,2],
				NA)); #if qreg.out exists, get p from there
		if(exists('scor.out')){
			summ.out['median','p']<-scor.out['0.5','p'];
			summ.out['90th percentile','p']<-scor.out['0.9','p'];
		}
		assign('summ.out',summ.out,env=penv); assign('lrnk.out',xylr,env=penv);
		prSummFunc();
		tkconfigure(svSummBut,state='normal',foreground='blue',activeforeground='blue');
		tkconfigure(svAllBut,state='normal',foreground='blue',activeforeground='blue');
		tkconfigure(prSummBut,state='normal',foreground='black',activeforeground='black');
		tkconfigure(prAllBut,state='normal',foreground='black',activeforeground='black');
		tkconfigure(doSummBut,foreground='black',activeforeground='black');
		tkconfigure(a,cursor='');
	}
	doScorFunc<-function(){
		tkconfigure(a,cursor='watch');
		doGeneralFunc();
		assign('scor.out',ezz(x,y,xynames,quant=seq(0,.95,.05)),env=penv);
		prScorFunc();
		tkconfigure(svScorBut,state='normal',foreground='blue',activeforeground='blue');
		tkconfigure(svAllBut,state='normal',foreground='blue',activeforeground='blue');
		tkconfigure(prScorBut,state='normal',foreground='black',activeforeground='black');
		tkconfigure(prAllBut,state='normal',foreground='black',activeforeground='black');
		tkconfigure(doSummBut,text='Update Report',foreground='blue',activeforeground='blue');
		tkconfigure(a,cursor='');
	}
	doQregFunc<-function(){
		tkconfigure(a,cursor='watch');
		doGeneralFunc(); qreg.sum<-0;class(qreg.sum)<-"try-error"; i<-1;
		while(class(qreg.sum)=="try-error"){
			if(i>length(qincr)){
				out<-"Unable to perform quantile regression.";
				next;
			}
			# qincr is a series of bandwidths, declared at the beginning
			# here, it is used to construct a sequence of tau values for rq
			# that don't crash it (and keep trying until it exhausts the
			# reasonable bandwidths)
			iqreg<-rq(log(xy)~group,tau=seq(qincr[i],1-qincr[i],qincr[i]));i<-i+1;
			qreg.sum<-try(summary(iqreg,se=tclvalue('qrse'))); cat('.');
		}
		if(class(qreg.sum)!="character"){
			qreg.out<-c(); qreg.sig=F;
			for(i in qreg.sum){
				qreg.out<-rbind(qreg.out,c(quantile=i$tau,i$coefficients[2,]));
			}
			tkconfigure(plQregBut,state='normal');
			tkconfigure(svQregBut,state='normal',foreground='blue',activeforeground='blue');
			tkconfigure(svAllBut,state='normal',foreground='blue',activeforeground='blue');
			tkconfigure(prQregBut,state='normal',foreground='black',activeforeground='black');
			tkconfigure(prAllBut,state='normal',foreground='black',activeforeground='black');
			assign('qreg.sum',qreg.sum,env=penv);
		}else{qreg.out<-qreg.sum;}
		assign('qreg.out',qreg.out,env=penv);
		tkconfigure(doSummBut,text='Update Report',foreground='blue',activeforeground='blue');
		prQregFunc();
		tkconfigure(a,cursor='');
	}
	doMortFunc<-function(){
		tkconfigure(a,cursor='watch');
		doGeneralFunc();
		assign('mort.out',findpars(xtab[,1],ytab[,1],cx=xtab[,2],cy=ytab[,2],summary=F),
		       env=penv);
		prMortFunc();
		tkconfigure(svMortBut,state='normal',foreground='blue',activeforeground='blue');
		tkconfigure(svAllBut,state='normal',foreground='blue',activeforeground='blue');
		tkconfigure(prMortBut,state='normal',foreground='black',activeforeground='black');
		tkconfigure(prAllBut,state='normal',foreground='black',activeforeground='black');
		tkconfigure(doSummBut,text='Update Report',foreground='blue',activeforeground='blue');
		tkconfigure(a,cursor='');
	}
	doAllFunc<-function(){
		doSummFunc(); doScorFunc(); doQregFunc(); doMortFunc();
	}
	plQregFunc<-function(){
		tkconfigure(a,cursor='watch');
		doGeneralFunc();
		plot(qreg.sum,parm=2,ols=F,ylim=c(-1,1),main=paste(tmain,"Quantile Regression"));
		qreg.sig<-na.exclude(qreg.out[qreg.out[,5]<.05,]);
		points(qreg.sig[,1],qreg.sig[,2],col='red',lwd=3);
		tkconfigure(svplSurvBut,state='disabled');
		tkconfigure(svplQregBut,state='normal',foreground='blue',activeforeground='blue');
		tkconfigure(a,cursor='');
	}
	plSurvFunc<-function(){
		tkconfigure(a,cursor='watch');
		doGeneralFunc();
		plotsurv(xtab[,1],ytab[,1],xtab[,2],ytab[,2],legend=c(xynames));
		tkconfigure(svplSurvBut,state='normal',foreground='blue',activeforeground='blue');
		tkconfigure(svplQregBut,state='disabled');
		tkconfigure(a,cursor='');
	}

	doSummBut<-tkbutton(frame2,text="Run",command=doSummFunc,font=ar10);
	doScorBut<-tkbutton(frame2,text="Run",command=doScorFunc,font=ar10);
	doQregBut<-tkbutton(frame2,text="Run",command=doQregFunc,font=ar10);
	doMortBut<-tkbutton(frame2,text="Run",command=doMortFunc,font=ar10);
	doAllBut<-tkbutton(frame2,text="Run All",command=doAllFunc,font=ar10);
	plQregBut<-tkbutton(frame2,text="Draw",
			    command=plQregFunc,font=ar10,state='disabled');
	plSurvBut<-tkbutton(frame2,text="Draw",command=plSurvFunc,font=ar10);

	delFunc<-function(){
		last<-as.numeric(tclvalue(tcl(textout,'index','end')));
		tcl(textout,'delete','1.0',paste(last,'end',sep='.'));
		assign('cursor',0,envir=penv);
	}

	prGeneralFunc<-function(d,title='',na='NA',dg=4){
		txt<-paste(capture.output(print(d,digits=dg)),collapse='\n');
		if(na!='NA'){txt<-gsub('\\bNA\\b',na,txt);}
		assign('cursor',as.numeric(tclvalue(tcl(textout,'index','end')))-2,envir=penv);
		tkinsert(textout,'end',paste(title,'\n\n'));
		tkinsert(textout,'end',txt);
		tkinsert(textout,'end','\n\n\n');
		tcl(textout,'yview',cursor);
	}

	prSummFunc<-function(){
		prGeneralFunc(summ.out,'Summary',na='  ');
		tempcursor<-cursor;
		tkinsert(textout,'end',paste('Log-rank','\n\n'));
		tkinsert(textout,'end',paste('Test statistic:',lrnk.out$stat,',','p =',lrnk.out$pval,
					     'based on',lrnk.out$nsim,'permutations.'));
		tkinsert(textout,'end','\n\n\n');
		if(exists('scor.out')){
			prGeneralFunc(scor.out[scor.out$p<0.05,],'Score Test: Significant Quantiles');
		}
		if(exists('qreg.out')){
			prGeneralFunc(qreg.out[!is.na(qreg.out[,5])&qreg.out[,5]<0.05,],
				      'Quantile Regression: Significant Quantiles');
		}
		if(exists('mort.out')){prGeneralFunc(mort.out,'Mortality Models');}
		tcl(textout,'yview',tempcursor);
	}
	prScorFunc<-function(){prGeneralFunc(scor.out,'Score Test');}
	prQregFunc<-function(){prGeneralFunc(qreg.out,'Quantile Regression');}
	prMortFunc<-function(){prGeneralFunc(mort.out,'Mortality Models');}
	prAllFunc<-function(){
 		tempcursor<-as.numeric(tclvalue(tcl(textout,'index','end')))-2;
		if(exists('summ.out')){prSummFunc();}
		if(exists('scor.out')){prScorFunc();}
		if(exists('qreg.out')){prQregFunc();}
		if(exists('mort.out')){prMortFunc();}
		tcl(textout,'yview',tempcursor);
	}

	prSummBut<-tkbutton(frame2,text="Show Output",command=prSummFunc,font=ar10,
			    state='disabled');
	prScorBut<-tkbutton(frame2,text="Show Output",command=prScorFunc,font=ar10,
			    state='disabled');
	prQregBut<-tkbutton(frame2,text="Show Output",command=prQregFunc,font=ar10,
			    state='disabled');
	prMortBut<-tkbutton(frame2,text="Show Output",command=prMortFunc,font=ar10,
			    state='disabled');
	prAllBut<-tkbutton(frame2,text="Show All Output",command=prAllFunc,font=ar10,
			    state='disabled');

	svGeneralFunc<-function(d,type,ext,filename=NULL,append=F){
		if(is.null(filename)){
			filetypes<-paste('{{',type,'} {',ext,'}}',sep='')
			filename<-tkgetSaveFile(initialdir=tclvalue('path'),filetypes=filetypes);
			filename<-tclvalue(filename);
		}
		if(filename!=""){
			#prGeneralFunc<-function(d,title='',na='NA',dg=4,to=NULL){
			switch(ext,
				'.txt'={write.table(d,file=filename,sep='\t',row.names=T,col.names=T,
					append=append);},
				'.rdata'={save(list=d,file=filename);},
				'.Rdata'={save(list=d,file=filename);},
				'.eps'={dev.copy2eps(file=filename);},
				'.dat'={write(d,file=filename,append=append)}
			)
			tkfocus(a);return(filename);
		} else {tkfocus(a);return(0)};
	}

	svSummFunc<-function(){
		if((filename<-svGeneralFunc(NULL,'Tab Delimited Test','.txt'))>0){
			svGeneralFunc('Summary\n','','.dat',filename=filename,append=T)
			svGeneralFunc(summ.out,'','.txt',filename=filename,append=T);
			svGeneralFunc(paste('\n\nLog-rank\n\nTest statistic:',lrnk.out$stat,',',
					    'p =',lrnk.out$pval,'based on',lrnk.out$nsim,
					    'permutations.'),'','.dat',filename=filename,append=T);
			if(exists('scor.out')){
				svGeneralFunc('\n\nScore Test: Significant Quantiles\n',
					      '','.dat',filename=filename,append=T);
				svGeneralFunc(scor.out[scor.out$p<0.05,],'','.txt',filename=filename,
					      append=T);
			};
			if(exists('qreg.out')){
				svGeneralFunc('\n\nQuantile Regression: Significant Quantiles\n',
					      '','.dat',filename=filename,append=T);
				svGeneralFunc(qreg.out[!is.na(qreg.out[,5])&qreg.out[,5]<0.05,],'','.txt',
					      filename=filename,append=T);
			}
			if(exists('mort.out')){
				svGeneralFunc('\n\nMortality Models\n',
					      '','.dat',filename=filename,append=T);
				svGeneralFunc(mort.out,'','.txt',filename=filename,append=T);
			};
			tkconfigure(svSummBut,foreground='black',activeforeground='black');
		}
	}

	svScorFunc<-function(){
		if(svGeneralFunc(scor.out,'Tab Delimited Text','.txt')>0){
			tkconfigure(svScorBut,foreground='black',activeforeground='black');
		}
	}

	svQregFunc<-function(){
		if(svGeneralFunc(qreg.out,'Tab Delimited Text','.txt')>0){
				tkconfigure(svQregBut,foreground='black',activeforeground='black');
		}
	}

	svMortFunc<-function(){
		if(svGeneralFunc(mort.out,'Tab Delimited Text','.txt')>0){
			tkconfigure(svMortBut,foreground='black',activeforeground='black');
		}
	}

	svAllFunc<-function(){
		doGeneralFunc();
		tosave<-c('xtab','ytab','xynames');
		if(exists('summ.out')){tosave<-c(tosave,'summ.out');}
		if(exists('lrnk.out')){tosave<-c(tosave,'lrnk.out');}
		if(exists('scor.out')){tosave<-c(tosave,'scor.out');}
		if(exists('qreg.out')){tosave<-c(tosave,'qreg.out');}
		if(exists('qreg.sum')){tosave<-c(tosave,'qreg.sum');}
		if(exists('mort.out')){tosave<-c(tosave,'mort.out');}
		if(svGeneralFunc(tosave,'R Data','.rdata')>0){
			tkconfigure(svAllBut,foreground='black',activeforeground='black');
		}
	}

	svplQregFunc<-function(){
		if(svGeneralFunc(NULL,'Encapsulated Postscript','.eps')>0){
			tkconfigure(svplQregBut,foreground='black',activeforeground='black');
		}
	}

	svplSurvFunc<-function(){
		if(svGeneralFunc(NULL,'Encapsulated Postscript','.eps')>0){
			tkconfigure(svplSurvBut,foreground='black',activeforeground='black');
		}
	}

	svSummBut<-tkbutton(frame2,text="Save Output (tab delimited spreadsheet)",command=svSummFunc,font=ar10,state='disabled');
	svScorBut<-tkbutton(frame2,text="Save Output (tab delimited spreadsheet)",command=svScorFunc,font=ar10,state='disabled');
	svQregBut<-tkbutton(frame2,text="Save Output (tab delimited spreadsheet)",command=svQregFunc,font=ar10,state='disabled');
	svMortBut<-tkbutton(frame2,text="Save Output (tab delimited spreadsheet)",command=svMortFunc,font=ar10,state='disabled');
	svAllBut<-tkbutton(frame2,text="Save All Output (.Rdata file)",command=svAllFunc,font=ar10,state='disabled');
	svplQregBut<-tkbutton(frame2,text="Save Image (.eps file)",command=svplQregFunc,font=ar10,state='disabled');
	svplSurvBut<-tkbutton(frame2,text="Save Image (.eps file)",command=svplSurvFunc,font=ar10,state='disabled');

	tkgrid(tklabel(frame2,text='Summary Report',font=ar10,anchor='e'),
	       doSummBut,prSummBut,svSummBut,sticky='we',padx=1,pady=3);
	tkgrid(tklabel(frame2,text='Survival Curves',font=ar10,anchor='e'),
	       tklabel(frame2,text=''),plSurvBut,svplSurvBut,sticky='we',padx=1,pady=3);
	tkgrid(tklabel(frame2,text='Score Test',font=ar10,anchor='e'),
	       doScorBut,prScorBut,svScorBut,sticky='we',padx=1,pady=3);
	tkgrid(tklabel(frame2,text='Quantile Regression',font=ar10,anchor='e'),
               doQregBut,prQregBut,svQregBut,sticky='we',padx=1,pady=3);
	tkgrid(tklabel(frame2,text='Quantile Regression Plot',font=ar10,anchor='e'),
	       tklabel(frame2,text=''),plQregBut,svplQregBut,sticky='we',padx=1,pady=3);
	tkgrid(tklabel(frame2,text='Mortality Models',font=ar10,anchor='e'),
	       doMortBut,prMortBut,svMortBut,sticky='we',padx=1,pady=3);
	tkgrid(tklabel(frame2,text='All Tests',font=ar10,anchor='e'),
	       doAllBut,prAllBut,svAllBut,sticky='we',padx=1,pady=3);

	# frame1 variables
	blank1<-tklabel(frame1,text='');blank2<-tklabel(frame1,text='');
#	blank1left<-tklabel(frame1left,text=' ');
# 	datatxt1<-tklabel(frame1left,text='Age at Death',font='Arial 8',padx=5,anchor='center');
# 	censtxt1<-tklabel(frame1left,text='Censor',font='Arial 8',padx=18,anchor='center');
# 	datatxt2<-tklabel(frame1right,text='Age at Death',font='Arial 8',padx=5,anchor='center');
# 	censtxt2<-tklabel(frame1right,text='Censor',font='Arial 8',padx=18,anchor='center');
	tcl("set","name1","Control"); tcl("set","name2","Experimental");
	name1txt<-tkentry(frame1,textvariable="name1",font=ar8,bg='white'); 
	name2txt<-tkentry(frame1,textvariable="name2",font=ar8,bg='white');
	vard1<-tclArray();vard2<-tclArray();
	vard1[[0,0]]<-'Age at Death'; vard1[[0,1]]<-'Censor';
	vard2[[0,0]]<-'Age at Death'; vard2[[0,1]]<-'Censor';
	for(i in 1:n){vard1[[i,0]]<-'-';vard1[[i,1]]<-1;vard2[[i,0]]<-'-';vard2[[i,1]]<-1;}
	# frame2 variables
	tcl("set","save",1);tcl("set","path",getwd());tcl("set","tscor",1); 
	tcl("set","tqreg",1);tcl("set","tmort",1);tcl("set","tsumm",1); 
	pathl<-tklabel(frame4,text='File Path:',anchor='e',font=ar10);
	# set up the input tables for the survival times in frame1
	d1<-tkwidget(frame1,'table',height=13,rows=n,cols=2,titlerows=0,titlecols=0,resizeborders='none',
		     colwidth=10,titlerows=1,yscrollcommand=function(...) tkset(yscr1,...));
	yscr1<-tkscrollbar(frame1,command=function(...)tkyview(d1,...));
	d2<-tkwidget(frame1,'table',height=13,rows=n,cols=2,titlerows=0,titlecols=0,resizeborders='none',
		     colwidth=10,titlerows=1,yscrollcommand=function(...) tkset(yscr2,...));
	yscr2<-tkscrollbar(frame1,command=function(...)tkyview(d2,...));
	# set up widgets in frame2
	pathtxt<-tkentry(frame4,textvariable="path",font=ar8,bg='white');
# 	savechk<-tkcheckbutton(frame2,variable="save", text='Save',font=ar10,anchor='w');
# 	scorchk<-tkcheckbutton(frame2,variable="tscor",text='Score Test',font=ar10,anchor='w');
# 	qregchk<-tkcheckbutton(frame2,variable="tqreg",text='Quant Reg',font=ar10,anchor='w');
# 	mortchk<-tkcheckbutton(frame2,variable="tmort",text='Mortality',font=ar10,anchor='w');
# 	summchk<-tkcheckbutton(frame2,variable="tsumm",text='Summary',font=ar10,anchor='w');
# 	runBut<-tkbutton(frame2,text="Run!",command=runFunc);
	dbBut<-tkbutton(frame4,text="Load from Database",command=function(){},padx=3,state='disabled');
	loadBut<-tkbutton(frame4,text="Load R Data",command=function() loadFunc(),padx=3);
	delBut<-tkbutton(frame4,text='Clear Output Window',command=delFunc,padx=3);
	quitBut<-tkbutton(frame4,text="Quit",command=quitFunc,padx=3);
	# set up frame 3
	xscr<-tkscrollbar(frame3,command=function(...)tkxview(textout,...),orient='horizontal');
	yscr3<-tkscrollbar(frame3,command=function(...)tkyview(textout,...));
	textout<-tktext(frame3,bg='white',font='courier 10',wrap='none',width=116,
			xscrollcommand=function(...) tkset(xscr,...),
			yscrollcommand=function(...) tkset(yscr3,...));
	# arrange the input tables in frame1
	tkgrid(name1txt,blank1,name2txt,blank1);
#	tkgrid(datatxt1,censtxt1,columnspan=2);
#	tkgrid(datatxt2,censtxt2,columnspan=1);tkgrid(frame1left,blank2,frame1right);
	tkgrid(d1,yscr1,d2,yscr2,columnspan=1);
	tkgrid.configure(yscr1,sticky='nsw');
	tkgrid.configure(yscr2,sticky='nsw');
	tkconfigure(d1,variable=vard1,background='white',selectmode='extended',
		    rowseparator="\"\n\"",multiline=0);
	tkconfigure(d2,variable=vard2,background='white',selectmode='extended',
		    rowseparator="\"\n\"",multiline=0);
	# arrange widgets in frame2
 	tkgrid(pathl,pathtxt,sticky='we',padx=2,pady=2,columnspan=1);
# 	tkgrid(savechk,sticky='we',columnspan=2);
# 	tkgrid(scorchk,sticky='we',columnspan=2);
# 	tkgrid(qregchk,sticky='we',columnspan=2);
# 	tkgrid(mortchk,sticky='we',columnspan=2);
# 	tkgrid(summchk,sticky='we',columnspan=2);
# 	tkgrid(runBut,sticky='we',padx=50,pady=2,columnspan=2);
	tkgrid(dbBut,loadBut,delBut,quitBut,pady=2,padx=5,columnspan=1);
	tkwm.title(a,'Survomatic v1.2.3');
	if(!is.null(dat)){loadFunc(dat);}
	# arrange widgets in frame3
	tkgrid(textout,yscr3); tkgrid(xscr); tkgrid.configure(textout,sticky='we');
	tkgrid.configure(xscr,sticky='we'); tkgrid.configure(yscr3,sticky='ns');
	#tkinsert(textout,'end',rep(' ',options('width')[[1]]));
	# arrange frames
	tkgrid(frame4,columnspan=4);
	tkgrid(frame1,frame2,padx=2,pady=2);
	tkgrid(frame3,columnspan=2,sticky='nwe');
	browser();
	options(width=oldwidth[[1]]);
}