plotkm<-function(data,surv,group=1,pos="bottomleft",units="months")
{
  if(class(group)=="numeric"){  
    kfit<-survfit(as.formula(paste("Surv(",surv[1],",",surv[2],")~1",sep="")),data=data)
    sk<-summary(kfit)$table
    levelnames<-paste("N=",sk[1], ", Events=",sk[4]," (",round(sk[4]/sk[1],2)*100,"%)",sep="")
    main<-paste("KM-Curve for ",nicename(surv[2]),sep="")
    
  }else if(length(group)>1){
    return("Currently you can only stratify by 1 variable")
  }else{      
    lr<-survdiff(as.formula(paste("Surv(",surv[1],",",surv[2],")~", paste(group,collapse="+"),sep="")),data=data)
    lrpv<-1-pchisq(lr$chisq, length(lr$n)- 1)
    levelnames<-levels(data[,group])
    kfit<-survfit(as.formula(paste("Surv(",surv[1],",",surv[2],")~", paste(group,collapse="+"),sep="")),data=data)
    main<-paste("KM-Curve for ",nicename(surv[2])," by ", nicename(group),sep="")
    levelnames<-sapply(1:length(levelnames), function(x){paste(levelnames[x]," n=",lr$n[x],sep="")})
    
  }  
  
  plot(kfit,mark.time=T, col=1:length(levelnames),xlab=paste("Time (",cap(units),")",sep=""),
       ylab="Suvival Probability ",conf.int=F,cex=1.1,
       main=main)
  
  if(class(group)=="numeric"){legend(pos,levelnames,col=1:length(levelnames),bty="n")
  }else{ legend(pos,c(levelnames,paste("p-value=",pvalue(lrpv)," (Log Rank)",sep="")),
                col=c(1:length(levelnames),"white"),lty=1,bty="n")}
}


etsum<- function(data,surv,group=1,times=c(12,24)){
  
  kfit<-summary(survfit(as.formula(paste("Surv(",surv[1],",",surv[2],")~",group,sep=""))  ,data=data))
  kfit2<-summary(survfit(as.formula(paste("Surv(",surv[1],",",surv[2],")~",group,sep="")) ,data=data),times=times)
  tab<-as.data.frame(cbind(strata=as.character(kfit2$strata),times=kfit2$time,SR=paste(round(kfit2$surv*100,0)," (",round(kfit2$lower*100,0),"-",round(kfit2$upper*100,0),")",sep="")))
  tbl<-kfit2$table    
  
  if(class(group)!="numeric"){
    med=by(data,data[,group],function(x) median(x[,surv[1]]))
    min=by(data,data[,group],function(x) min(x[,surv[1]]))
    max=by(data,data[,group],function(x) max(x[,surv[1]]))
    survtimes<-data.frame(strata=as.character(kfit$strata),kfit$time)   
    minst<-round(as.numeric(by(survtimes,survtimes$strata,function(x) min (x[,2]))),1)
    maxst<-round(as.numeric(by(survtimes,survtimes$strata,function(x) max (x[,2]))),1)    
    tab<-cast(tab, strata ~ times)
    names<-names(tab)
    tab<-data.frame(tab)
    names(tab)<-names
    tab[,1]<-levels(data[,group])
    if(length(times)>1){
      indx<-c(0,sapply(sort(as.numeric(names(tab)[-1])),function(x){which(as.numeric(names(tab)[-1])==x)}))+1
      tab<-tab[,indx]
      tab<-tab[c(2:length(tab),1)]
    }else{tab<-tab[c(2,1)]}
    noeventsindx<-ifelse(length(which(tbl[,4]==0))!=0,
                         which(tbl[,4]==0),NA)
    if(!is.na(noeventsindx)){
      for(i in noeventsindx){
        if(i==1){
          minst<-c(0,minst)
          maxst<-c(0,maxst)
        }else if(i>length(minst)){
          minst<-c(minst,0)
          maxst<-c(maxst,0)
        }else{
          minst<-c(minst[1:i-1],0,minst[i:length(minst)])
          maxst<-c(maxst[1:i-1],0,maxst[i:length(maxst)])
        }}}     
    
    
    tab<-cbind("n"=tbl[,1],"Events"=tbl[,4], "MedKM"=round(tbl[,5],1),
               "LCI"=round(tbl[,6],1), "UCI"=round(tbl[,7],1),
               "MedFU"=round(as.numeric(med),1),
               "MinFU"=round(as.numeric(min),1),"MaxFU"=round(as.numeric(max),1),
               "MinET"=minst,"MaxET"=maxst,tab)
    rownames(tab)<-NULL
  }else{
    med=median(data[,surv[1]])
    min=min(data[,surv[1]])
    max=max(data[,surv[1]])
    if(length(times)>1){
      tab<-data.frame(t(tab))
      rownames(tab)<-NULL      
      names(tab)<-as.numeric(as.matrix(tab[1,]))
      tab<-tab[-1,]      
    }else{
      rownames(tab)<-NULL
      names(tab)[2]<-times
      tab<-tab[-1]
    } 
    tab<-cbind("n"=tbl[1],"Events"=tbl[4],"MedKM"=round(tbl[5],1),"LCI"=round(tbl[6],1), "UCI"=round(tbl[7],1),
               "MedFU"=round(as.numeric(med),1),"MinFU"=round(as.numeric(min),1),"MaxFU"=round(as.numeric(max),1),
               "MinET"=round(min(kfit$time),1),"MaxET"=round(max(kfit$time),1),tab)
    rownames(tab)<-NULL
  }
  return(tab)
}

petsum<-function(data,surv,group="1",times=c(12,14),units="months"){
  t<-etsum(data,surv,group,times)
  
  #plotkm(nona,surv,group)
  
  names<-names(t)
  if("strata"%in% names){
    strta<-sapply(t[,"strata"], function(x) paste(x,": ",sep=""))
    offset<-2
    ofst<-1
  }else{
    strta=matrix(c("",""))
    offset<-1
    ofst<-0
  }
  
  
  out<-sapply(seq_len(nrow(t)),function(i){
    
    if(is.na(t[i,3])) {km<-paste("The KM median event time has not been achieved due to lack of events.",sep="")
    }else if (!is.na(t[i,5])){km<-paste("The KM median event time is ",t[i,3]," with 95",sanitizestr("%")," confidence Interval (",t[i,4],",",t[i,5],").",sep="")
    }else{km<-paste("The KM median event time is ",t[i,3]," with 95",sanitizestr("%")," confidence Interval (",t[i,4],",",t[i,10],").",sep="")}
    
    # if at least one event
    if(t[i,2]!=0){
      flet<-paste(" The first and last event times occurred at ",t[i,9],
                  " and ",t[i,10]," ",units," respectively. ",sep="")
      if(ncol(t)>11+ofst){
        ps<-paste("The ",paste(names[11:(ncol(t)-offset)], collapse=", ")," and ",names[ncol(t)-ofst], " " , substring(units,1,nchar(units)-1),
                  " probabilities of 'survival' and their 95",sanitizestr("%")," confidence intervals are ",
                  paste(sapply(t[i,11:(ncol(t)-offset)],function(x) paste(x)),collapse=", ")," and ", t[i,ncol(t)-ofst], " percent.",sep="")
        
      }else{
        ps<-paste("The ",names[11]," ", substring(units,1,nchar(units)-1),
                  " probability of 'survival' and 95",sanitizestr("%")," confidence interval is ",
                  t[i,11]," percent.",sep="")  
      }
      #if no events
    }else{
      km=""
      ps=""
      flet=""
    }
    
    
    out<-paste(lbld(sanitizestr(nicename(strta[i])))," There are ",t[i,1]," patients. There were ",t[i,2],
               " (",round(100*t[i,2]/t[i,1],0),sanitizestr("%"),") events. The median and range of the follow-up times is ",
               t[i,6]," (",t[i,7],"-",t[i,8],") ",units,". ", km, flet,ps,sep="")
    cat("\n",out,"\n")
  })
}
covsum<-function(data,covs,sigcov,numobs=NULL,markup=T,sanitize=T,nicenames=T){
  
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity    
  }
  if(!sanitize) sanitizestr<-identity  
  if(!nicenames) nicename<-identity
  
  levels<-levels(data[,sigcov])
  levels<-c(list(levels),as.list(levels))
  N=nrow(data)
  nsigcov<-c(sum(table(data[,sigcov])),table(data[,sigcov]))
  out<-sapply(covs,function(cov){
    n<-sum(table(data[,cov]))
    
    #Set up the first coulmn
    factornames<-NULL
    if(is.null(numobs[[cov]]))  numobs[[cov]]<-nsigcov
    if(numobs[[cov]][1]-n>0) factornames<-c(factornames,"Missing")
    
    
    
    
    #if the covariate is a factor
    if(is.factor(data[,cov])){
      factornames<-c(levels(data[,cov]),factornames)    
      p<-lpvalue(fisher.test(data[,sigcov],data[,cov])$p.value)
      
      
      
      #set up the main columns
      onetbl<-mapply(function(sublevel,N){
        missing<-NULL
        subdata<-subset(data,subset=data[,sigcov]%in%sublevel)
        table<-table(subdata[,cov])
        tbl<-table(subdata[,cov])
        n<-sum(tbl)
        prop<-round(tbl/n,2)*100
        prop<-sapply(prop,function(x){if(!is.nan(x)){x} else{0}})
        tbl<-mapply(function(num,prop){paste(num," (",prop,")",sep="")},tbl,prop)
        if("Missing" %in% factornames) missing<-N-n 
        tbl<-c(tbl,missing)    
        return(tbl)
      },levels,numobs[[cov]])
      
      #if the covariate is not a factor
    }else{
      #setup the first column
      factornames<-c("Mean (sd)", "Median (Min,Max)",factornames)
      p<-lpvalue(anova(lm(data[,cov]~data[,sigcov]))[5][[1]][1])
      
      #set up the main columns
      onetbl<-mapply(function(sublevel,N){
        missing<-NULL
        subdata<-subset(data,subset=data[,sigcov]%in%sublevel)
        summary<-round(summary(subdata[,cov]),1)
        meansd<-paste(summary[4]," (", round(sd(subdata[,cov],na.rm=T),1),")",sep="")
        mmm<-paste(summary[3]," (",summary[1],",",summary[6],")",sep="")
        
        #if there is a missing in the whole data
        if("Missing" %in% factornames){          
          n<-sum(table(subdata[,cov]))
          missing<-N-n }       
        tbl<-c(meansd,mmm,missing)
        
        return(tbl)}   
                     ,levels,numobs[[cov]])}
    
    
    #Add the first column to the main columns and get the matrix ready for later
    factornames<-addspace(sanitizestr(nicename(factornames)))    
    onetbl<-cbind(factornames,onetbl)
    onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),rep("",length(levels[[1]])+1)),onetbl)
    onetbl<-cbind(onetbl,c(p,rep("",nrow(onetbl)-1)))
    rownames(onetbl)<-NULL
    colnames(onetbl)<-NULL
    return(onetbl)})
  table<-do.call(rbind.data.frame, out)
  rownames(table)<-NULL
  colnames(table)<-c("Covariate",paste("Full Sample (n=",N,")",sep=""),
                     mapply(function(x,y){paste(x," (n=",y,")",sep="")},
                            levels(data[,sigcov]),table(data[,sigcov])),"p-value (indep)")
  return(table)
}

pcovsum<-function(data,covs,sigcov,numobs,latex=F){
  if(!latex){
    print.xtable(xtable(covsum(data,covs,sigcov,numobs)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{  
    
    print.xtable(xtable(covsum(data,covs,sigcov,numobs)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }}
uvsum<-function(data,surv,covs,strat="1",markup=T,sanitize=T,nicenames=T){
  
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity    
  }
  if(!sanitize) sanitizestr<-identity  
  if(!nicenames) nicename<-identity
  
  
  names<-names(data)
  if(class(strat)!="numeric") {strat<-sapply(strat,function(stra){paste("strata(",stra,")",sep="")})
  }else strat<-1
  m1<-coxph(as.formula(paste(paste("Surv(",surv[1],",",surv[2],")",sep=""),"~",paste(strat,collapse="+"),sep="")),data=data)  
  out<-lapply(covs,function(cov){
    #cat(cov,"\n")
    m2<-update(m1,as.formula(paste(".~.+",cov,sep="")))
    
    if(is.factor(data[,cov])){
      levelnames<-sapply(sapply(sapply(levels(factor(data[,cov])),nicename),sanitizestr),addspace)
      cov<-lbld(sanitizestr(nicename(cov)))
      title<-NULL
      body<-NULL    
      hazardratio<-c("Reference",apply(matrix(summary(m2)$conf.int[,c(1,3,4)],ncol=3),1,psthr))    
      pvalue<-c("",sapply(summary(m2)$coef[,5],lpvalue))
      title<-c(cov,"","",lpvalue(summary(m2)$logtest[3]))
      
      if(length(levelnames)==2){
        body<-cbind(levelnames,hazardratio,c("",""),c("",""))    
      }else{
        body<-cbind(levelnames,hazardratio,pvalue,rep("",length(levelnames)))      
      }
      out<-rbind(title,body)
      rownames(out)<-NULL
      colnames(out)<-NULL
      return(list(out,nrow(out)))
    }else
    {
      cov<-lbld(sanitizestr(nicename(cov)))    
      out<-matrix(c(cov,psthr(summary(m2)$conf.int[,c(1,3,4)]),"",lpvalue(summary(m2)$logtest[3])),ncol=4)
      return(list(out,nrow(out)))}})
  table<-lapply(out,function(x){return(x[[1]])})
  table<-do.call(rbind.data.frame, table)  
  colnames(table)<-sapply(c("Covariate","HR(95\\%CI)","p-value","Global p-value"),lbld)
  return(table)
  
  index<-unlist(lapply(out,function(x){return(x[[2]])}))  
  lineindex<-rep(-1,length(index))
  for(i in 1:length(index)){lineindex[i]<-sum(index[c(1:i)])}
  lineindex<-c(-1,0,lineindex)
  if(output=="latex") return(list(table,lineindex))
  else return(table)
}

puvsum<-function(data,surv,covs,strata=1,latex=F){
  if(!latex){
    print.xtable(xtable(uvsum(data,surv,covs,strata)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{
    print.xtable(xtable(uvsum(data,surv,covs,strata)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }
  
}

mvsum<-function(data,model,markup=T,sanitize=T,nicenames=T){
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity    
  }
  if(!sanitize) sanitizestr<-identity  
  if(!nicenames) nicename<-identity  
  
  call<-paste(deparse(summary(model)$call),collapse="")
  call<-unlist(strsplit(call,"~",fixed=T))[2]
  call<-unlist(strsplit(call,",",fixed=T))[1]
  call<-unlist(strsplit(call,"+",fixed=T))
  call<-unlist(strsplit(call,"*",fixed=T))
  call<-unique(call)
  call<-call[which(is.na(sapply(call,function(cov){charmatch("strata(",cov)}))==T)]
  call<-gsub("\\s","", call)
  betanames<-attributes(summary(model)$coefficients)$dimnames[[1]]
  indx<-as.vector(sapply(betanames,function(string){
    
    indx<-which(sapply(call,function(cov){charmatch(cov,string)})==1)
    if(length(indx)==1) return(indx)
    #If one  facorname is a subset of another
    indx2<-which.max(sapply(call[indx],nchar))
    if(length(indx2)==1) return(indx[indx2])
    indx3<-which(sapply(call[indx2],function(c){substr(betaname,1,nchar(c))==c}))
    if(length(indx3)==1)  return(call[indx[indx2[indx3]]])  
    return(-1)
  }))
  
  if(min(indx)==-1) return("error") 
  
  
  
  
  
  y<-betaindx(indx)
  
  out<-lapply(y,function(covariateindex){
    
    #Get attribute names and split by ineractions
    betaname<-betanames[covariateindex]
    betaname<-strsplit(betaname,":",fixed=T)
    #get the covariate names
    oldcovname<-covnm(betaname[[1]],call)
    
    #get the levelnames
    levelnames<-unlist(lapply(betaname,function(level){
      paste(mapply(function(lvl,cn){result<-unlist(strsplit(lvl,cn,fixed=T))[2]
                                    out<-ifelse(is.na(result),cn,result)},level,oldcovname),collapse=":")}))
    levelnames<-addspace(sanitizestr(nicename(levelnames)))
    covariatename<-lbld(sanitizestr(nicename(paste(oldcovname,collapse=":"))))
    globalpvalue<-try(wald.test(b=coef(model)[covariateindex],Sigma=vcov(model)[covariateindex,covariateindex],Terms=seq_along(covariateindex))$result$chi2[3]);
    if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
    #globalpvalue=NA
    globalpvalue<-lpvalue(globalpvalue)
    reference=NULL
    title=NULL
    body=NULL
    
    #if not interaction
    hazardratio<-c(apply(matrix(summary(model)$conf.int[covariateindex,c(1,3,4)],ncol=3),1,psthr))
    pvalues<-c(sapply(summary(model)$coef[covariateindex,5],lpvalue))
    
    if(length(betaname[[1]])==1){
      #if cts
      if(!is.factor(data[,oldcovname])){
        
        title<-c(nicename(covariatename),hazardratio,"",globalpvalue)
      }else if(length(levelnames)==1){
        title<-c(covariatename,"","",globalpvalue)        
        if(!is.null(data)) reference<-c(addspace(sanitizestr(names(table(data[,which(names(data)==oldcovname)]))[1])),"reference","","")
        body<-c(levelnames,hazardratio,"","")       
        
        
      }else{
        if(!is.null(data)) reference<-c(addspace(sanitizestr(names(table(data[,which(names(data)==oldcovname)]))[1])),"reference","","")
        title<-c(covariatename,"","",globalpvalue)
        body<-cbind(levelnames,hazardratio,pvalues,rep("",length(levelnames)))     
        
        #if interaction
      }}else{
        if(length(levelnames)!=1){
          title<-c(covariatename,"","",globalpvalue)
          body<-cbind(levelnames,hazardratio,pvalues,rep("",length(levelnames)))
        }else{
          title<-c(covariatename,hazardratio,"",globalpvalue)        
        }
      }    
    
    out<-rbind(title,reference,body)
    rownames(out)<-NULL
    colnames(out)<-NULL
    return(list(out,nrow(out)))
  })
  table<-lapply(out,function(x){return(x[[1]])})
  index<-unlist(lapply(out,function(x){return(x[[2]])}))  
  table<-do.call(rbind.data.frame, table)
  colnames(table)<-sapply(c("Covariate",
                            paste(sanitizestr("HR(95%CI)"),sep=""),"p-value","Global p-value"),lbld)
  return(table)
  lineindex<-rep(-1,length(index))
  #for(i in 1:length(index)){lineindex[i]<-sum(index[c(1:i)])}
  #lineindex<-c(-1,0,lineindex)
  if(output=="latex") return(list(table,index))
  else return(table)  
}

pmvsum<-function(data,model){
  print.xtable(xtable(mvsum(data,model)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
}

makedocx<-function(dir,fname,imwd,pdwd){
  oldwd<-getwd()
  setwd(imwd)
  command<-paste("mogrify -path ", dir,"figure\\ ", "-format png ", dir, "figure\\*.pdf",sep="" )
  shell(command)
  setwd(pdwd)
  command<-paste("pandoc -o ",dir,fname,".docx ",dir,fname,".tex ",
                 "--default-image-extension=png",sep="")
  shell(command)
  setwd(oldwd)
}