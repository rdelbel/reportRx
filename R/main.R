plotkm<-function(data,response,group=1,pos="bottomleft",units="months",CI=F,legend=T)
{
  if(class(group)=="numeric"){  
    kfit<-survfit(as.formula(paste("Surv(",response[1],",",response[2],")~1",sep="")),data=data)
    sk<-summary(kfit)$table
    levelnames<-paste("N=",sk[1], ", Events=",sk[4]," (",round(sk[4]/sk[1],2)*100,"%)",sep="")
    main<-paste("KM-Curve for ",nicename(response[2]),sep="")
    
  }else if(length(group)>1){
    return("Currently you can only stratify by 1 variable")
  }else{      
    lr<-survdiff(as.formula(paste("Surv(",response[1],",",response[2],")~", paste(group,collapse="+"),sep="")),data=data)
    lrpv<-1-pchisq(lr$chisq, length(lr$n)- 1)
    levelnames<-levels(data[,group])
    kfit<-survfit(as.formula(paste("Surv(",response[1],",",response[2],")~", paste(group,collapse="+"),sep="")),data=data)
    main<-paste("KM-Curve for ",nicename(response[2])," stratified by ", nicename(group),sep="")
    levelnames<-sapply(1:length(levelnames), function(x){paste(levelnames[x]," n=",lr$n[x],sep="")})
    
  }  
  
  
  plot(kfit,mark.time=T, lty=1:length(levelnames),xlab=paste("Time (",cap(units),")",sep=""),
       ylab="Suvival Probability ",cex=1.1, conf.int=CI,
       main=main)
  
  
  if(legend){
    if(class(group)=="numeric"){legend(pos,levelnames,lty=1:length(levelnames),bty="n")
    }else{ legend(pos,c(levelnames,paste("p-value=",pvalue(lrpv)," (Log Rank)",sep="")),
                  col=c(rep(1,length(levelnames)),"white"),lty=1:(length(levelnames)+1),bty="n")}
  }
}

etsum<- function(data,response,group=1,times=c(12,24)){
  kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep=""))  ,data=data))
  kfit2<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep="")) ,data=data),times=times)
  tab<-as.data.frame(cbind(strata=as.character(kfit2$strata),times=kfit2$time,SR=paste(round(kfit2$surv*100,0)," (",round(kfit2$lower*100,0),"-",round(kfit2$upper*100,0),")",sep="")))
  tbl<-kfit2$table    
  
  if(class(group)!="numeric"){
    med=by(data,data[,group],function(x) median(x[,response[1]],na.rm=T))
    min=by(data,data[,group],function(x) min(x[,response[1]],na.rm=T))
    max=by(data,data[,group],function(x) max(x[,response[1]],na.rm=T))
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
    med=median(data[,response[1]],na.rm=T)
    min=min(data[,response[1]],na.rm=T)
    max=max(data[,response[1]],na.rm=T)
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

petsum<-function(data,response,group=1,times=c(12,14),units="months"){
  t<-etsum(data,response,group,times)
  
  #plotkm(nona,response,group)
  
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

<<<<<<< HEAD
covsum<-function(data,covs,maincov=NULL,numobs=NULL,markup=T,sanitize=T,nicenames=T){
=======
covsum<-function(data,covs,maincov,numobs=NULL,markup=T,sanitize=T,nicenames=T){
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
  
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity    
  }
  if(!sanitize) sanitizestr<-identity  
  if(!nicenames) nicename<-identity
<<<<<<< HEAD
  if(!is.null(maincov)){
=======
  
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
  levels<-names(table(data[,maincov]))
  levels<-c(list(levels),as.list(levels))
  }else{
    levels<-"NOMAINCOVNULLNA"
  }
  N=nrow(data)
<<<<<<< HEAD
  if(!is.null(maincov)){
  nmaincov<-c(sum(table(data[,maincov])),table(data[,maincov]))
  }else{
    nmaincov<-N
    p<-NULL
  }
=======
  nmaincov<-c(sum(table(data[,maincov])),table(data[,maincov]))
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
  out<-lapply(covs,function(cov){
    ismiss=F
    n<-sum(table(data[,cov]))
    
    #Set up the first coulmn
    factornames<-NULL
    if(is.null(numobs[[cov]]))  numobs[[cov]]<-nmaincov
    if(numobs[[cov]][1]-n>0) {ismiss=T
      factornames<-c(factornames,"Missing")
    }
    #if the covariate is a factor
    if(is.factor(data[,cov])){
<<<<<<< HEAD
      factornames<-c(levels(data[,cov]),factornames)
      if(!is.null(maincov)){
      p<-try(lpvalue(fisher.test(data[,maincov],data[,cov])$p.value))
      if(class(p)=="try-error") p<-chisq.test(data[,maincov],data[,cov])$p.value
      p<-lpvalue(p)
      } 
=======
      factornames<-c(levels(data[,cov]),factornames)    
      p<-try(lpvalue(fisher.test(data[,maincov],data[,cov])$p.value))
      if(class(p)=="try-error") p<-chisq.test(data[,maincov],data[,cov])$p.value
      p<-lpvalue(p) 
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
        
      
      #set up the main columns
      onetbl<-mapply(function(sublevel,N){
        missing<-NULL
<<<<<<< HEAD
       if(sublevel[1]!="NOMAINCOVNULLNA"){
        subdata<-subset(data,subset=data[,maincov]%in%sublevel)
        }else{
          subdata<-data
        }
=======
        subdata<-subset(data,subset=data[,maincov]%in%sublevel)
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
        table<-table(subdata[,cov])
        tbl<-table(subdata[,cov])
        n<-sum(tbl)
        prop<-round(tbl/n,2)*100
        prop<-sapply(prop,function(x){if(!is.nan(x)){x} else{0}})
        tbl<-mapply(function(num,prop){paste(num," (",prop,")",sep="")},tbl,prop)       
        if(ismiss) missing<-N-n 
        tbl<-c(tbl,lbld(missing))       
        return(tbl)
      },levels,numobs[[cov]])
      
      #if the covariate is not a factor
    }else{
      #setup the first column
      factornames<-c("Mean (sd)", "Median (Min,Max)",factornames)
<<<<<<< HEAD
      if(!is.null(p)){
      p<-try(anova(lm(data[,cov]~data[,maincov]))[5][[1]][1])
      if(class(p)=="try-error") p<-NA
      p<-lpvalue(p)}
=======
      p<-try(anova(lm(data[,cov]~data[,maincov]))[5][[1]][1])
      if(class(p)=="try-error") p<-NA
      p<-lpvalue(p)
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
      
      
      #set up the main columns
      onetbl<-mapply(function(sublevel,N){
        missing<-NULL
<<<<<<< HEAD
       if(sublevel[1]!="NOMAINCOVNULLNA"){
        subdata<-subset(data,subset=data[,maincov]%in%sublevel)
        }else{subdata<-data}
=======
        subdata<-subset(data,subset=data[,maincov]%in%sublevel)
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
        summary<-round(summary(subdata[,cov]),1)
        meansd<-paste(summary[4]," (", round(sd(subdata[,cov],na.rm=T),1),")",sep="")
        mmm<-paste(summary[3]," (",summary[1],",",summary[6],")",sep="")
        
        
        #if there is a missing in the whole data
        if(ismiss){          
          n<-sum(table(subdata[,cov]))
          missing<-N-n 
        }
        tbl<-c(meansd,mmm,lbld(missing))
        
        return(tbl)}   
                     ,levels,numobs[[cov]])}
    
    #Add the first column to the main columns and get the matrix ready for later
    factornames<-addspace(sanitizestr(nicename(factornames)))    
    onetbl<-cbind(factornames,onetbl)
    
<<<<<<< HEAD
    if(!is.null(maincov)){
      onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),rep("",length(levels[[1]])+1)),onetbl)
      onetbl<-cbind(onetbl,c(p,rep("",nrow(onetbl)-1)))
    }else{
      onetbl<-rbind(lbld(sanitizestr(nicename(cov))),onetbl)
    }
=======
    onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),rep("",length(levels[[1]])+1)),onetbl)
    onetbl<-cbind(onetbl,c(p,rep("",nrow(onetbl)-1)))
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
    rownames(onetbl)<-NULL
    colnames(onetbl)<-NULL
    return(onetbl)})
  
  table<-do.call(rbind.data.frame, out)
 
  rownames(table)<-NULL
  if(!is.null(maincov)){
  colnames(table)<-c("Covariate",paste("Full Sample (n=",N,")",sep=""),
                     mapply(function(x,y){paste(x," (n=",y,")",sep="")},
                            names(table(data[,maincov])),table(data[,maincov])),"p-value (indep)")
<<<<<<< HEAD
  }else{
    colnames(table)<-c("Covariate",paste("n=",N,")",sep=""))
    
  }
=======
>>>>>>> 069d7ed3d5af49f9b8d8fa8f020db5886132e74e
  return(table)
}

pcovsum<-function(data,covs,maincov,numobs,latex=F){
  if(!latex){
    print.xtable(xtable(covsum(data,covs,maincov,numobs)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{  
    
    print.xtable(xtable(covsum(data,covs,maincov,numobs)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }}

uvsum<-function(data,response,covs,type,boxcox=F,strata=1,markup=T,sanitize=T,nicenames=T,testing=F){

  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity    
  }
  if(!sanitize) sanitizestr<-identity  
  if(!nicenames) nicename<-identity
  
  if(class(strata)!="numeric") {strata<-sapply(strata,function(stra){paste("strata(",stra,")",sep="")})
   }else{strata<-""}
  
  if(type=="coxph"){                                 
  beta<-"HR(95%CI)"
  }else if (type=="logistic"){
  beta<-"OR(95%CI)"  
  }else if (type=="linear"){
    beta<-"Estimate(95%CI)"
  }else{
    return("type must be either coxph, logisitc, linear")
  }
  out<-lapply(covs,function(cov){
    cov2<-cov
    if(testing) print(cov)
    if(is.factor(data[,cov])){
      levelnames<-sapply(sapply(sapply(levels(factor(data[,cov])),nicename),sanitizestr),addspace)
      cov<-lbld(sanitizestr(nicename(cov)))
      title<-NULL
      body<-NULL
      if(type=="coxph"){
      m2<-coxph(as.formula(paste(paste("Surv(",response[1],",",response[2],")",sep=""),"~",cov2,ifelse(strata=="","","+"),paste(strata,collapse="+"),sep="")),data=data)  
        
      hazardratio<-c("Reference",apply(matrix(summary(m2)$conf.int[,c(1,3,4)],ncol=3),1,psthr))    
      pvalue<-c("",sapply(summary(m2)$coef[,5],lpvalue))
      title<-c(cov,"","",lpvalue(summary(m2)$waldtest[3]))
      }else if(type=="logistic"){
        m2<-glm(as.formula(paste(response,"~",cov2,sep="")),family="binomial",data=data)   
        #globalpvalue<-1-pchisq(2*(summary(m2)$null.deviance-summary(m2)$deviance),summary(m2)$df.null-summary(m2)$df.residual)
        globalpvalue<-try(wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        
        m<-summary(m2)$coefficients
        hazardratio<-c("Reference",apply(cbind(exp(m[-1,1]),exp(m[-1,1]-1.96*m[-1,2]),exp(m[-1,1]+1.96*m[-1,2])),1,psthr))        
        pvalue<-c("",sapply(m[-1,4],lpvalue))
        title<-c(cov,"","",lpvalue(globalpvalue))
        
      }else if(type=="linear"){
        if(!boxcox){m2<-lm(as.formula(paste(response,"~",cov2,sep="")),data=data)
        }else{m2<-boxcoxlm(data[,response],data[,cov2])[[1]]}
        m<-summary(m2)$coefficients
        #globalpvalue<-anova(m2)[5][[1]][1])
        globalpvalue<-try(wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        
        
        hazardratio<-c("Reference",apply(cbind(m[-1,1],m[-1,1]-1.96*m[-1,2],m[-1,1]+1.96*m[-1,2]),1,psthr))        
        pvalue<-c("",sapply(m[-1,4],lpvalue))
        title<-c(cov,"","",lpvalue(globalpvalue))            
    }
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
      if(type=="coxph"){
        m2<-coxph(as.formula(paste(paste("Surv(",response[1],",",response[2],")",sep=""),"~",cov2,ifelse(strata=="","","+"),paste(strata,collapse="+"),sep="")),data=data)  
        
        out<-matrix(c(cov,psthr(summary(m2)$conf.int[,c(1,3,4)]),"",lpvalue(summary(m2)$waldtest[3])),ncol=4)
        
      }else if(type=="logistic"){
        m2<-glm(as.formula(paste(response,"~",cov2,sep="")),family="binomial",data=data)    
        
        m<-summary(m2)$coefficients
        
        #globalpvalue<-1-pchisq(2*(summary(m2)$null.deviance-summary(m2)$deviance),summary(m2)$df.null-summary(m2)$df.residual)
        globalpvalue<-try(wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        
        out<-matrix(c(cov,psthr(c(exp(m[-1,1]),exp(m[-1,1]-1.96*m[-1,2]),exp(m[-1,1]+1.96*m[-1,2]))),"",lpvalue(globalpvalue)),ncol=4)
        
        
      }else if(type=="linear"){
        if(!boxcox){m2<-lm(as.formula(paste(response,"~",cov2,sep="")),data=data)
        }else{m2<-boxcoxlm(data[,response],data[,cov2])[[1]]}
        #globalpvalue<-anova(m2)[5][[1]][1])
        globalpvalue<-try(wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        
        m<-summary(m2)$coefficients
        
        out<-matrix(c(cov,psthr(c(m[-1,1],m[-1,1]-1.96*m[-1,2],m[-1,1]+1.96*m[-1,2])),"",lpvalue(globalpvalue)),ncol=4)
                        
      }
      return(list(out,nrow(out)))}})
  table<-lapply(out,function(x){return(x[[1]])})
  table<-do.call(rbind.data.frame, table)  
  colnames(table)<-sapply(c("Covariate",sanitizestr(beta),"p-value","Global p-value"),lbld)
  return(table)
  
  index<-unlist(lapply(out,function(x){return(x[[2]])}))  
  lineindex<-rep(-1,length(index))
  for(i in 1:length(index)){lineindex[i]<-sum(index[c(1:i)])}
  lineindex<-c(-1,0,lineindex)
  if(output=="latex") return(list(table,lineindex))
  else return(table)
}

puvsum<-function(data,response,covs,type,boxcox=F,strata=1,latex=F){
  if(!latex){
    print.xtable(xtable(uvsum(data,response,covs,type,boxcox,strata)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{
    print.xtable(xtable(uvsum(data,response,covs,type,boxcox,strata)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }
  
}

mvsum<-function(data,model,type,markup=T,sanitize=T,nicenames=T){
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
  if(type=="linear"){
    betanames<-betanames[-1]
    beta<-"Estimate(95%CI)"    
  }else if(type=="logistic"){
    betanames<-betanames[-1]
    beta<-"OR(95%CI)"  
  }else if (type=="coxph"){
      beta<-"HR(95%CI)"
  }else{
      return("type must be linear, logistic, or coxph")
    }
  
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
  if(type%in%c("linear","logistic")){
    y<-lapply(y,function(x){ x+1})
    betanames<-c("intercept",betanames)
  }
  
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
    reference=NULL
    title=NULL
    body=NULL
    
    globalpvalue<-try(wald.test(b=coef(model)[covariateindex],Sigma=vcov(model)[covariateindex,covariateindex],Terms=seq_along(covariateindex))$result$chi2[3]);
    if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
    globalpvalue<-lpvalue(globalpvalue)
    
    if(type=="coxph"){
    hazardratio<-c(apply(matrix(summary(model)$conf.int[covariateindex,c(1,3,4)],ncol=3),1,psthr))
    pvalues<-c(sapply(summary(model)$coef[covariateindex,5],lpvalue))
    }else if (type=="logistic"){
      m<-summary(model)$coefficients      
      hazardratio<-apply(cbind(exp(m[covariateindex,1]),exp(m[covariateindex,1]-1.96*m[covariateindex,2]),
                                             exp(m[covariateindex,1]+1.96*m[covariateindex,2])),1,psthr)        
      pvalues<-c(sapply(m[covariateindex,4],lpvalue))
    }else if (type=="linear"){
      m<-summary(model)$coefficients      
      
      hazardratio<-apply(cbind(m[covariateindex,1],m[covariateindex,1]-1.96*m[covariateindex,2],
                               m[covariateindex,1]+1.96*m[covariateindex,2]),1,psthr)        
      pvalues<-sapply(m[covariateindex,4],lpvalue)
    }
    
    #if not interaction
    
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
  colnames(table)<-sapply(c("Covariate",sanitizestr(beta),"p-value","Global p-value"),lbld)
  return(table)
  lineindex<-rep(-1,length(index))
  #for(i in 1:length(index)){lineindex[i]<-sum(index[c(1:i)])}
  #lineindex<-c(-1,0,lineindex)
  if(output=="latex") return(list(table,index))
  else return(table)  
}

pmvsum<-function(data,model,type){
  print.xtable(xtable(mvsum(data,model,type)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
}

makedocx<-function(dir,fname,pdwd,imwd=""){
  oldwd<-getwd()
  if(imwd!=""){
  setwd(imwd)
  command<-paste("mogrify -path ", dir,"figure\\ ", "-format png ", dir, "figure\\*.pdf",sep="" )
  shell(command)
  }
  setwd(pdwd)
  command<-paste("pandoc -o ",dir,fname,".docx ",dir,fname,".tex ",
                 "--default-image-extension=png",sep="")
  shell(command)
  setwd(oldwd)
}

bwselect<-function(data,response,covs,strata=1,type,force=NULL,p=0.05,test=F,boxcox=F){
  while(T){
    
    if(type=="logistic"){
      pvalues<-sapply(covs,function(cov){
        m0<-glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
                data=subset(data,!is.na(data[,cov])),family="binomial")
        m1<-update(m0,as.formula(paste(".~.-",cov,sep="")))
        1-pchisq(summary(m1)$deviance-summary(m0)$deviance,
                 summary(m1)$df.residual-summary(m0)$df.residual)
      })          
      if(test) print(pvalues)
      if(length(covs)==1){
        if(max(pvalues)>p){
          return(glm(as.formula(paste(response,"~1",sep="")),
                     data=data,family="binomial"))
        }else{
          return(glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
                     data=data,family="binomial"))
        }}else{
          if(max(pvalues)>p){
            covs<-covs[-which.max(pvalues)]
          }else{
            return(glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
                       data=data,family="binomial"))
          }}
    }else if(type=="linear"){
      pvalues<-sapply(covs,function(cov){
        if(!boxcox){
          m0<-lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
                 data=subset(data,!is.na(data[,cov])))
          m1<-update(m0,as.formula(paste(".~.-",cov,sep="")))
        }else{
          m0<-boxcoxlm(data[,response],data[,covs])[[1]]
          m1<-boxcoxlm(data[!is.na(data[,cov]),response],data[!is.na(data[,cov]),setdiff(covs,cov)])[[1]]
        }
        1-pchisq(2*(logLik(m0)-logLik(m1)),
                 summary(m0)$df[1]-summary(m1)$df[1])
      })          
      if(test) print(pvalues)
      if(length(covs)==1){
        if(max(pvalues)>p){
          return(lm(as.formula(paste(response,"~1",sep="")),
                    data=data))
        }else{
          return(lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
                    data=data))
        }}else if(length(covs)==2 & boxcox & max(pvalues)>p){
          return(boxcoxlm(data[,response],data[,covs[which.min(pvalues)]]))
          
        }
      else{
        if(max(pvalues)>p){
          covs<-covs[-which.max(pvalues)]
        }else{
          if(!boxcox){
            return(lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
                      data=data))
          }else{
            return(boxcoxlm(data[,response],data[,covs]))
          }
        }}}
    
  }}

boxcoxlm<-function(y,x){
  require(geoR)
  missing<-unique(unlist(lapply(data.frame(y,x),function(xx) which(is.na(xx)))))
  notmissing<-setdiff(seq_len(length(y)),missing)
  a<-y[notmissing]
  if(class(x)=="data.frame"){
    b<-x[notmissing,]
    b<-sapply(b,function(bb){
      if(is.factor(bb)){
        bb<-as.numeric(bb)-1  
        if(length(table(bb))>=3){
          bb<-as.matrix(model.matrix(~factor(bb)-1)[,-1])
        }}
      return(bb)})
    if(class(b)=="list"){
      b<-as.matrix(do.call(cbind.data.frame, b))
    }else{
      b<-as.matrix(b)
    }}
  else{
    b<-x[notmissing]
    if (is.factor(b)) {
      b <- as.numeric(b) - 1
      if (length(table(x)) >= 3) {
        b <- as.matrix(model.matrix(~factor(b) - 1)[, -1])
      }
    }
    b<-as.matrix(b)   
  }
  
  
  
  
  lambda<-unlist(unlist(boxcoxfit(a,b,lambda2=T))[1:2])                     
  a<-((a+lambda[2])^lambda[1]-1)/lambda[1]
  out<-lm(a~b)
  return(list(out,lambda))
}

plotci<-function(data,response,units="months"){
  fita<-cuminc(data[,response[1]],data[,response[2]])
  plot(fita[[1]]$time,sapply(fita[[1]]$est+1.96*sqrt(fita[[1]]$var),function(x) min(x,1)),
       type="l",lty=2,main=paste("CI plot for ",sanitizestr(nicename(response[2])),sep=""),
       xlab=paste("Time (",cap(units),")",sep=""), ylim=c(0,1),
       ylab=paste("Incidence of ",sanitizestr(nicename(response[2])),sep=""))
  lines(fita[[1]]$time,fita[[1]]$est)
  lines(fita[[1]]$time,sapply(fita[[1]]$est-1.96*sqrt(fita[[1]]$var),
                              function(x)max(x,0)),lty=2)
}
