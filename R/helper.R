pstprn<-function(x){paste(x[1]," (",paste(x[-1],collapse=","),")",sep="")}
psthr<-function(x,y=2){paste(round(x[1],y)," (",round(x[2],2),",",round(x[3],2),")",sep="")}
covnm<-function(betanames,call){
  sapply(betanames,function(betaname){
    
    indx<-which(sapply(call,function(cov){charmatch(cov,betaname)})==1)
    if(length(indx)==1) return(call[indx])
    #If one  facorname is a subset of another
    indx2<-which.max(sapply(call[indx],nchar))
    if(length(indx2)==1) return(call[indx[indx2]])
    indx3<-which(sapply(call[indx2],function(c){substr(betaname,1,nchar(c))==c}))
    if(length(indx3)==1)  return(call[indx[indx2[indx3]]])                      
  })  
}

alleql<-function(x,y){
  !any((x==y)==F)
}


betaindx<-function(x){
  i=1
  out<-1
  result<-NULL
  while(TRUE){
    if(i+1>length(x)){
      result<-c(result,list(out))
      return(result)
    }
    else if(alleql(x[[i+1]],x[[i]])){
      out<-c(out,i+1)
    }
    else{
      result<-c(result,list(out))
      out<-i+1
    }
    i=i+1
  }
}


cap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

nicename<-function(strings){
  out<-sapply(strings,function(x){
    x<-chartr(".", " ",x)
    x<-chartr("_", " ",x)
    return(x)})
  return(out)
}
pvalue<-function(x){
  if(is.na(x)|class(x)=="character") return(x)
  else if (x<=0.001) return("<0.001")
  else return(signif(x,2))
}

sanitize <- function(str) {
  result <- str
  result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
  result <- gsub("$", "\\$", result, fixed = TRUE)
  result <- gsub(">", "$>$", result, fixed = TRUE)
  result <- gsub("<", "$<$", result, fixed = TRUE)
  result <- gsub("|", "$|$", result, fixed = TRUE)
  result <- gsub("{", "\\{", result, fixed = TRUE)
  result <- gsub("}", "\\}", result, fixed = TRUE)
  result <- gsub("%", "\\%", result, fixed = TRUE)
  result <- gsub("&", "\\&", result, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result <- gsub("#", "\\#", result, fixed = TRUE)
  result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
  result <- gsub("~", "\\~{}", result, fixed = TRUE)
  result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", 
                 result, fixed = TRUE)
  return(result)
}

sanitizestr<-function(str){
  as.vector(sapply(str,function(char){sanitize(char)}))
}

lbld<-function(strings){sapply(strings,function(x){
  if(is.null(x)) return(x)
  if(is.na(x)) return(x)
  return(paste("\\textbf{",x,"}",sep=""))})}

addspace<-function(x){
  paste("~~~",x,sep="")
}


lpvalue<-function(x){
  if(is.na(x)|class(x)=="character") return(x)
  else if (x<=0.001) return("\\textbf{$<$0.001}")
  else x=signif(x,2)
  if(x<=0.05) return(paste("\\textbf{",x,"}",sep=""))
  else return(x)
}

boxcoxlm<-function(y,x){
  require(geoR)
  missing<-unique(unlist(lapply(data.frame(y,x),function(xx) which(is.na(xx)))))
  notmissing<-setdiff(seq_len(length(y)),missing)
  a<-y[notmissing]
  if(class(x)=="data.frame"){
    b<-x[notmissing,]
  }else{
    b<-x[notmissing]
  }
  b<-sapply(b,function(bb){
    if(is.factor(bb)){
      bb<-as.numeric(bb)-1  
      if(length(table(bb))>=3){
        bb<-as.matrix(model.matrix(~factor(bb)-1)[,-1])
      }}
    return(bb)}
  )
  if(class(b)=="list"){
    b<-as.matrix(do.call(cbind.data.frame, b))
  }else{
    b<-as.matrix(b)              
  }
  
  lambda<-unlist(unlist(boxcoxfit(a,b,lambda2=T))[1:2])                     
  a<-((a+lambda[2])^lambda[1]-1)/lambda[1]
  out<-lm(a~b)
  return(list(out,lambda))
}