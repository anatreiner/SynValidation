########################################################################################################################
#                                                                                                                      #
#   Program for validating the results obtained from synthetic data:                                                   #
#   1. compare 5 synthetic sets to the real data:                                                                      #
#   2. compare the mean across estimates obtained from 1,000 synthetic sets, to the estimate obtained from real data   #
#                                                                                                                      #
#   Analysis - Survival:                                                                                               #
#   - Kaplan-Meier survival curves                                                                                     #  
#   - Cox proportional hazard regression                                                                               #
#                                                                                                                      #
########################################################################################################################


# define survival event
{
  outcome<-"Mortality" # name the type of survival event (e.g. "CHF", "Mortality_or_CHF")
  days.num<-180 # choose survival followup period
}

### assume 1000 files of synthetic data (named with prefix "Syn_" followed by the number) for a given file of real data ###

##########################
# analyze synthetic data #
##########################

library(readr)
km.data<-numeric()
hr.data<-numeric()
files.num<-1000 # number of repeatedly generated synthetic sets
for(k in 1:files.num)
{
  # read data:
  {  
    filename<-paste("Syn_",k,sep="")
    current.data<-read_csv(paste(filename,".csv",sep=""),col_names = TRUE,col_types = cols())
  }
  #  define variables for survival analysis:
  {
    # The following variables should be defined from the dataset:
    # status (1=event happened during the followup period, 0=event did not happen during the followup period)
    # surv.time (time to event)
    surv.obj <- with(current.data, Surv(surv.time, status == 1))
  }  
  # Kaplan-Meier curve 
  for (i in 1)
  {  
    km.SYN<-survfit(surv.obj~1,data=current.data,conf.type="log")
    km.summary<-summary(km.SYN)
    syn.num<-rep(k,length(km.summary$time))
    km.data<-rbind(km.data,data.frame(syn.num,km.summary$time,km.summary$surv,km.summary$lower,km.summary$upper))  # output 1 
    if (k<=5)
      assign(paste("km.SYN",k,sep=""),km.SYN)
  }
  # Cox proportional hazard regression model 
  for(i in 1)
  {  
    library(MASS)
    # coxfit<-coxph(Surv(surv.time,status==1)~#add names of explanatory variables#)
    sum.cox<-summary(coxfit)
    # output 2:
    hr.SYN<-sum.cox$conf.int[,1]
    hr.low.SYN<-sum.cox$conf.int[,3]
    hr.up.SYN<-sum.cox$conf.int[,4]
    syn.num<-rep(k,length(hr.SYN))
    # var.name<- c() # define a vector of characters for naming the explanatory variables
    hr.data<-rbind(hr.data,data.frame(syn.num,var.name,as.vector(hr.SYN),as.vector(hr.low.SYN),as.vector(hr.up.SYN)))
    # keep results for 5 synthetic datasets 
    if (k<=5)
    {
      assign(paste("hr.SYN",k,sep=""),hr.SYN)
      assign(paste("hr.low.SYN",k,sep=""),hr.low.SYN)
      assign(paste("hr.up.SYN",k,sep=""),hr.up.SYN)
    }
  }  
}

names(hr.data)<-c("syn.num","var.name","hr.SYN","hr.low.SYN","hr.up.SYN")
hr.MEAN.SYN<-tapply(hr.data$hr.SYN,hr.data$var.name,mean)
hr.min.SYN<-tapply(hr.data$hr.SYN,hr.data$var.name,min)
hr.max.SYN<-tapply(hr.data$hr.SYN,hr.data$var.name,max)
hr.sd.SYN<-tapply(hr.data$hr.SYN,hr.data$var.name,sd)
km.curve.MEAN.SYN<-tapply(km.data$km.summary.surv,km.data$km.summary.time,mean,na.rm=TRUE)

#####################
# analyze real data #
#####################

library(readr)
km.data<-numeric()
hr.data<-numeric()
{
  # read data:
  {
    filename<-"Real"
    current.data<-read_csv(paste(filename,".csv",sep=""),col_names = TRUE,col_types = cols())
  }
  #  variables for survival analysis:
  {
    # The following variables should be defined from the dataset:
    # status (1=event happened during the followup period, 0=event did not happen during the followup period)
    # surv.time (time from diagnosis/procedure to event)
    surv.obj <- with(current.data, Surv(surv.time, status == 1))
  }  
  # Kaplan-Meier curve (OUTPUT 1)
  for (i in 1)
  {  
    km.REAL<-survfit(surv.obj~1,data=current.data,conf.type="log")
    km.summary<-summary(km.REAL)
    km.data<-rbind(km.data,data.frame(km.summary$time,km.summary$surv,km.summary$lower,km.summary$upper))  # output 1 
  }
  # Cox proportional hazard regression model (OUTPUT 2)
  for(i in 1)
  {  
    library(MASS)
    # coxfit<-coxph(Surv(surv.time,status==1)~#add names of explanatory variables#)
    sum.cox<-summary(coxfit)
    # output 2:
    hr.REAL<-sum.cox$conf.int[,1]
    hr.low.REAL<-sum.cox$conf.int[,3]
    hr.up.REAL<-sum.cox$conf.int[,4]
    # var.name<- c() # define a vector of characters for naming the explanatory variables
    hr.data<-rbind(hr.data,data.frame(var.name,as.vector(hr.REAL),as.vector(hr.low.REAL),as.vector(hr.up.REAL)))
  }  
}  

########
# plot #
########

# compare hazard ratios:

{
  n.hr<-cbind(n.REAL.diff.hours.gt90," ",n.REAL.males," ",n.REAL.severe.card,n.REAL.prior.ischemic,n.REAL.high.bun,n.REAL.high.creat,n.REAL.low.hmg)
  pct.hr<-cbind(pct.REAL.diff.hours.gt90," ",pct.REAL.males," ",pct.REAL.severe.card,pct.REAL.prior.ischemic,pct.REAL.high.bun,pct.REAL.high.creat,pct.REAL.low.hmg)
  hr<-rbind(hr.REAL,hr.SYN1,hr.SYN2,hr.SYN3,hr.SYN4,hr.SYN5)
  hr.low<-rbind(hr.low.REAL,hr.low.SYN1,hr.low.SYN2,hr.low.SYN3,hr.low.SYN4,hr.low.SYN5)
  hr.up<-rbind(hr.up.REAL,hr.up.SYN1,hr.up.SYN2,hr.up.SYN3,hr.up.SYN4,hr.up.SYN5)
  windows(width=80,height=50)
  par(mfrow=c(1,1),mar=c(10,5,5,5))
  library(gplots)
  x<-barplot2(hr,beside=TRUE,ylab="Hazard Ratio",ylim=c(0,6),col=c("lightblue",rep("lightgreen",5)),las=3,plot.grid = FALSE,cex.names=1.2,cex.lab=1.5,cex.axis=1.5,cex.main=2,
              names.arg=var.name,
              main=paste("Hazard Ratio (by Cox Regression) - ",surv.type,sep=""),
              plot.ci = TRUE,ci.l = hr.low, ci.u = hr.up,ci.color="black",ci.lwd=1,ci.lty=2) 
  legend(50,5.8,c("real","synthetic","1000-synthetic mean"),col=c("lightblue","lightgreen","red"),lty=c(1,1,2),lwd=2,cex=1.3,bg="white")
}

# compare surival curves:

for (i in 1)
{  
  windows(width=50,height=30)
  par(mfrow=c(1,1),mar=c(5,5,7,5))
  plot(km.SYN1,mark.time = TRUE,conf.int = FALSE,col="green",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,lwd=2,lty=1,main=paste("Survival curve - ",surv.type,sep=""),xlab="Time (days)",ylab="survival (proportion)",ylim=c(y.lim.km.all,1))
  lines(km.SYN2,mark.time = TRUE,conf.int = FALSE,col="green",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,lwd=2,lty=1)
  lines(km.SYN3,mark.time = TRUE,conf.int = FALSE,col="green",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,lwd=2,lty=1)
  lines(km.SYN4,mark.time = TRUE,conf.int = FALSE,col="green",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,lwd=2,lty=1)
  lines(km.SYN5,mark.time = TRUE,conf.int = FALSE,col="green",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,lwd=2,lty=1)
  lines(km.REAL,mark.time = TRUE,conf.int = TRUE,col="blue",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,lwd=c(2,1,1))
  km.SYN.loess<-loess(km.curve.MEAN.SYN~as.numeric(names(km.curve.MEAN.SYN)),span=0.1)
  lines(km.SYN.loess$x,km.SYN.loess$fitted,col="red",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,lwd=2,lty=1)
  legend(120,1,c("real","real 95% conf. limits","synthetic","1000-synthetic mean (smoothed)"),col=c("blue","blue","green","red"),lty=c(1,2,1,1),lwd=2,cex=1.3)
  grid()
}



