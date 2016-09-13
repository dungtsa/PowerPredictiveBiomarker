#' Power calculation of predictive biomarker. Calculation requires a series of parameters to to determine subgroup proportion and subgroup censoring rate. Depending on study type (prospective or retrospective study), parameter setting is different (check the reference).
#'
#' @description Power calcuation of predictive predictive biomarker in prospective and retrospective studies.
#' @param n.list sample size
#' @param mst.data input of median survival time (MST) for 4 subgroups:
#' \itemize {
#' \item low (negative) biomarker in control group;
#' \item high (positive) biomarker in control group;
#' \item low (negative) biomarker in treatment group;
#' \item high (positive) biomarker in treatment group
#' }
#' @param my.p.treat prevalence in treatment group
#' @param my.p.pos.control prevalence of positive biomarker in control group
#' @param my.p.pos.treatment prevalence of positive biomarker in treatment group
#' @param my.target.event.rate event rate (For MST, it would be 0.5)
#' @param my.follow.time follow time for prospective study
#' @param my.total.time total study time for prospective study
#' @param my.mst.overall Overall MST for a retrospective study cohort if avaialble
#' @param my.mst.treatment MST in treatment group for a retrospective study cohort if avaialble
#' @param my.mst.control MST in control group for a retrospective study cohort if avaialble
#' @param my.target.event.rate.overall Overall event rate in a retrospective study cohort if available
#' @param my.target.event.rate.treatment  event rate of treatment group in a retrospective study cohort if available
#' @param my.target.event.rate.control  event rate of control group in a retrospective study cohort if available
#' @param alpha type I error (two-sided)
#' @param study.design Study design: 'prospective' for prospective study and 'retrospective' for retorpsective study
#' @param my.censoringAndMST 'Overall' for available overall MSt and event rate; otherwise for availabel MST and event rate in both tretament and control groups
#' @return statistical power
#' @details This statistcial tool is to caculate power of predictive biomarker in survival data based on statsitical interaction model in Cox PH using Peterson et al. and Lachin methods.
#' @references  Strategies for power calculations in predictive biomarker studies in survival data; Chen et al; Oncotarget; 2016; in press
#' @examples
#' #-----------------Example: Prospective Study ------------------------------
#' # A prospective study with a sample size of 300
#' n.tmp=300
#' # a two sided 5% type I error
#' alpha.tmp=0.05
#' # Median survival time for the 4 subgroups
#' MSTdata.tmp=data.frame(C_low=3,C_high=1,T_low=4,T_high=3)
#' # 50% patients in the treatment group
#' p.treat.tmp=0.5
#' # In the treatment group, assume 50% patients with positive biomarker
#' p.pos.treat.tmp=0.5
#' # In the control group, assume 50% patients with positive biomarker
#' p.pos.control.tmp=0.5
#' # 1 years of follow up time
#' follow.time.tmp=1
#' # A total study time: 5 years
#' total.time.tmp=5
#' predictive.power.fun(n.list=n.tmp,mst.data=MSTdata.tmp,my.p.treat=p.treat.tmp,my.p.pos.control=p.pos.control.tmp,my.p.pos.treat=p.pos.treat.tmp,my.follow.time=follow.time.tmp,my.total.time=total.time.tmp,alpha=alpha.tmp,study.design='prospective')
#'
#' #-----------------Example: Retrospective Study ------------------------------
#' # A retrospective study with a sample size of 135
#' n.tmp=135
#' # a two sided 5% type I error
#' alpha.tmp=0.05
#' # Median survival time for the 4 subgroups
#' MSTdata.tmp=data.frame(C_low=10.11,C_high=3.1,T_low=6.66,T_high=11.01)
#' # 48% patients in the treatment group
#' p.treat.tmp=0.48
#' # In the treatment group, assume 50% patients with positive biomarker
#' p.pos.treat.tmp=0.5
#' # In the control group, assume 50% patients with positive biomarker
#' p.pos.control.tmp=0.5
#' # Available MST and event rate for both treatment and control in study cohort
#' censoringAndMST.tmp='TreatmentAndControl'
#' # MST (7.8 years) and event rate (64% censoring rate) for treatment
#' MST.treatment.study.cohort.tmp=7.8
#' event.rate.treatment.study.cohort.tmp=1-0.64
#' # MST (4.8 years) and event rate (55% censoring rate) for control
#' MST.control.study.cohort.tmp=4.8
#' event.rate.control.study.cohort.tmp=1-0.55
#' predictive.power.fun(n.list=n.tmp,mst.data=MSTdata.tmp,my.p.treat=p.treat.tmp,my.p.pos.control=p.pos.control.tmp,my.p.pos.treat=p.pos.treat.tmp,alpha=alpha.tmp,my.mst.treatment=MST.treatment.study.cohort.tmp,my.mst.control=MST.control.study.cohort.tmp,my.target.event.rate.treatment=event.rate.treatment.study.cohort.tmp,my.target.event.rate.control=event.rate.control.study.cohort.tmp,study.design='retrospective',my.censoringAndMST=censoringAndMST.tmp)
#'
#' @export


predictive.power.fun<-function(n.list=c(300),mst.data=data.frame(C_low=c(6.2,6.2),C_high=c(3.1,6.2),T_low=c(6.2,6.2),T_high=c(7,6.2)),my.p.treat=0.5,my.p.pos.control=.5,my.p.pos.treat=0.5,my.target.event.rate=0.5,my.follow.time=2,my.total.time=5,my.mst.overall=NULL,my.mst.treatment=NULL,my.mst.control=NULL,my.target.event.rate.overall=NULL,my.target.event.rate.treatment=NULL,my.target.event.rate.control=NULL,alpha=0.05,study.design='prospective',my.censoringAndMST='Overall')
{
  require(shiny)
  require(knitr)
  require(rmarkdown)
  require(survival)

#mst.data<-data.frame(C_low=c(6.2,6.2),C_high=c(3.1,6.2),T_low=c(6.2,6.2),T_high=c(7,6.2))
#dimnames(mst.data)[[1]]<-c('H1','H0')
dimnames(mst.data)[[1]][1]<-c('H1')

MST.cutoff.fun<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,target.event.rate=.5,shape.par=1)
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker

  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))
  t.tmp<-as.numeric(as.vector(data1[1,]))
  t1.seq<-as.numeric(seq(min(t.tmp),max(t.tmp)*3,by=0.01))
#  if(my.dist.tmp=='exp')
#  {

    event.rate.list.overall<-(log(2)/t1.seq)*((w1[1]/((log(2)/t1.seq)+log(2)/data1[1,1]))+(w1[2]/((log(2)/t1.seq)+log(2)/data1[1,2]))+(w1[3]/((log(2)/t1.seq)+log(2)/data1[1,3]))+(w1[4]/((log(2)/t1.seq)+log(2)/data1[1,4])))


    event.rate.list.treatment<-
      (log(2)/t1.seq)*(((w1[3]/sum(w1[3:4]))/((log(2)/t1.seq)+log(2)/data1[1,3]))+((w1[4]/sum(w1[3:4]))/((log(2)/t1.seq)+log(2)/data1[1,4])))
    event.rate.list.control<-
      (log(2)/t1.seq)*(((w1[1]/sum(w1[1:2]))/((log(2)/t1.seq)+log(2)/data1[1,1]))+((w1[2]/sum(w1[1:2]))/((log(2)/t1.seq)+log(2)/data1[1,2])))
    event.rate.list.pos<-
      (log(2)/t1.seq)*(((w1[2]/sum(w1[c(2,4)]))/((log(2)/t1.seq)+log(2)/data1[1,2]))+((w1[4]/sum(w1[c(2,4)]))/((log(2)/t1.seq)+log(2)/data1[1,4])))
    event.rate.list.neg<-
      (log(2)/t1.seq)*(((w1[1]/sum(w1[c(1,3)]))/((log(2)/t1.seq)+log(2)/data1[1,1]))+((w1[3]/sum(w1[c(1,3)]))/((log(2)/t1.seq)+log(2)/data1[1,3])))



#    }

  MST.overall.cutoff<-t1.seq[abs(event.rate.list.overall-target.event.rate)==min(abs(event.rate.list.overall-target.event.rate))]
  MST.treatment.cutoff<-t1.seq[abs(event.rate.list.treatment-target.event.rate)==min(abs(event.rate.list.treatment-target.event.rate))]
  MST.control.cutoff<-t1.seq[abs(event.rate.list.control-target.event.rate)==min(abs(event.rate.list.control-target.event.rate))]
  MST.pos.cutoff<-t1.seq[abs(event.rate.list.pos-target.event.rate)==min(abs(event.rate.list.pos-target.event.rate))]
  MST.neg.cutoff<-t1.seq[abs(event.rate.list.neg-target.event.rate)==min(abs(event.rate.list.neg-target.event.rate))]

  event.rate.data.combined<-cbind(event.rate.list.overall, event.rate.list.treatment, event.rate.list.control, event.rate.list.pos,event.rate.list.neg)
  MST.cutoff.list<-c(MST.overall.cutoff,MST.treatment.cutoff, MST.control.cutoff, MST.pos.cutoff, MST.neg.cutoff)
  list(MST.overall.cutoff=MST.overall.cutoff,
       MST.treatment.cutoff=MST.treatment.cutoff,
       MST.control.cutoff=MST.control.cutoff,
       MST.positive.cutoff=MST.pos.cutoff,
       MST.negative.cutoff=MST.neg.cutoff)
}

censor.rate.fun<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,target.event.rate=.5,shape.par=1,sim.nn=10000000,my.follow.time=NULL,my.total.time=NULL)
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker
  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))

  time.truncated.list<-MST.cutoff.fun(data1=data1,p.treat=p.treat,p.pos.control=p.pos.control,p.pos.treat=p.pos.treat,target.event.rate=target.event.rate,shape.par=shape.par)

  time.truncated<-time.truncated.list$MST.overall.cutoff

  event.rate.list<-list()
  for(j in 1:dim(data1)[1])
  {
  t.tmp<-as.numeric(as.vector(data1[j,]))
  rate.tmp<-rate.tmp0<-numeric()

#  if(my.dist.tmp=='exp')
#  {
    for(i in 1:length(t.tmp))
    {

        rate.tmp<-c(rate.tmp, 1-1/((my.total.time-my.follow.time)*(log(2)/t.tmp[i]))*(exp(-my.follow.time*log(2)/t.tmp[i])-exp(-my.total.time*log(2)/t.tmp[i])))
    }

#  }

    event.rate.list[[j]]<-list(event.rate.ind=rate.tmp,event.rate.overall=sum(rate.tmp*w1))
  }
#  names(event.rate.list)<-c('H1','H0')
  names(event.rate.list)<-c('H1')
  event.rate.list<-list(event.rate=event.rate.list,prop.n=w1,MST.list=time.truncated.list)
  event.rate.list
}

censor.fun<-function(x,time.cutoff=10,interval.follow=1,total.freq=5)
{
  n<-dim(x)[1]
  n.k<-round(n/total.freq)
  data.tmp1<-numeric(0)
  for(i in 1:total.freq)
  {
    if(i!=total.freq) index1<-(1:n.k)+(i-1)*n.k else
      index1<-((i-1)*n.k+1):n
    data.tmp<-x[index1,]
    time.cutoff.tmp<-time.cutoff+interval.follow*(i-1)
    data.tmp$censor<- as.numeric(data.tmp$time<=time.cutoff.tmp)
    data.tmp$time[data.tmp$censor==0]<-time.cutoff.tmp
    data.tmp1<-rbind(data.tmp1,data.tmp)
  }
  data.tmp1
}

#--for retrospective study---
MST.cutoff.fun.scale.factor<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,target.event.rate=.5,shape.par=1,mst.overall=10,mst.treatment=15,mst.control=8)
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker
  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))
  t.tmp<-as.numeric(as.vector(data1[1,]))
  k.overall<-k.treatment<-k.control<-seq(0.01,10,by=0.01)
  t1.seq<-as.numeric(seq(min(t.tmp),max(t.tmp)*3,by=0.01))
#  if(my.dist.tmp=='exp')
#  {

    if(!is.null(mst.overall))   event.rate.list.overall<-(log(2)/mst.overall)*((w1[1]/((log(2)/mst.overall)+log(2)/(k.overall*data1[1,1])))+(w1[2]/((log(2)/mst.overall)+log(2)/(k.overall*data1[1,2])))+(w1[3]/((log(2)/mst.overall)+log(2)/(k.overall*data1[1,3])))+(w1[4]/((log(2)/mst.overall)+log(2)/(k.overall*data1[1,4])))) else event.rate.list.overall<-NULL

    if(!is.null(mst.treatment))   event.rate.list.treatment<-(log(2)/mst.treatment)*(((w1[3]/sum(w1[3:4]))/((log(2)/mst.treatment)+log(2)/(k.treatment*data1[1,3])))+((w1[4]/sum(w1[3:4]))/((log(2)/mst.treatment)+log(2)/(k.treatment*data1[1,4])))) else event.rate.list.treatment<-NULL
    if(!is.null(mst.control)) event.rate.list.control<-
        (log(2)/mst.control)*(((w1[1]/sum(w1[1:2]))/((log(2)/mst.control)+log(2)/(k.control*data1[1,1])))+((w1[2]/sum(w1[1:2]))/((log(2)/mst.control)+log(2)/(k.control*data1[1,2])))) else event.rate.list.control<-NULL


#  }


  if(!is.null(event.rate.list.overall))
  {
    index.overall<-abs(event.rate.list.overall-target.event.rate)==min(abs(event.rate.list.overall-target.event.rate))
    overall<-cbind(k.overall,event.rate.list.overall)[index.overall,]
  } else overall<-NULL
  if(!is.null(event.rate.list.treatment))
  {
    index.treatment<-abs(event.rate.list.treatment-target.event.rate)==min(abs(event.rate.list.treatment-target.event.rate))
    treatment<-cbind(k.treatment,event.rate.list.treatment)[index.treatment,]
  } else treatment<-NULL
  if(!is.null(event.rate.list.control))
  {
    index.control<-abs(event.rate.list.control-target.event.rate)==min(abs(event.rate.list.control-target.event.rate))
    control<-cbind(k.control,event.rate.list.control)[index.control,]
  } else control<-NULL
  list(overall=overall,treatment=treatment,control=control)

}

censoring.cutoff.fun.retro<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,target.event.rate.overall=NULL,target.event.rate.treatment=NULL,target.event.rate.control=NULL, shape.par=1)
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker

  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))
  t.tmp<-as.numeric(as.vector(data1[1,]))
  t1.seq<-as.numeric(seq(0.001,max(t.tmp)*3,by=0.001))
#  if(my.dist.tmp=='exp')
#  {
    #-------
    #-- prob(A>B)= r_A/(r_A+r_B) or (log2/MST_A)/(log2/MST_A+log2/MST_B)
    #-so event rate will be sum( (log2/MST_targeted)/(log2/MST_targeted+log2/MST_subgroup))

    if(!is.null(target.event.rate.overall))    event.rate.list.overall<-(log(2)/t1.seq)*((w1[1]/((log(2)/t1.seq)+log(2)/data1[1,1]))+(w1[2]/((log(2)/t1.seq)+log(2)/data1[1,2]))+(w1[3]/((log(2)/t1.seq)+log(2)/data1[1,3]))+(w1[4]/((log(2)/t1.seq)+log(2)/data1[1,4]))) else event.rate.list.overall<-NULL
    if(!is.null(target.event.rate.treatment))    event.rate.list.treatment<-
        (log(2)/t1.seq)*(((w1[3]/sum(w1[3:4]))/((log(2)/t1.seq)+log(2)/data1[1,3]))+((w1[4]/sum(w1[3:4]))/((log(2)/t1.seq)+log(2)/data1[1,4]))) else event.rate.list.treatment<-NULL
        if(!is.null(target.event.rate.control))    event.rate.list.control<-
            (log(2)/t1.seq)*(((w1[1]/sum(w1[1:2]))/((log(2)/t1.seq)+log(2)/data1[1,1]))+((w1[2]/sum(w1[1:2]))/((log(2)/t1.seq)+log(2)/data1[1,2]))) else event.rate.list.control<-NULL
#  }


  if(!is.null(event.rate.list.overall)) MST.overall.cutoff<-t1.seq[abs(event.rate.list.overall-target.event.rate.overall)==min(abs(event.rate.list.overall-target.event.rate.overall))] else MST.overall.cutoff<-NULL
  if(!is.null(event.rate.list.treatment))   MST.treatment.cutoff<-t1.seq[abs(event.rate.list.treatment-target.event.rate.treatment)==min(abs(event.rate.list.treatment-target.event.rate.treatment))] else MST.treatment.cutoff<-NULL
  if(!is.null(event.rate.list.control))     MST.control.cutoff<-t1.seq[abs(event.rate.list.control-target.event.rate.control)==min(abs(event.rate.list.control-target.event.rate.control))] else MST.control.cutoff<-NULL


  list(MST.overall.cutoff=MST.overall.cutoff,
       MST.treatment.cutoff=MST.treatment.cutoff,
       MST.control.cutoff=MST.control.cutoff)
}


censor.rate.fun.retro<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,shape.par=1,sim.nn=10000000,my.follow.time=NULL,my.total.time=NULL,mst.overall=10,mst.treatment=15,mst.control=8,target.event.rate.overall=NULL,target.event.rate.treatment=NULL,target.event.rate.control=NULL,censoringAndMST='Overall')
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker
  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))

  MST.scale<-MST.cutoff.fun.scale.factor(data1=data1,p.treat=p.treat,p.pos.control=p.pos.control,p.pos.treat=p.pos.treat,target.event.rate=0.5,shape.par=shape.par,mst.overall=mst.overall,mst.treatment=mst.treatment,mst.control=mst.control)

  data1.scale<-data1
  if(censoringAndMST=='Overall') data1.scale<-data1*MST.scale$overall[1] else
  {
    data1.scale[1,1:2]<-data1.scale[1,1:2]*MST.scale$control[1]
    data1.scale[1,3:4]<-data1.scale[1,3:4]*MST.scale$treatment[1]
  }

  MST.list.org<-MST.cutoff.fun(data1=data1,p.treat=p.treat,p.pos.control=p.pos.control,p.pos.treat=p.pos.treat,target.event.rate=0.5,shape.par=shape.par)

  MST.list.scale<-MST.cutoff.fun(data1=data1.scale,p.treat=p.treat,p.pos.control=p.pos.control,p.pos.treat=p.pos.treat,target.event.rate=0.5,shape.par=shape.par)

  MST.censoring.list<-censoring.cutoff.fun.retro(data1=data1.scale,p.treat=p.treat,p.pos.control=p.pos.control, target.event.rate.overall=target.event.rate.overall,target.event.rate.treatment=target.event.rate.treatment,target.event.rate.control=target.event.rate.control, shape.par=shape.par)

  if(censoringAndMST=='Overall') time.truncated<-MST.censoring.list$MST.overall.cutoff else time.truncated<-unlist(MST.censoring.list[3:2])
  event.rate.list<-list()
  for(j in 1:1)
  {
    t.tmp<-as.numeric(as.vector(data1.scale[j,]))
    rate.tmp<-rate.tmp1<-numeric()
#    if(my.dist.tmp=='exp')
#    {
      if(censoringAndMST=='Overall')
      {

       rate.tmp<-(log(2)/time.truncated)/((log(2)/t.tmp)+(log(2)/time.truncated))

      } else
      {

        rate.tmp<-c((log(2)/time.truncated[1])/((log(2)/t.tmp[1:2])+(log(2)/time.truncated[1])),(log(2)/time.truncated[2])/((log(2)/t.tmp[3:4])+(log(2)/time.truncated[2])))

      }
#    }

    event.rate.list[[j]]<-list(event.rate.ind=rate.tmp,event.rate.overall=sum(rate.tmp*w1))
  }
  names(event.rate.list)<-c('H1')
  event.rate.list.all<-list(event.rate=event.rate.list,prop.n=w1,censoring.MST.list=time.truncated, scale.factor=MST.scale,data.scale=data1.scale,MST.data.scale=MST.list.scale,MST.data.org=MST.list.org)
  event.rate.list.all
}

#-------
if(study.design=='prospective') censor.rate.est<-censor.rate.fun(data1=mst.data,p.treat=my.p.treat,p.pos.control=my.p.pos.control,p.pos.treat=my.p.pos.treat,target.event.rate=my.target.event.rate,shape.par=shape.par1,my.follow.time=my.follow.time,my.total.time=my.total.time) else
  censor.rate.est<-censor.rate.fun.retro(data1=mst.data,p.treat=my.p.treat,p.pos.control=my.p.pos.control,p.pos.treat=my.p.pos.treat,shape.par=shape.par1,mst.overall=my.mst.overall,mst.treatment=my.mst.treatment,mst.control=my.mst.control,target.event.rate.overall=my.target.event.rate.overall,target.event.rate.treatment=my.target.event.rate.treatment,target.event.rate.control=my.target.event.rate.control,censoringAndMST=my.censoringAndMST)


sd1<-sqrt(sum(1/(censor.rate.est$event.rate$H1$event.rate.ind*censor.rate.est$prop.n*n.list[1])))

z1<-qnorm(1-alpha/2)
delta<-log(prod(mst.data[1,c(1,4)])/prod(mst.data[1,c(2,3)]))
power.approx<-1 - (pnorm(z1 - abs(delta)/sd1) - pnorm(-z1 - abs(delta)/sd1))
Peterson.power<-power.approx


prop.n<-censor.rate.est$prop.n
HR1<-prod(mst.data[1,c(2,3)])/prod(mst.data[1,c(1,4)])

intpower=function(Nsample,eventA,eventB,LHRA,LHRB,alpha){
  A = c(1/eventA[1]+1/eventA[2])
  B = c(1/eventB[1]+1/eventB[2])
  VB = solve(solve(B)+solve(A))
  bet = VB%*%(solve(A)%*%LHRA+solve(B)%*%LHRB)
  a1 = t(LHRA-bet)%*%solve(A*Nsample)%*%(LHRA-bet)
  a2 = t(LHRB-bet)%*%solve(B*Nsample)%*%(LHRB-bet)
  ps = (a1+a2)*Nsample
  power = 1-pchisq(qchisq(1-alpha, df=1),1,ps)
  return(power)
}


eventNum = censor.rate.est$event.rate$H1$event.rate.ind*censor.rate.est$prop.n*n.list[1]
eventA=c(eventNum[1],eventNum[2])
eventB=c(eventNum[3],eventNum[4])
LHRA=c(log(prod(mst.data[1,c(1,4)])))#neg
LHRB=c(log(prod(mst.data[1,c(2,3)])))#pos
Lachin = intpower(n.list[1],eventA,eventB,LHRA,LHRB,alpha)


#---------summary: end----


if(study.design=='prospective') PowerAns<-list(mst.data=mst.data,freq=prop.n,HRR=HR1,event.rate=censor.rate.est$event.rate$H1,Peterson.power=Peterson.power,Lachin.power=Lachin,MST.list=censor.rate.est$MST.list) else

  PowerAns<-list(mst.data=mst.data,freq=prop.n,HRR=HR1,event.rate=censor.rate.est$event.rate$H1,Peterson.power=Peterson.power,Lachin.power=Lachin,MST.list=censor.rate.est$MST.list,censoring.MST=censor.rate.est$censoring.MST.list, scale.factor=censor.rate.est$scale.factor,data.scale=censor.rate.est$data.scale,MST.data.scale=censor.rate.est$MST.data.scale,MST.data.org=censor.rate.est$MST.data.org)

 PowerAns$Lachin.power
}

