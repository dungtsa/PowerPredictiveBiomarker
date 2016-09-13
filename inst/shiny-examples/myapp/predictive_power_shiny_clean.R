library(shiny)
library(knitr)
library(rmarkdown)
require(survival)

predictive.power.fun<-function(n.list=c(300),sim.n=10,mst.data=data.frame(C_low=c(6.2,6.2),C_high=c(3.1,6.2),T_low=c(6.2,6.2),T_high=c(7,6.2)),my.p.treat=0.5,my.p.pos.control=.5,my.p.pos.treat=0.5,my.target.event.rate=0.5,my.follow.time=NULL,my.total.time=5,shape.par1=1,my.dist='exp',my.mst.overall=NULL,my.mst.treatment=NULL,my.mst.control=NULL,my.target.event.rate.overall=NULL,my.target.event.rate.treatment=NULL,my.target.event.rate.control=NULL,alpha=0.05,study.design='prospective',my.censoringAndMST='Overall')
{

#mst.data<-data.frame(C_low=c(6.2,6.2),C_high=c(3.1,6.2),T_low=c(6.2,6.2),T_high=c(7,6.2))
dimnames(mst.data)[[1]]<-c('H1','H0')

MST.cutoff.fun<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,target.event.rate=.5,shape.par=1,my.dist.tmp='exp')
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker

  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))
  t.tmp<-as.numeric(as.vector(data1[1,]))
  t1.seq<-as.numeric(seq(min(t.tmp),max(t.tmp)*3,by=0.01))
  if(my.dist.tmp=='exp')
  {

    event.rate.list.overall<-(log(2)/t1.seq)*((w1[1]/((log(2)/t1.seq)+log(2)/data1[1,1]))+(w1[2]/((log(2)/t1.seq)+log(2)/data1[1,2]))+(w1[3]/((log(2)/t1.seq)+log(2)/data1[1,3]))+(w1[4]/((log(2)/t1.seq)+log(2)/data1[1,4])))


    event.rate.list.treatment<-
      (log(2)/t1.seq)*(((w1[3]/sum(w1[3:4]))/((log(2)/t1.seq)+log(2)/data1[1,3]))+((w1[4]/sum(w1[3:4]))/((log(2)/t1.seq)+log(2)/data1[1,4])))
    event.rate.list.control<-
      (log(2)/t1.seq)*(((w1[1]/sum(w1[1:2]))/((log(2)/t1.seq)+log(2)/data1[1,1]))+((w1[2]/sum(w1[1:2]))/((log(2)/t1.seq)+log(2)/data1[1,2])))
    event.rate.list.pos<-
      (log(2)/t1.seq)*(((w1[2]/sum(w1[c(2,4)]))/((log(2)/t1.seq)+log(2)/data1[1,2]))+((w1[4]/sum(w1[c(2,4)]))/((log(2)/t1.seq)+log(2)/data1[1,4])))
    event.rate.list.neg<-
      (log(2)/t1.seq)*(((w1[1]/sum(w1[c(1,3)]))/((log(2)/t1.seq)+log(2)/data1[1,1]))+((w1[3]/sum(w1[c(1,3)]))/((log(2)/t1.seq)+log(2)/data1[1,3])))



    }

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

censor.rate.fun<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,target.event.rate=.5,shape.par=1,my.dist.tmp='exp',sim.nn=10000000,my.follow.time=NULL,my.total.time=NULL)
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker
  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))

  time.truncated.list<-MST.cutoff.fun(data1=data1,p.treat=p.treat,p.pos.control=p.pos.control,p.pos.treat=p.pos.treat,target.event.rate=target.event.rate,shape.par=shape.par,my.dist.tmp=my.dist.tmp)

  time.truncated<-time.truncated.list$MST.overall.cutoff

  event.rate.list<-list()
  for(j in 1:dim(data1)[1])
  {
  t.tmp<-as.numeric(as.vector(data1[j,]))
  rate.tmp<-rate.tmp0<-numeric()

  if(my.dist.tmp=='exp')
  {
    for(i in 1:length(t.tmp))
    {

        rate.tmp<-c(rate.tmp, 1-1/((my.total.time-my.follow.time)*(log(2)/t.tmp[i]))*(exp(-my.follow.time*log(2)/t.tmp[i])-exp(-my.total.time*log(2)/t.tmp[i])))
    }

  }

    event.rate.list[[j]]<-list(event.rate.ind=rate.tmp,event.rate.overall=sum(rate.tmp*w1))
  }
  names(event.rate.list)<-c('H1','H0')
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
MST.cutoff.fun.scale.factor<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,target.event.rate=.5,shape.par=1,my.dist.tmp='exp',mst.overall=10,mst.treatment=15,mst.control=8)
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker

  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))
  t.tmp<-as.numeric(as.vector(data1[1,]))
  k.overall<-k.treatment<-k.control<-seq(0.01,10,by=0.01)
  t1.seq<-as.numeric(seq(min(t.tmp),max(t.tmp)*3,by=0.01))
  if(my.dist.tmp=='exp')
  {

    if(!is.null(mst.overall))   event.rate.list.overall<-(log(2)/mst.overall)*((w1[1]/((log(2)/mst.overall)+log(2)/(k.overall*data1[1,1])))+(w1[2]/((log(2)/mst.overall)+log(2)/(k.overall*data1[1,2])))+(w1[3]/((log(2)/mst.overall)+log(2)/(k.overall*data1[1,3])))+(w1[4]/((log(2)/mst.overall)+log(2)/(k.overall*data1[1,4])))) else event.rate.list.overall<-NULL

    if(!is.null(mst.treatment))   event.rate.list.treatment<-(log(2)/mst.treatment)*(((w1[3]/sum(w1[3:4]))/((log(2)/mst.treatment)+log(2)/(k.treatment*data1[1,3])))+((w1[4]/sum(w1[3:4]))/((log(2)/mst.treatment)+log(2)/(k.treatment*data1[1,4])))) else event.rate.list.treatment<-NULL
    if(!is.null(mst.control)) event.rate.list.control<-
        (log(2)/mst.control)*(((w1[1]/sum(w1[1:2]))/((log(2)/mst.control)+log(2)/(k.control*data1[1,1])))+((w1[2]/sum(w1[1:2]))/((log(2)/mst.control)+log(2)/(k.control*data1[1,2])))) else event.rate.list.control<-NULL
        if(1>2)
        {
          event.rate.list.pos<-
            (log(2)/t1.seq)*(((w1[2]/sum(w1[c(2,4)]))/((log(2)/t1.seq)+log(2)/data1[1,2]))+((w1[4]/sum(w1[c(2,4)]))/((log(2)/t1.seq)+log(2)/data1[1,4])))
          event.rate.list.neg<-
            (log(2)/t1.seq)*(((w1[1]/sum(w1[c(1,3)]))/((log(2)/t1.seq)+log(2)/data1[1,1]))+((w1[3]/sum(w1[c(1,3)]))/((log(2)/t1.seq)+log(2)/data1[1,3])))
        }


  }


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

censoring.cutoff.fun.retro<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,target.event.rate.overall=NULL,target.event.rate.treatment=NULL,target.event.rate.control=NULL, shape.par=1,my.dist.tmp='exp')
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker

  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))
  t.tmp<-as.numeric(as.vector(data1[1,]))
  t1.seq<-as.numeric(seq(0.001,max(t.tmp)*3,by=0.001))
  if(my.dist.tmp=='exp')
  {
    #-------
    #-- prob(A>B)= r_A/(r_A+r_B) or (log2/MST_A)/(log2/MST_A+log2/MST_B)
    #-so event rate will be sum( (log2/MST_targeted)/(log2/MST_targeted+log2/MST_subgroup))

    if(!is.null(target.event.rate.overall))    event.rate.list.overall<-(log(2)/t1.seq)*((w1[1]/((log(2)/t1.seq)+log(2)/data1[1,1]))+(w1[2]/((log(2)/t1.seq)+log(2)/data1[1,2]))+(w1[3]/((log(2)/t1.seq)+log(2)/data1[1,3]))+(w1[4]/((log(2)/t1.seq)+log(2)/data1[1,4]))) else event.rate.list.overall<-NULL
    if(!is.null(target.event.rate.treatment))    event.rate.list.treatment<-
        (log(2)/t1.seq)*(((w1[3]/sum(w1[3:4]))/((log(2)/t1.seq)+log(2)/data1[1,3]))+((w1[4]/sum(w1[3:4]))/((log(2)/t1.seq)+log(2)/data1[1,4]))) else event.rate.list.treatment<-NULL
        if(!is.null(target.event.rate.control))    event.rate.list.control<-
            (log(2)/t1.seq)*(((w1[1]/sum(w1[1:2]))/((log(2)/t1.seq)+log(2)/data1[1,1]))+((w1[2]/sum(w1[1:2]))/((log(2)/t1.seq)+log(2)/data1[1,2]))) else event.rate.list.control<-NULL
  }


  if(!is.null(event.rate.list.overall)) MST.overall.cutoff<-t1.seq[abs(event.rate.list.overall-target.event.rate.overall)==min(abs(event.rate.list.overall-target.event.rate.overall))] else MST.overall.cutoff<-NULL
  if(!is.null(event.rate.list.treatment))   MST.treatment.cutoff<-t1.seq[abs(event.rate.list.treatment-target.event.rate.treatment)==min(abs(event.rate.list.treatment-target.event.rate.treatment))] else MST.treatment.cutoff<-NULL
  if(!is.null(event.rate.list.control))     MST.control.cutoff<-t1.seq[abs(event.rate.list.control-target.event.rate.control)==min(abs(event.rate.list.control-target.event.rate.control))] else MST.control.cutoff<-NULL


  list(MST.overall.cutoff=MST.overall.cutoff,
       MST.treatment.cutoff=MST.treatment.cutoff,
       MST.control.cutoff=MST.control.cutoff)
}


censor.rate.fun.retro<-function(data1=mst.data,p.treat=0.5,p.pos.control=.5,p.pos.treat=.5,shape.par=1,my.dist.tmp='exp',sim.nn=10000000,my.follow.time=NULL,my.total.time=NULL,mst.overall=10,mst.treatment=15,mst.control=8,target.event.rate.overall=NULL,target.event.rate.treatment=NULL,target.event.rate.control=NULL,censoringAndMST='Overall')
{
  #data1[1,1]:(1-p.treat)*(1-p.pos.control): control group with negative biomarker
  #data1[1,2]:(1-p.treat)*(p.pos.control): control group with positive biomarker
  #data1[1,3]:(p.treat)*(1-p.pos.treat): treatment group with negative biomarker
  #data1[1,4]:(p.treat)*(p.pos.treat): treatment group with positive biomarker
  w1<-c((1-p.treat)*(1-p.pos.control),(1-p.treat)*(p.pos.control),(p.treat)*(1-p.pos.treat),(p.treat)*(p.pos.treat))

  MST.scale<-MST.cutoff.fun.scale.factor(data1=data1,p.treat=p.treat,p.pos.control=p.pos.control,p.pos.treat=p.pos.treat,target.event.rate=0.5,shape.par=shape.par,my.dist.tmp=my.dist.tmp,mst.overall=mst.overall,mst.treatment=mst.treatment,mst.control=mst.control)

  data1.scale<-data1
  if(censoringAndMST=='Overall') data1.scale<-data1*MST.scale$overall[1] else
  {
    data1.scale[,1:2]<-data1.scale[,1:2]*MST.scale$control[1]
    data1.scale[,3:4]<-data1.scale[,3:4]*MST.scale$treatment[1]
  }

  MST.list.org<-MST.cutoff.fun(data1=data1,p.treat=p.treat,p.pos.control=p.pos.control,p.pos.treat=p.pos.treat,target.event.rate=0.5,shape.par=shape.par,my.dist.tmp=my.dist.tmp)

  MST.list.scale<-MST.cutoff.fun(data1=data1.scale,p.treat=p.treat,p.pos.control=p.pos.control,p.pos.treat=p.pos.treat,target.event.rate=0.5,shape.par=shape.par,my.dist.tmp=my.dist.tmp)

  MST.censoring.list<-censoring.cutoff.fun.retro(data1=data1.scale,p.treat=p.treat,p.pos.control=p.pos.control, target.event.rate.overall=target.event.rate.overall,target.event.rate.treatment=target.event.rate.treatment,target.event.rate.control=target.event.rate.control, shape.par=shape.par,my.dist.tmp=my.dist.tmp)

  if(censoringAndMST=='Overall') time.truncated<-MST.censoring.list$MST.overall.cutoff else time.truncated<-unlist(MST.censoring.list[3:2])
  event.rate.list<-list()
  for(j in 1:1)
  {
    t.tmp<-as.numeric(as.vector(data1.scale[j,]))
    rate.tmp<-rate.tmp1<-numeric()
    if(my.dist.tmp=='exp')
    {
      if(censoringAndMST=='Overall')
      {

       rate.tmp<-(log(2)/time.truncated)/((log(2)/t.tmp)+(log(2)/time.truncated))

      } else
      {

        rate.tmp<-c((log(2)/time.truncated[1])/((log(2)/t.tmp[1:2])+(log(2)/time.truncated[1])),(log(2)/time.truncated[2])/((log(2)/t.tmp[3:4])+(log(2)/time.truncated[2])))

      }
    }

    event.rate.list[[j]]<-list(event.rate.ind=rate.tmp,event.rate.overall=sum(rate.tmp*w1))
  }
  names(event.rate.list)<-c('H1')
  event.rate.list.all<-list(event.rate=event.rate.list,prop.n=w1,censoring.MST.list=time.truncated, scale.factor=MST.scale,data.scale=data1.scale,MST.data.scale=MST.list.scale,MST.data.org=MST.list.org)
  event.rate.list.all
}

#-------
if(study.design=='prospective') censor.rate.est<-censor.rate.fun(data1=mst.data,p.treat=my.p.treat,p.pos.control=my.p.pos.control,p.pos.treat=my.p.pos.treat,target.event.rate=my.target.event.rate,shape.par=shape.par1,my.dist.tmp=my.dist,my.follow.time=my.follow.time,my.total.time=my.total.time) else
  censor.rate.est<-censor.rate.fun.retro(data1=mst.data,p.treat=my.p.treat,p.pos.control=my.p.pos.control,p.pos.treat=my.p.pos.treat,shape.par=shape.par1,my.dist.tmp=my.dist,mst.overall=my.mst.overall,mst.treatment=my.mst.treatment,mst.control=my.mst.control,target.event.rate.overall=my.target.event.rate.overall,target.event.rate.treatment=my.target.event.rate.treatment,target.event.rate.control=my.target.event.rate.control,censoringAndMST=my.censoringAndMST)


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


if(study.design=='prospective') list(mst.data=mst.data,freq=prop.n,HRR=HR1,event.rate=censor.rate.est$event.rate$H1,Peterson.power=Peterson.power,Lachin.power=Lachin,MST.list=censor.rate.est$MST.list) else

  list(mst.data=mst.data,freq=prop.n,HRR=HR1,event.rate=censor.rate.est$event.rate$H1,Peterson.power=Peterson.power,Lachin.power=Lachin,MST.list=censor.rate.est$MST.list,censoring.MST=censor.rate.est$censoring.MST.list, scale.factor=censor.rate.est$scale.factor,data.scale=censor.rate.est$data.scale,MST.data.scale=censor.rate.est$MST.data.scale,MST.data.org=censor.rate.est$MST.data.org)

}



ui <- fluidPage(
  headerPanel('Statistical Power Calculation for Predictive Biomarker'),

  tabsetPanel(

    tabPanel("Prospective Study: Statistical Power Calculation for Predictive Biomarker",
             sidebarLayout(
               sidebarPanel(
                 numericInput("num", label ="Sample size", value = 300),

                 numericInput("MstNegControl", label ="Median survival time (MST) in control group with negative biomarker", value = 3),
                 numericInput("MstPosControl", label ="Median survival time (MST) in control group with positive biomarker", value = 4),
                 numericInput("MstNegTreatment", label ="Median survival time (MST) in treatment group with negative biomarker", value = 2),
                 numericInput("MstPosTreatment", label ="Median survival time (MST) in treatment group with positive biomarker", value = 5),

                 numericInput("alpha",label="Type I error",value = 0.05),


                 sliderInput("pTreatment",label="Percentage of tretament",min = 0, max = 1, value = 0.5),

                 sliderInput("pPosTreatment",label="prevalance of positive biomarker in treatment group",min = 0, max = 1, value =0.5),

                 sliderInput("pPosControl",label="prevalance of positive biomarker in control group",min = 0, max = 1, value =0.5),

                 selectInput("timeunit",label="Time unit",choices=c("years","months","days")),

                 numericInput("timeFollow", label ="Follow up time", value = 1),
                 numericInput("timeTotal", label ="Total time of study", value = 5),

                 actionButton("Submit2","Calculate"),
                 downloadButton('downloadReport')

               ),

               mainPanel(
                 h2("Sample Size Justification"),
                textOutput("summary"),
                tags$hr(),


                 h4("Table 1A: Subgroup Median Survival Time (MST)"),
                 tableOutput("MSTtable"),
                 tags$hr(),

                h4("Table 1B: Overall Median Survival Time (MST)"),
                tableOutput("MSTmarginal"),
                tags$hr(),


                 h4("Table 2: Subgroup Proportion"),
                 tableOutput("Nproportion"),
                 tags$hr(),


                 h4("Table 3: Subgroup Sample Size"),
                 tableOutput("samplesize"),
                 tags$hr(),


                 h4("Table 4: Subgroup Censoring Rate"),
                 tableOutput("eventrate"),
                 tags$hr(),

                 h4("Table 5: Power"),
                 tableOutput("powerTable"),
                 tags$hr(),

                h4("Reference"),
                verbatimTextOutput("reference"),
                tags$hr()

               ))),

tabPanel("Retrospective Study: Statistical Power Calculation for Predictive Biomarker",
         sidebarLayout(
           sidebarPanel(
             numericInput("numRetro", label ="Sample size", value = 174),
             numericInput("alphaRetro",label="Type I error",value = 0.05),


             selectInput("timeunitRetro",label="Time unit",choices=c("months","years","days")),
             numericInput("MstNegControlRetro", label ="Preliminary Data: Median survival time (MST) in control group with negative biomarker", value = 31),
             numericInput("MstPosControlRetro", label ="Preliminary Data: Median survival time (MST) in control group with positive biomarker", value = 9),
             numericInput("MstNegTreatmentRetro", label ="Preliminary Data: Median survival time (MST) in treatment group with negative biomarker", value = 44),
             numericInput("MstPosTreatmentRetro", label ="Preliminary Data: Median survival time (MST) in treatment group with positive biomarker", value = 80),

             selectInput("censoringAndMST",label="Type of Available MST and Censoring Information",choices=c("Overall","Treatment and Control")),
             numericInput("MstOverall", label ="Stud cohort: Overall Median survival time (if choose 'Overall')", value = 80),
             numericInput("censoringOverallRetro", label ="Overall censoring rate (if choose 'Overall')", value = 0.5,min=0,max=1),


             numericInput("MstControl", label ="Stud cohort: Median survival time (MST) in Control group (if choose 'treatment and control')", value = 80),
             numericInput("MstTreatment", label ="Stud cohort: Median survival time (MST) in Treatment group (if choose 'treatment and control')", value = 80),

             numericInput("censoringControlRetro", label ="Censoring rate in Control group (if choose 'treatment and control')", value = 0.5,min=0,max=1),
             numericInput("censoringTreatmentRetro", label ="Censoring rate in Treatment (if choose 'treatment and control')", value = 0.5,min=0,max=1),



             sliderInput("pTreatmentRetro",label="Percentage of tretament",min = 0, max = 1, value = 0.52),

             sliderInput("pPosTreatmentRetro",label="prevalance of positive biomarker in treatment group",min = 0, max = 1, value =0.5),

             sliderInput("pPosControlRetro",label="prevalance of positive biomarker in control group",min = 0, max = 1, value =0.5),

             actionButton("Submit2Retro","Calculate"),
             downloadButton('downloadReportRetro')

           ),

           mainPanel(
             h2("Sample Size Justification"),
             textOutput("summaryRetro"),
             tags$hr(),


             h4("Table 1A: Subgroup Median Survival Time (Raw)"),
             tableOutput("MSTtableRetroOrg"),
             tags$hr(),

             h4("Table 1B: Overall Median Survival Time (Raw)"),
             tableOutput("MSTmarginalRetroOrg"),
             tags$hr(),

             h4("Table 1C: Subgroup Median Survival Time (Scaled)"),
             tableOutput("MSTtableRetroScale"),
             tags$hr(),

             h4("Table 1D: Overall Median Survival Time (Scaled)"),
             tableOutput("MSTmarginalRetroScale"),
             tags$hr(),


             h4("Table 2: Subgroup Proportion"),
             tableOutput("NproportionRetro"),
             tags$hr(),


             h4("Table 3: Subgroup Sample Size"),
             tableOutput("samplesizeRetro"),
             tags$hr(),


             h4("Table 4: Subgroup Censoring Rate"),
             tableOutput("eventrateRetro"),
             tags$hr(),

             h4("Table 5: Power"),
             tableOutput("powerTableRetro"),
             tags$hr(),

             h4("Reference"),
             verbatimTextOutput("referenceRetro"),
             tags$hr()

           )))
)
)




server <- function(input,output){

  #---Prospective study--------

  get.result.prospective <- eventReactive(input$Submit2, {
    nn<-input$num
    pTreatment<-input$pTreatment
    pPosTreatment<-input$pPosTreatment
    pPosControl<-input$pPosControl
    alpha<-input$alpha
    time.unit<-input$timeunit

    MstNegControl<-input$MstNegControl
    MstPosControl<-input$MstPosControl
    MstNegTreatment<-input$MstNegTreatment
    MstPosTreatment<-input$MstPosTreatment

    timeFollow<-input$timeFollow
    timeTotal<-input$timeTotal

    mst.data=data.frame(C_low=c(MstNegControl,MstNegControl),C_high=c(MstPosControl,MstPosControl),T_low=c(MstNegTreatment,MstNegTreatment),T_high=c(MstPosTreatment,MstPosTreatment))

    tmp99<-predictive.power.fun(n.list=nn,mst.data=mst.data,
                                my.p.treat=pTreatment,
                                my.p.pos.control=pPosTreatment,
                                my.p.pos.treat=pPosControl,
                                my.follow.time=timeFollow,
                                my.total.time==timeTotal,
                                alpha=alpha)


    mst.data.tmp<-as.matrix(tmp99$mst.data)[1,]

    mst.table<-cbind(mst.data.tmp[1:2],mst.data.tmp[3:4])
    HRR<-round(prod(mst.data.tmp[c(2,3)])/prod(mst.data.tmp[c(1,4)]),2)
    MST.marginal.tmp<-tmp99$MST.list
    MST.marginal<-as.matrix(unlist(MST.marginal.tmp))
    dimnames(MST.marginal)<-list(sub('.cutoff','',sub('MST.','',names(MST.marginal.tmp))),'MST')
    event.rate.table<-round(cbind(tmp99$event.rate$event.rate.ind[1:2],tmp99$event.rate$event.rate.ind[3:4]),2)
    prop.subgroup<-cbind(tmp99$freq[1:2],tmp99$freq[3:4])
   dimnames(mst.table)<-dimnames(event.rate.table)<-dimnames(prop.subgroup)<-list(c('Negative','Positive'),c('Control','Treatment'))
   mst.table.more<-cbind(mst.table,c(round(mst.data.tmp[1]/mst.data.tmp[3],2),round(mst.data.tmp[2]/mst.data.tmp[4],2)),c(HRR,NA))
   dimnames(mst.table.more)[[2]][3:4]<-c('HR','HRR')
    n.subgroup<-round(prop.subgroup*nn)
    power.table<-round(rbind(tmp99$Peterson.power,tmp99$Lachin.power),2)
     dimnames(power.table)<-list(c('Peterson','Lachin'),'Power')
       a1<-paste('We plan to use a total sample size of ',nn,' to validate the predictive effect of the biomarker. Justification of the sample size is based on the statistical interaction model with a biomarker variable (positive and negative), a treatment variable (treatment and control), and the interaction term of the two variables in the Cox proportional hazard model: $h(t)=h0(t) \\times exp(\\beta_1 \\times biomarker+\\beta_2 \\times treatment+\\beta_3 \\times biomarker \\times treatment)$ where $h(t)$ is a hazard at time $t$ and $h0(t)$ is the baseline hazard. The goal is to test $H_0: \\beta_3=0$ where $\\beta_3$ can be expressed as the log scale of hazard ratio\'s ratio (HRR; described later). The methodology we used is from Peterson and Lachin methods (ref 1 and 2).  Below is the justification of the sample size.',sep='')


     a2<-paste('We assume the median survival time (MST) is ', mst.data.tmp[4],' and ', mst.data.tmp[2],' ', time.unit,' in the treatment and control groups, respectively, for patients with high (or positive) biomarker. The corresponding hazard ratio (HR) is ',round(mst.data.tmp[2]/mst.data.tmp[4],2),'. For patients with the low (or negative) biomarker, we assume their MSTs are ',mst.data.tmp[3],' and ',mst.data.tmp[1],' ',time.unit,' in the treatment and control groups, respectively. The corresponding HR is ',round(mst.data.tmp[1]/mst.data.tmp[3],2),'. Therefore, the hazard ratio\'s ratio (HRR) of high versus low biomarker is ',HRR,'. Table 1A summarizes the MST for each subgroup. In addition, the overall MST (combination of the 4 subgroups) is also obtained as ',round(MST.marginal.tmp$MST.overall.cutoff,2) ,' ',time.unit,' accordingly (Table 1B). Similarly, MSTs for the treatment, control, positive, and negative groups are ', paste(round(as.vector(MST.marginal)[-1],2),collapse=', '),' ',time.unit,', respectively. ',sep='')
     a3<-paste('We also assume that (1) the percentage of patients in the treatment group is ',round(pTreatment*100),'%, and (2) the prevalence of high biomarker is ',round(pPosTreatment*100),'% and ',round(pPosControl*100),'% in the treatment and control groups, respectively. With both assumptions, the subgroup proportion ranges ',min(round(prop.subgroup*100)),'% to ',max(round(prop.subgroup*100)),'% (Table 2). The sample size of subgroup is between ',min(round(n.subgroup)),' and ',max(round(n.subgroup)),' (Table 3).',sep='')
     a4<-paste('Total time of the study will be ',timeTotal, ' ',time.unit,' with ',timeFollow,' ', time.unit,' of follow-up. Under the assumption of uniform distribution for the censoring time, the censoring time follows a uniform distribution between ',timeFollow,' and ',timeTotal,' ', time.unit,'. By comparing with MST in Table 1A (assuming exponential distribution for survival time), we are able to calculate subgroup censoring rate. Table 4 lists censoring rate for each subgroup, ranging from ',min(1-event.rate.table),' to ',max(1-event.rate.table),'.',sep='')
     a5<-paste('By taking all together for consideration with a two-sided ',round(alpha*100),'% type I error, the sample size of ',nn,' will have ',round(power.table[2]*100) ,'% power to detect a HRR of ',HRR,'.',sep='')
     a6<-paste(a1,a2,a3,a4,a5,sep ='\n\n')
     a7<-paste('1. Peterson, B. and S.L. George, Sample size requirements and length of study for testing interaction in a 2 x k factorial design when time-to-failure is the outcome [corrected]. Control Clin Trials, 1993. 14(6): p. 511-22. \n\n 2. Lachin, J.M., Sample size and power for a logrank test and Cox proportional hazards model with multiple groups and strata, or a quantitative covariate with multiple strata. Stat Med, 2013. 32(25): p. 4413-25.')

     censoring.rate.table<-1-event.rate.table

    tmp98<-list(mst.table=mst.table,mst.table.more=mst.table.more,event.rate.table=event.rate.table,censoring.rate.table=censoring.rate.table,prop.subgroup=prop.subgroup,n.subgroup=n.subgroup,power.table=power.table,a1=a1,a2=a2,a3=a3,a4=a4,a5=a5,a6=a6,a7=a7,MST.marginal=MST.marginal)

    return(tmp98)
  })


  output$summary <- renderText({
    get.result.prospective()$a6
  })


  output$MSTtable <- renderTable({
    res <- get.result.prospective()$mst.table.more
  })

  output$MSTmarginal <- renderTable({
    res <- get.result.prospective()$MST.marginal
  })


  output$Nproportion <- renderTable({
    get.result.prospective()$prop.subgroup
  })

  output$samplesize <- renderTable({
    get.result.prospective()$n.subgroup
  })


  output$eventrate <- renderTable({
    get.result.prospective()$censoring.rate.table
  })

  output$powerTable <- renderTable({
    get.result.prospective()$power.table
  })

  output$reference <- renderText({
    get.result.prospective()$a7
  })



  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('predictive_power_prospective.docx',sep='')
    },

    content = function(file) {
      out <- render(paste(system.file("shiny-examples", "myapp", package = "PowerPredictiveBiomarker"),'/predictive_power_prospective.Rmd',sep=''))
      file.rename(out, file)
    }
  )


  #---Retrospective study--------

  get.result <- eventReactive(input$Submit2Retro, {
    nn<-input$numRetro
    pTreatment<-input$pTreatmentRetro
    pPosTreatment<-input$pPosTreatmentRetro
    pPosControl<-input$pPosControlRetro
    alpha<-input$alphaRetro
    time.unit<-input$timeunitRetro

    MstNegControl<-input$MstNegControlRetro
    MstPosControl<-input$MstPosControlRetro
    MstNegTreatment<-input$MstNegTreatmentRetro
    MstPosTreatment<-input$MstPosTreatmentRetro

    censoringAndMST<-input$censoringAndMST
    MstOverall<-input$MstOverall
    censoringOverallRetro<-input$censoringOverallRetro
    MstTreatment<-input$MstTreatment
    MstControl<-input$MstControl
    censoringControlRetro<-input$censoringControlRetro
    censoringTreatmentRetro<-input$censoringTreatmentRetro

    mst.data=data.frame(C_low=c(MstNegControl,MstNegControl),C_high=c(MstPosControl,MstPosControl),T_low=c(MstNegTreatment,MstNegTreatment),T_high=c(MstPosTreatment,MstPosTreatment))

    tmp99<-predictive.power.fun(n.list=nn,mst.data=mst.data,
                                my.p.treat=pTreatment,
                                my.p.pos.control=pPosTreatment,
                                my.p.pos.treat=pPosControl,
                                my.target.event.rate=censoringOverall,
                                my.mst.overall=MstOverall,
                                my.mst.treatment=MstTreatment,
                                my.mst.control=MstControl,
                                my.target.event.rate.overall=1-censoringOverallRetro,
                                my.target.event.rate.treatment=1-censoringTreatmentRetro,
                                my.target.event.rate.control=1-censoringControlRetro,
                                alpha=alpha,
                                study.design='retrospective',
                                my.censoringAndMST=censoringAndMST)

    mst.data.tmp.org<-as.matrix(tmp99$mst.data)[1,]
    HRR<-round(prod(mst.data.tmp.org[c(2,3)])/prod(mst.data.tmp.org[c(1,4)]),2)

    mst.table.org<-cbind(mst.data.tmp.org[1:2],mst.data.tmp.org[3:4])
    dimnames(mst.table.org)[[2]]<-c('Control','Treatment')
    mst.table.more.org<-cbind(mst.table.org,c(round(mst.data.tmp.org[1]/mst.data.tmp.org[3],2),round(mst.data.tmp.org[2]/mst.data.tmp.org[4],2)),c(HRR,NA))
    dimnames(mst.table.more.org)[[2]][3:4]<-c('HR','HRR')

    MST.marginal.tmp<-tmp99$MST.data.org
    MST.marginal.org<-as.matrix(unlist(MST.marginal.tmp))
    dimnames(MST.marginal.org)<-list(sub('.cutoff','',sub('MST.','',dimnames(MST.marginal.org)[[1]])),'MST')
    MST.marginal.org<-MST.marginal.org[1:3,,drop=F]


    mst.data.tmp.scale<-as.matrix(tmp99$data.scale)[1,]
    mst.table.scale<-cbind(mst.data.tmp.scale[1:2],mst.data.tmp.scale[3:4])
    dimnames(mst.table.scale)[[2]]<-c('Control','Treatment')

    mst.table.more.scale<-cbind(mst.table.scale,c(round(mst.data.tmp.scale[1]/mst.data.tmp.scale[3],2),round(mst.data.tmp.scale[2]/mst.data.tmp.scale[4],2)),c(HRR,NA))
    dimnames(mst.table.more.scale)[[2]][3:4]<-c('HR','HRR')
    MST.marginal.scale<-as.matrix(unlist(tmp99$MST.data.scale))
    dimnames(MST.marginal.scale)<-list(sub('.cutoff','',sub('MST.','',dimnames(MST.marginal.scale)[[1]])),'MST')
    MST.marginal.scale<-MST.marginal.scale[1:3,,drop=F]

    event.rate.table<-round(cbind(tmp99$event.rate$event.rate.ind[1:2],tmp99$event.rate$event.rate.ind[3:4]),2)

    prop.subgroup<-cbind(tmp99$freq[1:2],tmp99$freq[3:4])
    dimnames(mst.table.org)<-dimnames(mst.table.scale)<-dimnames(event.rate.table)<-dimnames(prop.subgroup)<-list(c('Negative','Positive'),c('Control','Treatment'))
    n.subgroup<-round(prop.subgroup*nn)
    power.table<-round(rbind(tmp99$Peterson.power,tmp99$Lachin.power),2)
    dimnames(power.table)<-list(c('Peterson','Lachin'),'Power')

    censoring.MST<-tmp99$censoring.MST
    scale.factor<-round(unlist(sapply(tmp99$scale.factor,function(x)x[1])),4)
    data.scale<-tmp99$data.scale
    MST.data.scale<-tmp99$MST.data.scale
    MST.data.org<-tmp99$MST.data.org

    a1<-paste('We plan to use a total sample size of ',nn,' to validate the predictive effect of the biomarker. Justification of the sample size is based on the statistical interaction model with a biomarker variable (positive and negative), a treatment variable (treatment and control), and the interaction term of the two variables in the Cox proportional hazard model: $h(t)=h0(t) \\times exp(\\beta_1 \\times biomarker+\\beta_2 \\times treatment+\\beta_3 \\times biomarker \\times treatment)$ where $h(t)$ is a hazard at time $t$ and $h0(t)$ is the baseline hazard. The goal is to test $H_0: \\beta_3=0$ where $\\beta_3$ can be expressed as the log scale of hazard ratio\'s ratio (HRR; described below). The methodology we used is from Peterson and Lachin methods (ref 1 and 2).  Below is the justification of the sample size.',sep='')
    a20<-paste('Our preliminary data indicate the biomarker is able to detect an HRR of ', HRR,' (effect size) with a hazard ratio (HR) of ',round(mst.data.tmp.org[2]/mst.data.tmp.org[4],2), ' in the positive biomarker (',mst.data.tmp.org[4],' and ',mst.data.tmp.org[2],' ', time.unit,' in the treatment and control groups, respectively) and a HR of ',round(mst.data.tmp.org[1]/mst.data.tmp.org[3],2),' in the negative biomarker (',mst.data.tmp.org[3],' and ',mst.data.tmp.org[1],' ', time.unit,' in the treatment and control groups, respectively) in Table 1A.',sep='')

    a201<-paste('The overall MST (combination of the 4 subgroups) is ',round(MST.marginal.org['overall',1],2) ,' ',time.unit,' accordingly (Table 1B). Similarly, MSTs for the treatment and control groups are ',paste(round(MST.marginal.org[c('treatment','control'),1],2),collapse=' and ') ,' ',time.unit,' respectively.',sep='')
    a20<-paste(a20,a201)


    if(censoringAndMST=='Overall')
      {
        a21<-paste('To apply the preliminary data to the study cohort to testing the effect size (HRR= ', HRR,') for the predictive biomarker, we compare both overall MSTs and find that the preliminary data had a ',if(MST.marginal.org['overall',1]>MstOverall) 'higher' else 'lower',' MST than the study cohort ','(',MST.marginal.org['overall',1],' vs ',MstOverall , ' ', time.unit,'). To make bot data comparable, we rescale the preliminary data with a scale factor of ',scale.factor[1],'. As a result,  the overall MST in Table 1 (after rescale; Table 1C) is same as the overall MST in the study cohort while retaining the same HRR.',sep='')
      } else
      {
        a21<-paste('To apply the preliminary data to the study cohort to testing the effect size (HRR= ', HRR,') for the predictive biomarker, we compare MST of the treatment (and control) between the preliminary data and the study cohort. Results show that for the treatment group, the preliminary data had a ',if(MST.marginal.org['treatment',1]>MstTreatment) 'higher' else 'lower',' MST than the study cohort ','(',MST.marginal.org['treatment',1],' vs ',MstTreatment , ' ', time.unit,'). For the control group, the preliminary data had a ',if(MST.marginal.org['control',1]>MstControl) 'higher' else 'lower',' MST than the study cohort ','(',MST.marginal.org['control',1],' vs ',MstControl , ' ', time.unit,'). To make both data comparable, we rescale the preliminary data with a scale factor of ',paste(scale.factor[2:3],collapse=' and '),' for the treatment and control, respectively. As a result,  MSTs in treatment and control from Table 1 (after rescale; Table 1C) match well to the ones in the study cohort while retaining the same HRR.',sep='')
      }

    a2<-paste('We assume the median survival time (MST) is ', mst.data.tmp.org[4],' and ', mst.data.tmp.org[2],' ', time.unit,' in the treatment and control groups, respectively, for patients with high (or positive) biomarker. The corresponding hazard ratio (HR) is ',round(mst.data.tmp.org[2]/mst.data.tmp.org[4],2),'. For patients with the low (or negative) biomarker, we assume their MSTs are ',mst.data.tmp.org[3],' and ',mst.data.tmp.org[1],' ',time.unit,' in the treatment and control groups, respectively. The corresponding HR is ',round(mst.data.tmp.org[1]/mst.data.tmp.org[3],2),'. Therefore, the hazard ratio\'s ratio (HRR) of high versus low biomarker is ',HRR,'. Table 1A summarizes the MST for each subgroup. In addition, the overall MST (combination of the 4 subgroups) is also obtained as ',round(MST.marginal.tmp$MST.overall.cutoff,2) ,' ',time.unit,' accordingly (Table 1B). Similarly, MSTs for the treatment, control, positive, and negative groups are ', paste(round(as.vector(MST.marginal.org)[-1],2),collapse=', '),' ',time.unit,', respectively. ',sep='')
    a3<-paste('We also assume that (1) the percentage of patients in the treatment group is ',round(pTreatment*100),'%, and (2) the prevalence of positive biomarker is ',round(pPosTreatment*100),'% and ',round(pPosControl*100),'% in the treatment and control groups, respectively. With both assumptions, the subgroup proportion ranges ',min(round(prop.subgroup*100)),'% to ',max(round(prop.subgroup*100)),'% (Table 2). The sample size of subgroup is between ',min(round(n.subgroup)),' and ',max(round(n.subgroup)),' (Table 3).',sep='')

    if(censoringAndMST=='Overall') a4<-paste('Since the study cohort has an overall censoring rate of ',censoringOverallRetro,', the exponential distribution with MST of ',censoring.MST,' was used for the censoring time to ensure the same censoring rates for Table 1C. Table 4 lists censoring rate for each subgroup, ranging from ',min(1-event.rate.table),' to ',max(1-event.rate.table),sep='') else
      a4<-paste('Since the study cohort has a censoring rate of ',censoringTreatmentRetro,' and ', censoringControlRetro,' for the treatment and control, respectively, in order to match the same censoring rates for Table 1C, the exponential distribution with MST of ',paste(round(censoring.MST[2:1],2),collapse=' and '),' was used for the censoring time in the treatment and control groups, respectively. Table 4 lists censoring rate for each subgroup, ranging from ',min(1-event.rate.table),' to ',max(1-event.rate.table),sep='')

    a5<-paste('By taking all together for consideration with a two-sided ',round(alpha*100),'% type I error, the sample size of ',nn,' will have ',round(power.table[2]*100) ,'% power to detect a HRR of ',HRR,'.',sep='')
    #a6<-paste(a1,a20,a21,a2,a3,a4,a5,sep ='\n\n')
    a6<-paste(a1,a20,a21,a3,a4,a5,sep ='\n\n')
    a7<-paste('1. Peterson, B. and S.L. George, Sample size requirements and length of study for testing interaction in a 2 x k factorial design when time-to-failure is the outcome [corrected]. Control Clin Trials, 1993. 14(6): p. 511-22. \n\n 2. Lachin, J.M., Sample size and power for a logrank test and Cox proportional hazards model with multiple groups and strata, or a quantitative covariate with multiple strata. Stat Med, 2013. 32(25): p. 4413-25.')

    censoring.rate.table<-1-event.rate.table

    tmp98<-list(mst.table.org=mst.table.org,mst.table.more.org=mst.table.more.org,
                event.rate.table=event.rate.table,censoring.rate.table=censoring.rate.table,prop.subgroup=prop.subgroup,n.subgroup=n.subgroup,power.table=power.table,a1=a1,a2=a2,a3=a3,a4=a4,a5=a5,a6=a6,a7=a7,MST.marginal.org=MST.marginal.org,MST.marginal.scale=MST.marginal.scale,
    censoring.MST=censoring.MST,
    scale.factor=scale.factor,
    mst.table.more.scale=mst.table.more.scale,
    data.scale=data.scale,
    MST.data.scale=MST.data.scale,
    MST.data.org=MST.data.org)

    mst.data.tmp.scale<-as.matrix(tmp99$data.scale)[1,]
    mst.table.scale<-cbind(mst.data.tmp.scale[1:2],mst.data.tmp.scale[3:4])
    MST.marginal.scale<-as.matrix(unlist(tmp99$MST.data.scale))
    dimnames(MST.marginal.scale)<-list(sub('.cutoff','',sub('MST.','',names(MST.marginal.scale))),'MST')


    return(tmp98)
  })


  output$summaryRetro <- renderText({
    get.result()$a6
  })


  output$MSTtableRetroOrg <- renderTable({
    res <- get.result()$mst.table.more.org
  })

  output$MSTmarginalRetroOrg <- renderTable({
    res <- get.result()$MST.marginal.org
  })

  output$MSTtableRetroScale <- renderTable({
    res <- get.result()$mst.table.more.scale
  })

  output$MSTmarginalRetroScale <- renderTable({
    res <- get.result()$MST.marginal.scale
  })



  output$NproportionRetro <- renderTable({
    get.result()$prop.subgroup
  })

  output$samplesizeRetro <- renderTable({
    get.result()$n.subgroup
  })


  output$eventrateRetro <- renderTable({
    get.result()$censoring.rate.table
  })

  output$powerTableRetro <- renderTable({
    get.result()$power.table
  })

  output$referenceRetro <- renderText({
    get.result()$a7
  })



  output$downloadReportRetro <- downloadHandler(
    filename = function() {
      paste('predictivePowerRetro.docx')
    },

    content = function(file) {
      out <- render(paste(system.file("shiny-examples", "myapp", package = "PowerPredictiveBiomarker"),'/predictive_power_retro.Rmd',sep=''))
      file.rename(out, file)
    }
  )


}

shinyApp(ui=ui,server=server)
