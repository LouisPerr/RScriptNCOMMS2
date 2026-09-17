
#rm(list = ls())

library(deSolve)
library(readxl)
library(bbmle)
library(ggplot2)
library(tidyverse)
library(MASS)
library(tictoc)
library(numbers)
library(plyr); library(dplyr);library(lubridate)
library(metR)
library(plotly)
library(RColorBrewer)
library(geomtextpath)
library(tidyquant)

#change system language
#Sys.setlocale("LC_ALL", "English")

#vaccine waning = 5.5 × 10−3 week-1
#vaccine efficacy rate = 0.94

#let's start at week 24 of 2012, because before then, reporting was bad because of a lack of ressources at the lab! 

#Data, functions and parameters


#VaxData2012<-read_excel("C:/Users/perrlo/Documents/Stuff to Import/Data/Zinsstag2017yellowincasmodelepourrage.xlsx",sheet="2012Vax")  
#VaxData2013<-read_excel("C:/Users/perrlo/Documents/Stuff to Import/Data/Zinsstag2017yellowincasmodelepourrage.xlsx",sheet="2013Vax")  

VaxData2012<-read_excel("/Users/ljperriens/Documents/Documents SWISS TPH/Stuff to Import/Data/Zinsstag2017yellowincasmodelepourrage.xlsx",sheet="2012Vax")  
VaxData2013<-read_excel("/Users/ljperriens/Documents/Documents SWISS TPH/Stuff to Import/Data/Zinsstag2017yellowincasmodelepourrage.xlsx",sheet="2013Vax")  

load("Data2025NDJOnly")

EndemicEQ<-function( Eps_in,Beta_in,N_in )
{
  # parameters
  # defaults : epsilon= 0.1362; beta=3.45e-5; N=30000;
  epsilon=Eps_in;beta=Beta_in;N=N_in;
  mu=6.66e-3;d=6.66e-3;sigma=0.239;delta=1.23;
  vax_efficacy=0.94;
  lambda=5.5e-3; #vaccine waning from Zinsstag 2017
  # quick alphas
  U=1*N; #ie. 0% background vax
  V=N-U;
  alpha=mu*(V/U);alpha=alpha*vax_efficacy
  
  #new
  a=-beta*(sigma+d)*delta
  b=mu*N*beta*sigma - (alpha+d)*(sigma+d)*delta - beta*epsilon*sigma + lambda*( alpha/(d+lambda) )*(sigma+d)*delta
  c= -lambda*( alpha/(d+lambda) )*epsilon*sigma + (alpha+d)*epsilon*sigma
  
  #population to return
  I= max(((-b + sqrt(b^2-(4*a*c))) / (2*a))          , ((-b - sqrt(b^2-(4*a*c))) / (2*a))   )  ;
  S=  (  (sigma+d)*(delta) -( epsilon*sigma)/I  ) /  ( beta*sigma  )
  E= (delta/sigma)*I
  V=(alpha/(d+lambda))*S
  
  Pop0=c(S0=S,E0=E,I0=I,V0=V,N0=S+E+I+V);Pop0
  
  return(Pop0)
}

# INITIAL VALUES AND PARAMETERS
parameters<-c(epsilon=0.1362,beta=3.45e-5,mu=6.66e-3,d=6.66e-3,sigma=0.239,delta=1.23,N=35000)
Phi=100 #for vaccination rate estimation
#initial alpha= 2.93e-3
initialPop<-EndemicEQ(Eps_in = as.numeric(parameters["epsilon"]),Beta_in =as.numeric(parameters["beta"]) , N_in = as.numeric(parameters["N"]))[1:4]; initialPop<-c(initialPop,0);names(initialPop) <- c("S","E","I","V","C");

# #DISEASE-INDUCED DEATH RATE
# delta_0=1.23
# delta_1=7
# Phi_delta=0.3

################################
# VACCINATION RATE SHENANIGANS
################################

VaxDiff<-function( x,y,U,V ){((x-y*exp(-100))*U-V)^2} #function to minimise for estimation of alpha 0.

######### 2012 ######### 

# first week->  October 08, 2012
# last week-> December 31, 2012
Marked_i_12<-VaxData2012$Vdogs2012
V12<-cumsum(Marked_i_12)
U12<-sum(initialPop)-V12
VU12<-matrix(0,nrow = 3,ncol = length(U12))
rownames(VU12) <- c("U","V","Marked_i")
colnames(VU12)<-c(as.character(1:13))
for (i in 1:length(U12))
{ 
  VU12["U",i]<-U12[i]
  VU12["V",i]<-V12[i]
  VU12["Marked_i",i]<-Marked_i_12[i]
}

alphas2012<-matrix(0,nrow = 3,ncol = length(U12))
rownames(alphas2012) <- c("alpha_0","alpha_1","alpha_i")
colnames(alphas2012)<-c(as.character(1:13))

#FIRST WEEK 
alphas2012["alpha_1",1]<- VU12["Marked_i",1]/(VU12["U",1]*(exp(-Phi)-1))
alphas2012["alpha_0",1]<- -alphas2012["alpha_1",1]
alphas2012["alpha_i",1]<- alphas2012["alpha_0",1] + alphas2012["alpha_1",1] *exp(-Phi)

#RECURSION
for (i in 1:12)
{
  cons<-alphas2012["alpha_0",i]+alphas2012["alpha_1",i]*exp(-Phi) #in thus case just alpha_t
  alphas2012["alpha_1",i+1]<-(VU12["Marked_i",i+1]/VU12["U",i+1]-cons)/(exp(-Phi)-1)
  alphas2012["alpha_0",i+1]<-optimize(f=VaxDiff,c(-0.5,0.5),tol=0.0001,y=alphas2012["alpha_1",i+1],U=VU12["U",i+1],VU12["Marked_i",i+1])$minimum
  alphas2012["alpha_i",i+1]<-alphas2012["alpha_0",i+1]+alphas2012["alpha_1",i+1]*exp(-Phi)
}

#rm(VU12,U12,V12,Marked_i_12)
vaxDF2012<-data.frame(matrix(nrow=15,ncol=0))
vaxDF2012$weeks<-c(0:14)
vaxDF2012$vaxdogz<-as.numeric(c(0,alphas2012[3,],0))

ggplot(vaxDF2012,aes(weeks,vaxdogz))+geom_step()+scale_x_continuous(breaks=c(0:14),labels=0:14)+xlab("Vaccination campaign week, 2012")+ylab("Vaccination rate (α)")+ 
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)   , legend.title = element_text(size=14), legend.text = element_text(size=12))   

rm(VU12,U12,V12,Marked_i_12)
######### 2013 ######### 
#first week-> week 91 -> September 30, 2013
#last week-> week 104 -> December 23, 2013
Marked_i_13<-VaxData2013$VDogs2013
V13<-cumsum(Marked_i_13)
U13<-sum(initialPop)-V13
VU13<-matrix(0,nrow = 3,ncol = length(U13))
rownames(VU13) <- c("U","V","Marked_i")
colnames(VU13)<-c(as.character(1:13))

for (i in 1:length(U13))
{
  VU13["U",i]<-U13[i]
  VU13["V",i]<-V13[i]
  VU13["Marked_i",i]<-Marked_i_13[i]
}

alphas2013<-matrix(0,nrow = 3,ncol = length(U13))
rownames(alphas2013) <- c("alpha_0","alpha_1","alpha_i")
colnames(alphas2013)<-c(as.character(1:13))

#FIRST WEEK
alphas2013["alpha_1",1]<- VU13["Marked_i",1]/(VU13["U",1]*(exp(-Phi)-1))
alphas2013["alpha_0",1]<- -alphas2013["alpha_1",1]
alphas2013["alpha_i",1]<- alphas2013["alpha_0",1] + alphas2013["alpha_1",1] *exp(-Phi)

#Recursion 

for (i in 1:12)
{
  cons<-alphas2013["alpha_0",i]+alphas2013["alpha_1",i]*exp(-Phi) #in thus case just alpha_t
  alphas2013["alpha_1",i+1]<-(VU13["Marked_i",i+1]/VU13["U",i+1]-cons)/(exp(-Phi)-1)
  alphas2013["alpha_0",i+1]<-optimize(f=VaxDiff,c(-0.5,0.5),tol=0.0001,y=alphas2013["alpha_1",i+1],U=VU13["U",i+1],VU13["Marked_i",i+1])$minimum
  alphas2013["alpha_i",i+1]<-alphas2013["alpha_0",i+1]+alphas2013["alpha_1",i+1]*exp(-Phi)
}

rm(VU13,U13,V13,Marked_i_13)


#ggplot(as.data.frame(c(alphas2013[3,],0)),aes(1:14,c(alphas2013[3,],0) ))+geom_step()


# 2012
# first week->  October 08, 2012 -> week 19 
#  last week-> December 31, 2012 -> week 31

# 2012
#first week-> September 30, 2013 -> 70
# last week->  December 23, 2013 -> 83




#############################################################################
#  ANALYSIS BEGINNING
#############################################################################
#first week for this analysis: June 2012 
load("Data2025NDJOnly")

NDJ_full<-NDJ

NDJ<-NDJ[which(NDJ$Date>"2012-06-01"),]; NDJ$Week<-1:dim(NDJ)[1] #Oct 08 2012:19->31; Sept 30 2012: 74->86, end 2015:136
#falpha for start on 1st June 2012


#initial_alpha= 2.93e-3
initial_alpha= 0
#falpha <- approxfun(x = c(1:(18),(19):(31),(32):(70),(71):(83),(84:136),137), y = as.numeric(c(rep(0,(18)),alphas2012["alpha_i",],rep(0,39),alphas2013["alpha_i",],rep(0,53),initial_alpha)), method = "constant", rule = 2)
falpha <- approxfun(x = c(1:(18),(19):(31),(32):(70),(71):(83),(84:136),137), y = as.numeric(c(rep(0,(18)),alphas2012["alpha_i",],rep(0,39),alphas2013["alpha_i",],rep(initial_alpha,53),initial_alpha)), method = "constant", rule = 2)


rm(cons,VaxData2012,VaxData2013,vaxDF2012,vaxDF2013,i,Phi,alphas2012,alphas2013,initial_alpha,VaxDiff)

#Remove period containing follow up intervention from back end of the data
NDJ<-NDJ[which(NDJ$Date < "2024-06-11"),];


NDJ_pre<-NDJ[which(NDJ$Date<"2022-06-01"),]
NDJ_during<-NDJ[which(NDJ$Date > "2022-06-01" & NDJ$Date < "2023-07-01"),]
NDJ_post<-NDJ[which(NDJ$Date>"2023-07-01"),] #intervnetion period, 


Key_Info<-matrix(nrow=3,ncol=6)
colnames(Key_Info) <- c("Period","Epsilon","Beta","Recorded Cases","Predicted Cases","Average Weekly Incidence")
Key_Info[,1]<-c("Pre","Intervention","Post")
Key_Info[1,"Recorded Cases"]<-sum(NDJ_pre$TestResult)
Key_Info[2,"Recorded Cases"]<-sum(NDJ_during$TestResult)
Key_Info[3,"Recorded Cases"]<-sum(NDJ_post$TestResult)

Key_Info<-as.data.frame(Key_Info)

############################
# Pre-intervention period  #
############################


SEIV_Incidence_Pre<-function(parameters,initialPop,time_values)
{
  SEIV_equations<-function(time,vars,params,falpha)
  {
    S=vars[1]
    E=vars[2]
    I=vars[3]
    V=vars[4]
    C=vars[5]
    
    with(as.list(params,falpha,time),
         {dS=mu*N-beta*I*S   -falpha(time)*S*0.94 -d*S +(5.5e-3)*V
         dE=beta*I*S- (sigma+d)*E+epsilon  
         dC=sigma*E
         dI=sigma*E -  (delta_0)*I
         dV=-(d+5.5e-3)*V+falpha(time)*S*0.94                  
         return(list(c(dS,dE,dI,dV,dC)))
         })
  }
  
  initial_values <-c(S=initialPop[1],E=initialPop[2],I=initialPop[3],V=initialPop[4],C=initialPop[5])
  
  parameter_values<-c(epsilon=as.numeric(parameters["epsilon"]),beta=as.numeric(parameters["beta"]),mu=as.numeric(parameters["mu"]),d=as.numeric(parameters["d"]),sigma=as.numeric(parameters["sigma"]),delta_0=as.numeric(parameters["delta_0"]),N=as.numeric(sum(initial_values[1:4])))
  out<-ode(y=initial_values,times=time_values,func=SEIV_equations,parms=parameter_values,falpha=falpha)
  as.data.frame(out)
}


LogL_Pre<-function(ML_pars)
{
  Data<-NDJ_pre    #Pre_data
  ML_pars<-exp(ML_pars)
  initialPop_LL<-EndemicEQ(Eps_in = ML_pars["epsilon"],Beta_in = ML_pars["beta"],N=35000)[1:4];initialPop_LL<-c(initialPop_LL,0);names(initialPop_LL) <- c("S","E","I","V","C");
  LL_pars<-c(ML_pars["epsilon"],ML_pars["beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(initialPop_LL))
  times<-seq(1:dim(Data)[1])
  Cumu_incidence<-SEIV_Incidence_Pre(LL_pars,initialPop_LL,times)$C.C
  Incidence<-diff(Cumu_incidence)
  observations<-Data$TestResult[-1]  
  logL<-sum(dpois(x=observations,lambda = Incidence,log=TRUE))
  return(- logL )
}


ML_pars<-c(epsilon=log(0.1362),beta=log(3.45e-5))
MLE_pre_est<-optim(ML_pars,LogL_Pre);exp(MLE_pre_est$par)
Key_Info[1,"Epsilon"]<-exp(MLE_pre_est$par[1]);Key_Info[1,"Beta"]<-exp(MLE_pre_est$par[2])

#simulation period 1:
Input_initialPop<-EndemicEQ(Eps_in = as.numeric(exp(MLE_pre_est$par[1])),Beta_in =as.numeric(exp(MLE_pre_est$par[2])),N=35000)[1:4];Input_initialPop<-c(Input_initialPop,0);names(Input_initialPop) <- c("S","E","I","V","C");
Input_pars<-c(epsilon=as.numeric(exp(MLE_pre_est$par[1]))-0.0,beta=as.numeric(exp(MLE_pre_est$par[2])),mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(Input_initialPop))
#add one in the time horizon to get stating pop for next period
time_values<-as.numeric(seq(1,1+dim(NDJ_pre)[1]))
Pred_Incidence_pre<-SEIV_Incidence_Pre(Input_pars,Input_initialPop,time_values)

Pop_end_P1<- as.data.frame(tail(   Pred_Incidence_pre[c("S.S","E.E","I.I","V.V","C.C")] ,2  )) ;names(Pop_end_P1) <- c("S","E","I","V","C");

Key_Info[1,"Predicted Cases"]<-Pop_end_P1[1,"C"]
Key_Info[1,"Average Weekly Incidence"]<-as.numeric(Key_Info[1,"Predicted Cases"])/dim(NDJ_pre)[1]


##############################################
# Intervention, Beta fixed, estimate epsilon #
##############################################
initial_alpha=0
falpha_post <- approxfun(x = c(1:2,3), y = as.numeric(c(rep(initial_alpha,(2)),initial_alpha)), method = "constant", rule = 2)


#remove the falphas from ODE set
SEIV_Incidence_Intervention<-function(parameters,initialPop,time_values)
{
  SEIV_equations<-function(time,vars,params,falpha_post)
  {
    S=vars[1]
    E=vars[2]
    I=vars[3]
    V=vars[4]
    C=vars[5]
    
    with(as.list(params,falpha_post,time),
         {dS=mu*N -beta*I*S -falpha_post(time)*S*0.94  -d*S + (5.5e-3)*V
         dE=beta*I*S- (sigma+d)*E+epsilon  
         dC=sigma*E
         dI=sigma*E -  (delta_0)*I
         dV=-(d+5.5e-3)*V +falpha_post(time)*S*0.94              
         return(list(c(dS,dE,dI,dV,dC)))
         })
  }
  
  parameter_values<-c(epsilon=as.numeric(parameters["epsilon"]),beta=as.numeric(parameters["beta"]),mu=as.numeric(parameters["mu"]),d=as.numeric(parameters["d"]),sigma=as.numeric(parameters["sigma"]),delta_0=as.numeric(parameters["delta_0"]),N=as.numeric(parameters["N"]))
  initial_values <-c(S=initialPop[1],E=initialPop[2],I=initialPop[3],V=initialPop[4],C=initialPop[5])
  out<-ode(y=initial_values,times=time_values,func=SEIV_equations,parms=parameter_values,falpha_post=falpha_post)
  as.data.frame(out)
}

LogL_Intervention_exp<-function(y)
{
  Data<-NDJ_during
  initialPop_LL<-as.numeric(Pop_end_P1[2,]);names(initialPop_LL)<-c("S","E","I","V","C")
  LL_pars<-c(epsilon=as.numeric(exp(y)),beta=as.numeric(Key_Info$Beta)[1],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(initialPop_LL[1:4]))
  times<-seq(1:dim(Data)[1])
  Cumu_incidence<-SEIV_Incidence_Intervention(LL_pars,initialPop_LL,times)$C.C
  Incidence<-diff(Cumu_incidence)
  observations<-Data$TestResult[-1]
  logL<-sum(dpois(x=observations,lambda = Incidence,log=TRUE))
  return( - logL )
}
ymin<-optimise(LogL_Intervention_exp,c(-10,0)) ;  exp(ymin$minimum)

Key_Info[2,"Epsilon"]<-exp(ymin$minimum);Key_Info[2,"Beta"]<-exp(MLE_pre_est$par[2])

#simulation intervention period
Input_initialPop<-as.numeric(Pop_end_P1[2,]);names(Input_initialPop)<-c("S","E","I","V","C")
Input_pars<-c(epsilon=exp(ymin$minimum),beta=as.numeric(exp(MLE_pre_est$par[2])),mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(Input_initialPop[1:4]))
time_values<-as.numeric(seq(1,dim(NDJ_during)[1]+1))

Pred_Incidence_Intervention<-SEIV_Incidence_Intervention(Input_pars,Input_initialPop,time_values)
Pop_end_P2<- as.data.frame(tail(   Pred_Incidence_Intervention[c("S.S","E.E","I.I","V.V","C.C")] ,2  )) ;names(Pop_end_P2) <- c("S","E","I","V","C");

Key_Info[2,"Predicted Cases"]<- Pop_end_P2[1,"C"]-Pop_end_P1[2,"C"]
Key_Info[2,"Average Weekly Incidence"]<-as.numeric(Key_Info[2,"Predicted Cases"])/length(time_values)

##########################################################################
# Post intervention, Beta and epsilon from pre intervention- Validation #
##########################################################################


SEIV_Incidence_Post<-function(parameters,initialPop,time_values)
{
  SEIV_equations<-function(time,vars,params,falpha_post)
  {
    S=vars[1]
    E=vars[2]
    I=vars[3]
    V=vars[4]
    C=vars[5]
    
    with(as.list(params,time,falpha_post),
         {dS=mu*N -beta*I*S  -d*S + (5.5e-3)*V -falpha_post(time)*S*0.94
         dE=beta*I*S- (sigma+d)*E+epsilon  
         dC=sigma*E
         dI=sigma*E -  (delta_0)*I
         dV=-(d+5.5e-3)*V            +falpha_post(time)*S*0.94  
         return(list(c(dS,dE,dI,dV,dC)))
         })
  }
  
  parameter_values<-c(epsilon=as.numeric(parameters["epsilon"]),beta=as.numeric(parameters["beta"]),mu=as.numeric(parameters["mu"]),d=as.numeric(parameters["d"]),sigma=as.numeric(parameters["sigma"]),delta_0=as.numeric(parameters["delta_0"]),N=as.numeric(parameters["N"]))
  initial_values <-c(S=initialPop[1],E=initialPop[2],I=initialPop[3],V=initialPop[4],C=initialPop[5])
  out<-ode(y=initial_values,times=time_values,func=SEIV_equations,parms=parameter_values,falpha_post=falpha_post)
  as.data.frame(out)
}

Key_Info[3,"Epsilon"]<-Key_Info[1,"Epsilon"] ;Key_Info[3,"Beta"]<-Key_Info[2,"Beta"]

Input_initialPop<-as.numeric(Pop_end_P2[2,]);names(Input_initialPop)<-c("S","E","I","V","C")
Input_pars<-c(epsilon=as.numeric(Key_Info[3,"Epsilon"]),beta=as.numeric(Key_Info[3,"Beta"]),mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(Input_initialPop[1:4]))
time_values<-as.numeric(seq(1,dim(NDJ_post)[1]))

Pred_Incidence_post<-SEIV_Incidence_Post(Input_pars,Input_initialPop,time_values)
Pop_end_P3<- as.data.frame(tail(   Pred_Incidence_post[c("S.S","E.E","I.I","V.V","C.C")] ,1 )) ;names(Pop_end_P3) <- c("S","E","I","V","C");

Key_Info[3,"Predicted Cases"]<- Pop_end_P3["C"]-Pop_end_P2[2,"C"]
Key_Info[3,"Average Weekly Incidence"]<-as.numeric(Key_Info[3,"Predicted Cases"])/length(time_values)

print(Key_Info)
All_results<- cbind( rbind(                       Pred_Incidence_pre[ c(1:(nrow(Pred_Incidence_pre)-1)),],  Pred_Incidence_Intervention[c(1:(nrow(Pred_Incidence_Intervention)-1)),] , Pred_Incidence_post       )       ,                                      NDJ[,c("Date","Week","TestResult")])
trajectories_incidence<-as.data.frame(All_results[,c("Week","Date","TestResult","C.C")])
trajectories_incidence$incidence<-c(0,diff(trajectories_incidence$C.C))

###########################################
# Summary, plots and likelihood heat maps #
###########################################

ggplot(trajectories_incidence  ,aes(x=as.Date(Date)))+
  geom_line(aes(y = C.C),show.legend = FALSE,size=1,color="#3288bd")  +
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = cumsum(TestResult),colour=TestResult),show.legend = FALSE,size=1,color="black") +
  geom_vline(xintercept=as.Date(trajectories_incidence[c(522,578),"Date"]), color="purple", linetype="dashed", size=0.8) +
  ylab("Cumulative simulated incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c(trajectories_incidence[c(19,31,71,83),"Date"])), color="blue", linetype="dashed", size=0.8)+
  theme(axis.text=element_text(size=26), axis.title=element_text(size=26)   , legend.title = element_text(size=14), legend.text = element_text(size=12))  +
  annotate("text",x=c(as.Date(NDJ[25,"Date"])),y=110, label="2012 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date(NDJ[c(76),"Date"])),y=100, label="2013 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date(NDJ[c(550),"Date"])),y=100, label="Border vaccination", angle=0,size=9)


# #CUMULATIVE INCIDENCE NEAR INTERVENTION - GOOD
# ggplot(trajectories_incidence[465:nrow(trajectories_incidence),]  ,aes(x=as.Date(Date)))+
#   geom_line(aes(y = C.C   ) , size=1,color="#3288bd")  + 
#   scale_x_date(date_labels="%b %y",date_breaks="3 month")+
#   geom_line(aes(y =cumsum(trajectories_incidence[,"TestResult"])[464]+ cumsum(TestResult)),show.legend=FALSE,size=1,color="black") +
#   geom_vline(xintercept=as.Date(trajectories_incidence[c(522,578),"Date"]), color="purple", linetype="dashed", size=0.8) +
#   ylab("Cumulative rabies incidence")+xlab("Date")+
#   theme(axis.text=element_text(size=18), axis.title=element_text(size=22)   , legend.title = element_text(size=14), legend.text = element_text(size=12))  +
#   annotate("text",x=c(as.Date(NDJ[c(555),"Date"])),y=190, label="Campaign period", angle=0,size=7)
# 


# INCIDENCE NEAR INTERVENTION    - GOOD                                                                                                        
ggplot(trajectories_incidence[465:nrow(trajectories_incidence),]  ,aes(x=as.Date(Date))) +
  scale_x_date(date_labels="%b %y",date_breaks="3 month")+
  geom_line(aes(y = incidence),show.legend = FALSE,size=1,color="#3288bd")   + 
  geom_point(aes(y = TestResult),show.legend = FALSE) +
  xlab("Date")+ylab("Cases vs modelled weekly incidence")+
  geom_vline(xintercept=as.Date(trajectories_incidence[c(522,578),"Date"]), color="purple", linetype="dashed", size=0.8)+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)   , legend.title = element_text(size=14), legend.text = element_text(size=12))  +
  annotate("text",x=c(as.Date(NDJ[c(550),"Date"])),y=1.5, label="Border vaccination", angle=0,size=7)


#Incidence near intervention, on a monthly time step:
Monthly_TS<-c(round(dim(trajectories_incidence)[1]/4))

Monthly_TS<-as.data.frame(Monthly_TS)
Monthly_TS$monthly_true_cases<-0
Monthly_TS$monthly_sim_incidence<-0
head(trajectories_incidence)
Monthly_TS[1,]<-trajectories_incidence[1,c("Date","TestResult","C.C")]
head(Monthly_TS)
for (i in 1:dim(trajectories_incidence)[1])
{
  temp=1
  if(i%%4==0)
  {
    Monthly_TS[i/4,1]<-trajectories_incidence[i,"Date"]
    Monthly_TS[i/4,3]<-trajectories_incidence[i,"C.C"]
    Monthly_TS[i/4,2]<-sum(trajectories_incidence[c(temp:i),"TestResult"])
    temp=i+1
  }
}
Monthly_TS[,1]<-as.Date(Monthly_TS[,1])
Monthly_TS[,2]<-c(0,diff(Monthly_TS$monthly_true_cases))
Monthly_TS[,3]<-c(0,diff(Monthly_TS$monthly_sim_incidence))
colnames(Monthly_TS)<-c("Date","TestResult","incidence")


n_ma=1
ggplot(Monthly_TS[-c(1),],aes(x=as.Date(Date))) +
  geom_line(aes(y = incidence),show.legend = FALSE,size=1,color="#3288bd",linetype="solid")   + 
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_ma(aes(y = TestResult),n=n_ma,show.legend = FALSE,size=0.7,color="black",linetype="dashed") +
  xlab("Date")+ylab("Simulated weekly incidence vs Quarterly data")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  theme(axis.text=element_text(size=26), axis.title=element_text(size=26)   , legend.title = element_text(size=14), legend.text = element_text(size=12))  +
  annotate("text",x=c(as.Date("2012-11-22")),y=4, label="2012 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date("2013-11-10")),y=3.8, label="2013 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date(Monthly_TS[c(138),"Date"])),y=4, label="Border vaccination", angle=0,size=9)


####### Incidence with moving average, regular weekly timestep




Quick_cols<-c(
  "Simulated" = "#1966AD",
  "Data" = "black")



n_ma_weekly=12
ggplot(trajectories_incidence[-1,],aes(x=as.Date(Date))) +
  geom_line(aes(y = incidence,color="Simulated"),show.legend = TRUE,size=1)   + 
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_ma(aes(y = TestResult,color="Data"),n=n_ma_weekly,show.legend = TRUE,size=0.7,linetype="solid") +
  xlab("Date")+ylab("Weekly simulated incidence")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=1.1, label="2012 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date("2013-11-10")),y=1.00, label="2013 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date(trajectories_incidence[c(550),"Date"])),y=1.1, label="Ring vaccination", angle=00,size=9)+
  scale_color_manual(values = Quick_cols)+
  theme_bw()+
  theme(legend.title = element_blank(), legend.text = element_text(size=26),legend.position='bottom',legend.key.size = unit(1, 'cm')) +
  theme(axis.text=element_text(size=26), axis.title=element_text(size=26) )













#FINAL CUMULATIVE INCIDENCE 2012 onwards - GOOD
ggplot(trajectories_incidence  ,aes(x=as.Date(Date)))+
  geom_line(aes(y = C.C,color="Simulated"),show.legend = TRUE,size=1)  +
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = cumsum(TestResult),colour="Data"),show.legend = TRUE,size=1) +
  geom_vline(xintercept=as.Date(trajectories_incidence[c(522,578),"Date"]), color="purple", linetype="dashed", size=0.8) +
  ylab("Simulated cumulative weekly incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c(trajectories_incidence[c(19,31,71,83),"Date"])), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date(NDJ[25,"Date"])),y=120, label="2012 MDV", angle=0,size=8)+
  annotate("text",x=c(as.Date(NDJ[c(76),"Date"])),y=100, label="2013 MDV", angle=0,size=8)+
  annotate("text",x=c(as.Date(NDJ[c(550),"Date"])),y=100, label="Campaign period", angle=0,size=8)+
  scale_color_manual(values = Quick_cols)+
  theme(legend.title = element_blank(), legend.text = element_text(size=26),legend.position='bottom',legend.key.size = unit(1, 'cm')) +
  theme(axis.text=element_text(size=26), axis.title=element_text(size=26) )






#INCIDENCE VS DATE 2012 onwards - GOOD
ggplot(trajectories_incidence[-1,]  ,aes(x=Date))+
  geom_line(aes(y = incidence),size=1,color="#3288bd")   +
  geom_point(aes(y = TestResult),size=1) +
  geom_vline(xintercept=trajectories_incidence[c(510,583),"Date"], color="purple", linetype="dashed", size=0.8)+
  ylab("Simulated weekly incidence")+
  geom_vline(xintercept=as.Date(c(trajectories_incidence[c(19,31,70,83),"Date"])), color="blue", linetype="dashed", size=0.8)+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22)   , legend.title = element_text(size=14), legend.text = element_text(size=12))  +
  annotate("text",x=c(as.Date(NDJ[25,"Date"])),y=2.5, label="2012 MDV", angle=90,size=7)+
  annotate("text",x=c(as.Date(NDJ[c(76),"Date"])),y=2.5, label="2013 MDV", angle=90,size=7)+
  annotate("text",x=c(as.Date(NDJ[c(547),"Date"])),y=2.5, label="Campaign period", angle=90,size=7)

##############################################################
# Incidence vs Data 2012: Data Moving average Monthly (n=5)  #
##############################################################

###############
# HEATMAP PRE #
############### 

#reset SEIV_Incidence with non-zero alphas, and use same LogL as above


space_length<-317
n_heat_grid <- space_length^2
beta_val_heat_grid<-seq(from=0.70*as.numeric(Key_Info[1,"Beta"]),to=1.17*as.numeric(Key_Info[1,"Beta"]),le=space_length);   max(beta_val_heat_grid);min(beta_val_heat_grid)
eps_val_heat_grid<-seq(from=0.20*as.numeric(Key_Info[1,"Epsilon"]),to=2.5*as.numeric(Key_Info[1,"Epsilon"]),le=space_length);   max(eps_val_heat_grid);min(eps_val_heat_grid)
#colourtesterzz<-cbind(beta_val_heat_grid,eps_val_heat_grid); colnames(colourtesterzz)<-c("beta","epsilon")

HeatGrid <- expand.grid(beta_val_heat_grid, eps_val_heat_grid);names(HeatGrid)<-c("beta","epsilon");
HeatGrid$Likelihood<-0

library(svMisc)

j=1;

#j=i;j

for (i in j:n_heat_grid)
{
  HeatGrid[i,"Likelihood"]<-LogL_Pre(c(epsilon=log(HeatGrid[i,"epsilon"]),beta=log(HeatGrid[i,"beta"])))
  
  progress(i,n_heat_grid)
}


#save(HeatGrid,file="Heatmap_2026_no_background_vax_final")

load("Heatmap_2026_no_background_vax_final")
argmin_HG<-HeatGrid[which.min(HeatGrid[,"Likelihood"]),] ; argmin_HG #this should match previous output


# ggplot(HeatGrid,aes(x=beta,y=epsilon,z=Likelihood,fill=Likelihood,show.legend=FALSE),show.legend=FALSE)+  guides(fill=FALSE)+
#   geom_tile()+
#   geom_contour(aes(z=Likelihood,colour = stat(level)),binwidth=1)+
#   #theme(panel.background = element_rect(fill="lightblue"))+
#   # metR::geom_text_contour(aes(z = Likelihood),size=6) +
#   theme(axis.text=element_text(size=14), axis.title=element_text(size=14)   , legend.title = element_text(size=14), legend.text = element_text(size=12))  +
#   scale_colour_gradient(low="red", high="yellow") +
#   labs(color = "Negative
# Log-Likelihood") +
#   ylab("Importation rate (Ɛ)")+xlab("Transmission rate (β)")+
#   geom_point(aes(x = as.numeric(Key_Info[1,"Beta"]), y = as.numeric(Key_Info[1,"Epsilon"]) ), colour="red",shape=19)+
#   annotate("text",x = as.numeric(Key_Info[1,"Beta"])-0.35e-6, y = as.numeric(Key_Info[1,"Epsilon"])+0.002, label=as.character( 416.1 ), angle=0,size=5)  +
#   scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))



#FINAL HEATMAP

ggplot(HeatGrid[HeatGrid$Likelihood<460,],aes(x=beta,y=epsilon,z=Likelihood,show.legend=FALSE),show.legend=FALSE)+  guides(fill=FALSE)+
  scale_x_continuous(breaks=seq(2e-5,4e-5,0.25e-5))+
  geom_contour(aes(z=Likelihood,colour = stat(level)),binwidth=1)+
  metR::geom_text_contour(skip=0,size=8) +
  scale_colour_gradient(low="red", high="yellow") +
  labs(color = "Negative
Log-Likelihood") +
  ylab("Importation rate (Ɛ)")+xlab("Transmission rate (β)")+
  geom_point(aes(x = as.numeric(Key_Info[1,"Beta"]), y = as.numeric(Key_Info[1,"Epsilon"]) ), color="red",fill="red",shape=21,size=5)+
  annotate("text",x = as.numeric(Key_Info[1,"Beta"])-0.25e-6, y = as.numeric(Key_Info[1,"Epsilon"])+0.005, label=as.character( 416.14 ), angle=0,size=8) + 
  theme_bw()+ 
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28)   , legend.title = element_text(size=28), legend.text = element_text(size=24))




#############################################
# SENSITIVITY ANALYSIS (1) outer ellipses
#############################################

time_plus_1<-as.numeric(seq(1,1+dim(NDJ_pre)[1]))

tolerance=0.02
#HG416<-HeatGrid[which(HeatGrid$Likelihood<416+tolerance & HeatGrid$Likelihood>416-tolerance),]
HG417<-HeatGrid[which(HeatGrid$Likelihood<417+tolerance & HeatGrid$Likelihood>417-tolerance),]
HG418<-HeatGrid[which(HeatGrid$Likelihood<418+tolerance & HeatGrid$Likelihood>418-tolerance),]
HG419<-HeatGrid[which(HeatGrid$Likelihood<419+tolerance & HeatGrid$Likelihood>419-tolerance),]
HG420<-HeatGrid[which(HeatGrid$Likelihood<420+tolerance & HeatGrid$Likelihood>420-tolerance),]
HG421<-HeatGrid[which(HeatGrid$Likelihood<421+tolerance & HeatGrid$Likelihood>421-tolerance),]
HG422<-HeatGrid[which(HeatGrid$Likelihood<422+tolerance & HeatGrid$Likelihood>422-tolerance),]
HG423<-HeatGrid[which(HeatGrid$Likelihood<423+tolerance & HeatGrid$Likelihood>423-tolerance),]
HG424<-HeatGrid[which(HeatGrid$Likelihood<424+tolerance & HeatGrid$Likelihood>424-tolerance),]
HG425<-HeatGrid[which(HeatGrid$Likelihood<425+tolerance & HeatGrid$Likelihood>425-tolerance),]


m2 <- function(dataset, n){
  
  x<-(dataset$beta-mean(dataset$beta))/sd(dataset$beta)
  y<-(dataset$epsilon-mean(dataset$epsilon))/sd(dataset$epsilon)
  subset<-cbind(x,y)
  
  alldist <- as.matrix(dist(subset))
  
  while (nrow(subset) > n) {
    cdists = rowSums(alldist)
    closest <- which(cdists == min(cdists))[1]
    subset <- subset[-closest,]
    alldist <- alldist[-closest,-closest]
  }
  subset[,1]<-(subset[,1]*sd(dataset$beta))+mean(dataset$beta)
  subset[,2]<-(subset[,2]*sd(dataset$epsilon))+mean(dataset$epsilon)
  
  return(subset)
}

#pointz<-as.data.frame(rbind(m2(HG416,2),m2(HG417,2),m2(HG418,2),m2(HG419,2),m2(HG420,2),m2(HG421,2),m2(HG422,2),m2(HG423,2),m2(HG424,2),m2(HG425,2)))
pointz<-as.data.frame(rbind(m2(HG417,2),m2(HG418,2),m2(HG419,2),m2(HG420,2),m2(HG421,2),m2(HG422,2),m2(HG423,2),m2(HG424,2),m2(HG425,2)))


initial_HM_w_pointz<-ggplot(HeatGrid[HeatGrid$Likelihood<460,],aes(x=beta,y=epsilon,z=Likelihood,show.legend=FALSE),show.legend=FALSE)+ 
  scale_x_continuous(breaks=seq(2e-5,4e-5,0.25e-5))+
  geom_contour(aes(z=Likelihood,colour = stat(level)),binwidth=1)+
  metR::geom_text_contour(skip=0,size=8) +
  labs(color = "Negative
Log-Likelihood") +
  theme_bw()+
  theme(axis.text=element_text(size=26), axis.title=element_text(size=26)   , legend.title = element_text(size=26), legend.text = element_text(size=24))  +
  scale_colour_gradient(low="red", high="yellow") +
  ylab("Importation rate (Ɛ)")+xlab("Transmission rate (β)")+
  geom_point(aes(x = as.numeric(Key_Info[1,"Beta"]), y = as.numeric(Key_Info[1,"Epsilon"]) ), colour="gold",shape=16,size=5)+
  annotate("text",x = as.numeric(Key_Info[1,"Beta"])-0.25e-6, y = as.numeric(Key_Info[1,"Epsilon"])+0.005, label=as.character( 416.14 ), angle=0,size=8)  +
  # geom_point(aes(x = pointz[1,1], y = pointz[1,2]),   fill="#75B3D8",shape=21,size=5)+ 
  # geom_point(aes(x = pointz[2,1], y = pointz[2,2]),   fill="#FCA588",shape=21,size=5)+  
  geom_point(aes(x = pointz[1,1], y = pointz[1,2]),   fill="#62A8D2", color="#62A8D2",  shape=21,size=4)+
  geom_point(aes(x = pointz[2,1], y = pointz[2,2]),   fill="#FB8363", color="#FB8363",  shape=21,size=4)+
  geom_point(aes(x = pointz[3,1], y = pointz[3,2]),   fill="#519CCB", color="#519CCB",  shape=21,size=4)+ 
  geom_point(aes(x = pointz[4,1], y = pointz[4,2]),   fill="#FB7252", color="#FB7252",  shape=21,size=4)+
  geom_point(aes(x = pointz[5,1], y = pointz[5,2]),   fill="#4090C5", color="#4090C5",  shape=21,size=4)+
  geom_point(aes(x = pointz[6,1], y = pointz[6,2]),   fill="#F86043", color="#F86043",  shape=21,size=4)+
  geom_point(aes(x = pointz[7,1], y = pointz[7,2]),   fill="#3282BD", color="#3282BD",  shape=21,size=4)+
  geom_point(aes(x = pointz[8,1], y = pointz[8,2]),   fill="#F34C37", color="#F34C37",  shape=21,size=4)+
  geom_point(aes(x = pointz[9,1], y = pointz[9,2]),   fill="#2474B6", color="#2474B6",  shape=21,size=4)+
  geom_point(aes(x = pointz[10,1], y = pointz[10,2]), fill="#ED392B", color="#ED392B",  shape=21,size=4)+
  geom_point(aes(x = pointz[11,1], y = pointz[11,2]), fill="#1966AD", color="#1966AD",  shape=21,size=4)+
  geom_point(aes(x = pointz[12,1], y = pointz[12,2]), fill="#DD2A24", color="#DD2A24",  shape=21,size=4)+
  geom_point(aes(x = pointz[13,1], y = pointz[13,2]), fill="#0E59A2", color="#0E59A2",  shape=21,size=4)+
  geom_point(aes(x = pointz[14,1], y = pointz[14,2]), fill="#CE1B1E", color="#CE1B1E",  shape=21,size=4)+
  geom_point(aes(x = pointz[15,1], y = pointz[15,2]), fill="#084B94", color="#084B94",  shape=21,size=4)+
  geom_point(aes(x = pointz[16,1], y = pointz[16,2]), fill="#AF1117", color="#AF1117",  shape=21,size=4)+
  geom_point(aes(x = pointz[17,1], y = pointz[17,2]), fill="#083D7F", color="#083D7F",  shape=21,size=4)+
  geom_point(aes(x = pointz[18,1], y = pointz[18,2]), fill="#67000D", color="#67000D",  shape=21,size=4)

initial_HM_w_pointz




#########################
# Plots from M2 Heatmap #
#########################
chosen_sample_size=dim(pointz)[1]

time_values<-NDJ$Week
Various_Sim_Curves<- matrix(NA,nrow =length(time_values),ncol=1)
Various_Sim_Curves[,1]<-time_values
Various_Sim_Curves<-as.data.frame(Various_Sim_Curves)
names(Various_Sim_Curves)[1]<-"time"


pointz$LogL<-0
for (i in 1:dim(pointz)[1]){
  ML_pars<-c(epsilon=log(pointz[i,2]),beta=log(pointz[i,1]))
  pointz$LogL[i]<-LogL_Pre(ML_pars)
}
colnames(pointz)<-c("beta","epsilon","LogL")


pointz$eps2<-0
pointz$eps2LogL<-0


for (i in 1:dim(pointz)[1])
{
  
  #period 1
  initialPop_HM<-EndemicEQ(Eps_in = pointz[i,"epsilon"],Beta_in = pointz[i,"beta"],N=35000)[1:4];initialPop_HM<-c(initialPop_HM,0);names(initialPop_HM) <- c("S","E","I","V","C");
  Simulation_HM<-SEIV_Incidence_Pre(parameters<-c(epsilon=pointz[i,"epsilon"] ,beta=pointz[i,"beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,delta_1=7,Phi_delta=0.3,N=   sum(initialPop_HM[1:4])     ),initialPop_HM,  time_plus_1   )
  Pop_end_P1_HM<-as.numeric(tail(Simulation_HM,n=1)[c("S.S","E.E","I.I","V.V","C.C")]);names(Pop_end_P1_HM) <- c("S","E","I","V","C");
  
  #period 2
  LogL_Intervention_exp<-function(y)
  {
    Data<-NDJ_during
    initialPop_LL<-Pop_end_P1_HM
    LL_pars<-c(epsilon=as.numeric(exp(y)),beta=as.numeric(pointz[i,"beta"]),mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(initialPop_LL[1:4]))
    times<-seq(1:dim(Data)[1])
    Cumu_incidence<-SEIV_Incidence_Intervention(LL_pars,initialPop_LL,times)$C.C
    Incidence<-diff(Cumu_incidence)
    observations<-Data$TestResult[-1]
    logL<-sum(dpois(x=observations,lambda = Incidence,log=TRUE))
    return( - logL )
  }
  ymin<-optimise(LogL_Intervention_exp,c(-10,0))
  
  pointz[i,"eps2"]<-exp(ymin$minimum)
  pointz[i,"eps2LogL"]<-ymin$objective
  
  Input_initialPop_HM<-Pop_end_P1_HM
  Input_pars_HM<-c(epsilon=pointz[i,"eps2"],beta=pointz[i,"beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N= sum(Pop_end_P1_HM[1:4])    )
  time_values_HM_P2<-as.numeric(seq(1,1+dim(NDJ_during)[1]))
  Pred_Incidence_Intervention_HM<-SEIV_Incidence_Intervention(Input_pars_HM,Input_initialPop_HM,time_values_HM_P2)
  Pop_end_P2_HM<-as.numeric(tail(Pred_Incidence_Intervention_HM,n=1)[c("S.S","E.E","I.I","V.V","C.C")]);names(Pop_end_P2_HM) <- c("S","E","I","V","C");
  
  
  #period 3
  Input_pars_HM_post<-c(epsilon=pointz[i,"epsilon"],beta=pointz[i,"beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N= sum(Pop_end_P2_HM[1:4])    )
  time_values_HM_P3<-as.numeric(seq(1,dim(NDJ_post)[1]))
  Pred_Incidence_Post_HM<-SEIV_Incidence_Post(Input_pars_HM_post,Pop_end_P2_HM,time_values_HM_P3)
  
  run_sims<-cbind(c(Simulation_HM[1:length(time_plus_1)-1,"C.C"],Pred_Incidence_Intervention_HM[1:length(time_values_HM_P2)-1,"C.C"],Pred_Incidence_Post_HM[,"C.C"]))
  Various_Sim_Curves<-cbind(Various_Sim_Curves,run_sims)
  names(Various_Sim_Curves)[i+1]<-paste("Eps1=",round(pointz[i,"epsilon"],4),"Beta1=",round(pointz[i,"beta"],7),"Eps2=",round(pointz[i,"eps2"],4)      )
  
  
}

Various_Sim_Curves<-cbind(Various_Sim_Curves,NDJ$Date,NDJ$TestResult)
names(Various_Sim_Curves)[chosen_sample_size+2]<-"Date"
names(Various_Sim_Curves)[chosen_sample_size+3]<-"Recorded cases"
Various_Sim_Curves$"Model best estimation"<-trajectories_incidence$C.C


Brewed_colours<-c(
  "B1" = "#62A8D2",
  "B2" = "#519CCB",
  "B3" = "#4090C5",
  "B4" = "#3282BD",
  "B5" = "#2474B6",
  "B6" = "#1966AD",
  "B7" = "#0E59A2",
  "B8" = "#084B94",
  "B9" = "#083D7F",
  "E1" =  "#FB8363",
  "E2" = "#FB7252",
  "E3" = "#F86043",
  "E4" = "#F34C37",
  "E5" =  "#ED392B",
  "E6" =  "#DD2A24",
  "E7" =  "#CE1B1E",
  "E8" =  "#AF1117",
  "E9" = "#67000D",
  "MLE"="gold",
  "Data" = "black")

Sensitivity_analysis_plot<-ggplot(Various_Sim_Curves,aes(x=Date))+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = Various_Sim_Curves[,2],color=   "B1"),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,3],color=  "E1" ),size=1) +  ####this 
  geom_line(aes(y = Various_Sim_Curves[,4],color= "B2"   ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,5],color=   "E2" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,6],color= "B3"  ),size=1)  +    
  geom_line(aes(y = Various_Sim_Curves[,7],color=  "E3"),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,8],color= "B4"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,9],color= "E4"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,10],color="B5" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,11],color="E5" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,12],color= "B6"),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,13],color="E6" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,14],color="B7" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,15],color="E7" ),size=1)+   ## this
  geom_line(aes(y = Various_Sim_Curves[,16],color="B8"),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,17],color="E8"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,18],color="B9"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,19],color="E9"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,22],color="MLE"),size=1.5) +
  geom_line(aes(y = cumsum(Various_Sim_Curves[,21])),colour="black",size=1.2) +
  scale_color_manual(values = Brewed_colours)+xlab("Date")+ylab("Cumulative simulated incidence")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=100, label="2012 MDV", angle=0,size=8)+
  annotate("text",x=c(as.Date("2013-11-10")),y=100, label="2013 MDV", angle=0,size=8)+
  annotate("text",x=c(as.Date(NDJ[c(550),"Date"])),y=100, label="Border vaccination", angle=0,size=8)+
  theme(axis.text=element_text(size=26), axis.title=element_text(size=26)   , legend.title = element_text(size=26), legend.text = element_text(size=24))  


Sensitivity_analysis_plot 



############################################
# SENSITIVITY ANALYSIS (2) INNER ellipses
#############################################

HGOpti<-HeatGrid[which(HeatGrid$Likelihood<argmin_HG$Likelihood+1+tolerance),]


ggplot(HGOpti,aes(x=beta,y=epsilon,z=Likelihood,show.legend=FALSE),show.legend=FALSE)+  
  geom_point()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)   , legend.title = element_text(size=14), legend.text = element_text(size=12))  +
  scale_colour_gradient(low="red", high="yellow") +
  labs(color = "Negative
Log-Likelihood") +
  ylab("Importation rate (Ɛ)")+xlab("Transmission rate (β)")+
  geom_point(aes(x = as.numeric(Key_Info[1,"Beta"]), y = as.numeric(Key_Info[1,"Epsilon"]) ), colour="gold",shape=17,size=7)+
  annotate("text",x = as.numeric(Key_Info[1,"Beta"])-0.1e-6, y = as.numeric(Key_Info[1,"Epsilon"])+0.001, label=as.character( round(argmin_HG$Likelihood,2) ), angle=0,size=6)  

fit2<-lm(epsilon~beta,data=HGOpti)
summary(fit2)


Gen_betas<-as.data.frame(seq(from=as.numeric(min(HGOpti$beta)),to=as.numeric(max(HGOpti$beta)),le=10));names(Gen_betas)<-"beta"
Gen_betas$epsilon<-0
Gen_betas$NegLogL<-0
for (i in 1:dim(Gen_betas)[1])
{
  Gen_betas[i,"epsilon"]<-fit2$coefficients[1]+Gen_betas[i,"beta"]*fit2$coefficients[2]
  Gen_betas[i,"index"]<-i
  temp_pair<-c(epsilon=log(Gen_betas[i,"epsilon"]),beta=log(Gen_betas[i,"beta"]))
  Gen_betas[i,"NegLogL"]<-LogL_Pre(temp_pair)
}


# Gen_betas[,"index"]<-1:10
Gen_betas[,"eps2"]<-0
Gen_betas[,"eps2LogL"]<-0

###########################
#Near output colour coded

ggplot(HGOpti, aes(x = beta, y = epsilon)) +
  scale_x_continuous(breaks=seq(2.9e-5,3.4e-5,0.05e-5))+
  geom_point(size=1) +
  stat_smooth(method = "lm", se=FALSE,colour="black",size=1) +
  geom_point(aes(x = Gen_betas[1,1], y = Gen_betas[1,2]), fill="#9e0142",color="#9e0142",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[2,1], y = Gen_betas[2,2]), fill="#d53e4f",color="#d53e4f",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[3,1], y = Gen_betas[3,2]), fill="#f46d43",color="#f46d43",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[4,1], y = Gen_betas[4,2]), fill="#fdae61",color="#fdae61",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[5,1], y = Gen_betas[5,2]), fill="#fee08b",color="#fee08b",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[6,1], y = Gen_betas[6,2]), fill="#e6f598",color="#e6f598",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[7,1], y = Gen_betas[7,2]), fill="#abdda4",color="#abdda4",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[8,1], y = Gen_betas[8,2]), fill="#66c2a5",color="#66c2a5",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[9,1], y = Gen_betas[9,2]), fill="#3288bd",color="#3288bd",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[10,1], y = Gen_betas[10,2]), fill="#5e4fa2", color="#5e4fa2",shape=21,size=6)+
  geom_point(aes(x = argmin_HG$beta, y = argmin_HG$epsilon), color="red",shape=17, size=6)+
  # annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[1,2]), label=as.character( 1 ), angle=0,size=5, colour="black")  +
  # annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[2,2]), label=as.character( 2), angle=0,size=5, colour="black")  +
  # annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[3,2]), label=as.character( 3), angle=0,size=5, colour="black")  +
  # annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[4,2]), label=as.character( 4), angle=0,size=5, colour="black")  +
  #  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[5,2]), label=as.character( 5), angle=0,size=5, colour="black")  +
  #  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[6,2]), label=as.character( 6), angle=0,size=5, colour="black")  +
  #  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[7,2]), label=as.character( 7), angle=0,size=5, colour="black")  +
  #  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[8,2]), label=as.character( 8), angle=0,size=5, colour="black"  )+
  #  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[9,2]), label=as.character( 9), angle=0,size=5, colour="black")  +
  #  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[10,2]), label=as.character( 10 ), angle=0,size=5, colour="black")  +
  annotate("text",x = as.numeric( argmin_HG$beta)+0.2e-6, y = argmin_HG$epsilon, label=as.character( round(argmin_HG$Likelihood,2 )), angle=0,size=7,colour="red")  +
  ylab("Importation rate (Ɛ)")+xlab("Transmission rate (β)")+ 
  theme_bw()+ 
  theme(axis.text=element_text(size=32), axis.title=element_text(size=32)   , legend.title = element_text(size=32), legend.text = element_text(size=24))

#rwb <- colorRampPalette(colors = c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'))

chosen_sample_size_Opti=dim(Gen_betas)[1]
time_values<-NDJ$Week
Various_Sim_Curves_Opti<- matrix(NA,nrow =length(time_values),ncol=1)
Various_Sim_Curves_Opti[,1]<-time_values
Various_Sim_Curves_Opti<-as.data.frame(Various_Sim_Curves_Opti)
names(Various_Sim_Curves_Opti)[1]<-"time"

for (i in 1:dim(Gen_betas)[1])
{
  
  #period 1
  time_plus_1<-as.numeric(seq(1,1+dim(NDJ_pre)[1]))
  initialPop_HM<-EndemicEQ(Eps_in = Gen_betas[i,"epsilon"],Beta_in = Gen_betas[i,"beta"],N=35000)[1:4];initialPop_HM<-c(initialPop_HM,0);names(initialPop_HM) <- c("S","E","I","V","C");
  Simulation_HM<-SEIV_Incidence_Pre(parameters<-c(epsilon=Gen_betas[i,"epsilon"] ,beta=Gen_betas[i,"beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,delta_1=7,Phi_delta=0.3,N=   sum(initialPop_HM[1:4])     ),initialPop_HM,  time_plus_1   )
  
  Pop_end_P1_HM<-as.numeric(tail(Simulation_HM,n=1)[c("S.S","E.E","I.I","V.V","C.C")]);names(Pop_end_P1_HM) <- c("S","E","I","V","C");
  
  
  #period 2
  LogL_Intervention_exp<-function(y)
  {
    Data<-NDJ_during
    initialPop_LL<-Pop_end_P1_HM
    LL_pars<-c(epsilon=as.numeric(exp(y)),beta=as.numeric(Gen_betas[i,"beta"]),mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(initialPop_LL[1:4]))
    times<-seq(1:dim(Data)[1])
    Cumu_incidence<-SEIV_Incidence_Intervention(LL_pars,initialPop_LL,times)$C.C
    Incidence<-diff(Cumu_incidence)
    observations<-Data$TestResult[-1]
    logL<-sum(dpois(x=observations,lambda = Incidence,log=TRUE))
    return( - logL )
  }
  ymin<-optimise(LogL_Intervention_exp,c(-10,0)) ;  exp(ymin$minimum)
  
  Gen_betas[i,"eps2"]<-exp(ymin$minimum)
  Gen_betas[i,"eps2LogL"]<-ymin$objective
  
  
  Input_initialPop_HM<-Pop_end_P1_HM
  Input_pars_HM<-c(epsilon=Gen_betas[i,"eps2"],beta=Gen_betas[i,"beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N= sum(Pop_end_P1_HM[1:4])    )
  time_values_HM_P2<-as.numeric(seq(1,1+dim(NDJ_during)[1]))
  Pred_Incidence_Intervention_HM<-SEIV_Incidence_Intervention(Input_pars_HM,Input_initialPop_HM,time_values_HM_P2)
  Pop_end_P2_HM<-as.numeric(tail(Pred_Incidence_Intervention_HM,n=1)[c("S.S","E.E","I.I","V.V","C.C")]);names(Pop_end_P2_HM) <- c("S","E","I","V","C");
  
  #period 3
  Input_pars_HM_post<-c(epsilon=Gen_betas[i,"epsilon"],beta=Gen_betas[i,"beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N= sum(Pop_end_P2_HM[1:4])    )
  time_values_HM_P3<-as.numeric(seq(1,dim(NDJ_post)[1]))
  Pred_Incidence_Post_HM<-SEIV_Incidence_Post(Input_pars_HM_post,Pop_end_P2_HM,time_values_HM_P3)
  
  run_sims<-cbind(c(Simulation_HM[1:length(time_plus_1)-1,"C.C"],Pred_Incidence_Intervention_HM[1:length(time_values_HM_P2)-1,"C.C"],Pred_Incidence_Post_HM[,"C.C"]))
  Various_Sim_Curves_Opti<-cbind(Various_Sim_Curves_Opti,run_sims)
  names(Various_Sim_Curves_Opti)[i+1]<-paste("Eps1=",round(Gen_betas[i,"epsilon"],4),"Beta1=",round(Gen_betas[i,"beta"],7),"Eps=2",round(Gen_betas[i,"eps2"],4)      )
  
}

Various_Sim_Curves_Opti<-cbind(Various_Sim_Curves_Opti,NDJ$Date,NDJ$TestResult)
names(Various_Sim_Curves_Opti)[chosen_sample_size_Opti+2]<-"Date"
names(Various_Sim_Curves_Opti)[chosen_sample_size_Opti+3]<-"Recorded cases"
Various_Sim_Curves_Opti$"Model best estimation"<-trajectories_incidence$C.C


Brewed_colours_Opti<-c(
  "01"="#9e0142",
  "02"="#d53e4f",
  "03"="#f46d43",
  "04"="#fdae61",
  "05"="#fee08b",
  "06"="#e6f598",
  "07"="#abdda4",
  "08"="#66c2a5",
  "09"="#3288bd",
  "10"="#5e4fa2",
  "True cases" = "black",
  "MLE"= "black"
)


ggplot(Various_Sim_Curves_Opti,aes(x=Date))+
  geom_line(aes(y = Various_Sim_Curves_Opti[,2],color="01"),size=1)+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = Various_Sim_Curves_Opti[,3],color="02"),size=1) +  
  geom_line(aes(y = Various_Sim_Curves_Opti[,4],color="03"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,5],color="04"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,6],color="05"),size=1)  +    
  geom_line(aes(y = Various_Sim_Curves_Opti[,7],color="06"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,8],color="07"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,9],color="08"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,10],color="09"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,11],color="10"),size=1)+
  geom_line(aes(y = cumsum(Various_Sim_Curves_Opti[,13])),color="black",size=1.0) +
  geom_line(aes(y = Various_Sim_Curves_Opti[,14],color="MLE"),size=1.1) +
  scale_color_manual(values = Brewed_colours_Opti)+xlab("Date")+ylab("Simulated cumulative incidence")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=100, label="2012 MDV", angle=0,size=8)+
  annotate("text",x=c(as.Date("2013-11-10")),y=100, label="2013 MDV", angle=0,size=8)+
  #annotate("text",x=c(as.Date(Various_Sim_Curves_Opti[c(138),"Date"])),y=100, label="Ring vaccination campaign", angle=0,size=8)+
  theme(axis.text=element_text(size=26), axis.title=element_text(size=26)   , legend.title = element_text(size=26), legend.text = element_text(size=24))  



#########################################################
# Importation period I vs II: test and relative diff
#########################################################
#near output
t.test(Gen_betas$epsilon,Gen_betas$eps2,alternative = c("greater"),paired=TRUE)
mean(Gen_betas$epsilon)/mean(Gen_betas$eps2)

#larger frame
t.test(pointz$epsilon,pointz$eps2,alternative = c("greater"),paired=TRUE)
mean(pointz$epsilon)/mean(pointz$eps2)

#######################################
#  BASIC REPRODUCTIVE NUMBER
#######################################
R_0<- ( (as.numeric(Key_Info$Beta[1])*35000)/ parameters["delta"] ) * ( parameters["delta"]/ (parameters["delta"]+parameters["d"])   )



###############
# FINAL PLOTS #
###############

Quick_cols<-c(
  "Simulated" = "#1966AD",
  "Data" = "black")
#IINCIDENCE
n_ma_weekly=12
ggplot(trajectories_incidence[-1,],aes(x=Date)) +
  geom_line(aes(y = incidence,color="Simulated"),show.legend = TRUE,size=1)   + 
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_ma(aes(y = TestResult,color="Data"),n=n_ma_weekly,show.legend = TRUE,size=0.7,linetype="solid") +
  xlab("Date")+ylab("Weekly incidence")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=1.1, label="2012 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2013-11-10")),y=1.00, label="2013 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date(trajectories_incidence[c(550),"Date"])),y=1.1, label="Border vaccination", angle=00,size=10)+
  scale_color_manual(values = Quick_cols)+
  theme_bw()+
  theme(legend.title = element_blank(), legend.text = element_text(size=36),legend.position='bottom',legend.key.size = unit(1, 'cm')) +
  theme(axis.text=element_text(size=36), axis.title=element_text(size=36) )


#CUMULATIVE
ggplot(trajectories_incidence  ,aes(x=Date))+
  geom_line(aes(y = C.C,color="Simulated"),show.legend = TRUE,size=1)  +
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = cumsum(TestResult),colour="Data"),show.legend = TRUE,size=1) +
  geom_vline(xintercept=as.Date(trajectories_incidence[c(522,578),"Date"]), color="purple", linetype="dashed", size=0.8) +
  ylab("Cumulative weekly incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c(trajectories_incidence[c(19,31,71,83),"Date"])), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date(NDJ[25,"Date"])),y=120, label="2012 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date(NDJ[c(76),"Date"])),y=100, label="2013 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date(NDJ[c(550),"Date"])),y=100, label="Border vaccination", angle=0,size=10)+
  scale_color_manual(values = Quick_cols)+
  theme_bw()+
  theme(legend.title = element_blank(), legend.text = element_text(size=36),legend.position='bottom',legend.key.size = unit(1, 'cm')) +
  theme(axis.text=element_text(size=36), axis.title=element_text(size=36) )





#BASE HEATMAP -> RENDER W BIGGER TEXT
ggplot(HeatGrid[HeatGrid$Likelihood<460,],aes(x=beta,y=epsilon,z=Likelihood,show.legend=FALSE),show.legend=FALSE)+  guides(fill=FALSE)+
  scale_x_continuous(breaks=seq(2e-5,4e-5,0.25e-5))+
  geom_contour(aes(z=Likelihood,colour = stat(level)),binwidth=1)+
  metR::geom_text_contour(skip=0,size=8) +
  scale_colour_gradient(low="red", high="yellow") +
  labs(color = "Negative
Log-Likelihood") +
  ylab("Importation rate (Ɛ)")+xlab("Transmission rate (β)")+
  geom_point(aes(x = as.numeric(Key_Info[1,"Beta"]), y = as.numeric(Key_Info[1,"Epsilon"]) ), color="red",fill="red",shape=21,size=5)+
  annotate("text",x = as.numeric(Key_Info[1,"Beta"])-0.25e-6, y = as.numeric(Key_Info[1,"Epsilon"])+0.005, label=as.character( 416.14 ), angle=0,size=8) + 
  theme_bw()+ 
  theme(axis.text=element_text(size=32), axis.title=element_text(size=32)   , legend.title = element_text(size=32), legend.text = element_text(size=24))


#HEATMAP W FAR POINTZ WHITE BACKGROUND
initial_HM_w_pointz<-ggplot(HeatGrid[HeatGrid$Likelihood<460,],aes(x=beta,y=epsilon,z=Likelihood,show.legend=FALSE),show.legend=FALSE)+ 
  scale_x_continuous(breaks=seq(2e-5,4e-5,0.25e-5))+
  geom_contour(aes(z=Likelihood,colour = stat(level)),binwidth=1)+
  metR::geom_text_contour(skip=0,size=8) +
  labs(color = "Negative
Log-Likelihood") +
  theme_bw()+ 
  theme(axis.text=element_text(size=32), axis.title=element_text(size=32)   , legend.title = element_text(size=32), legend.text = element_text(size=24))  +
  scale_colour_gradient(low="red", high="yellow") +
  ylab("Importation rate (Ɛ)")+xlab("Transmission rate (β)")+
  geom_point(aes(x = as.numeric(Key_Info[1,"Beta"]), y = as.numeric(Key_Info[1,"Epsilon"]) ), colour="gold",shape=16,size=5)+
  annotate("text",x = as.numeric(Key_Info[1,"Beta"])-0.25e-6, y = as.numeric(Key_Info[1,"Epsilon"])+0.005, label=as.character( 416.14 ), angle=0,size=8)  +
  # geom_point(aes(x = pointz[1,1], y = pointz[1,2]),   fill="#75B3D8",shape=21,size=5)+ 
  # geom_point(aes(x = pointz[2,1], y = pointz[2,2]),   fill="#FCA588",shape=21,size=5)+  
  geom_point(aes(x = pointz[1,1], y = pointz[1,2]),   fill="#62A8D2", color="#62A8D2",  shape=21,size=4)+
  geom_point(aes(x = pointz[2,1], y = pointz[2,2]),   fill="#FB8363", color="#FB8363",  shape=21,size=4)+
  geom_point(aes(x = pointz[3,1], y = pointz[3,2]),   fill="#519CCB", color="#519CCB",  shape=21,size=4)+ 
  geom_point(aes(x = pointz[4,1], y = pointz[4,2]),   fill="#FB7252", color="#FB7252",  shape=21,size=4)+
  geom_point(aes(x = pointz[5,1], y = pointz[5,2]),   fill="#4090C5", color="#4090C5",  shape=21,size=4)+
  geom_point(aes(x = pointz[6,1], y = pointz[6,2]),   fill="#F86043", color="#F86043",  shape=21,size=4)+
  geom_point(aes(x = pointz[7,1], y = pointz[7,2]),   fill="#3282BD", color="#3282BD",  shape=21,size=4)+
  geom_point(aes(x = pointz[8,1], y = pointz[8,2]),   fill="#F34C37", color="#F34C37",  shape=21,size=4)+
  geom_point(aes(x = pointz[9,1], y = pointz[9,2]),   fill="#2474B6", color="#2474B6",  shape=21,size=4)+
  geom_point(aes(x = pointz[10,1], y = pointz[10,2]), fill="#ED392B", color="#ED392B",  shape=21,size=4)+
  geom_point(aes(x = pointz[11,1], y = pointz[11,2]), fill="#1966AD", color="#1966AD",  shape=21,size=4)+
  geom_point(aes(x = pointz[12,1], y = pointz[12,2]), fill="#DD2A24", color="#DD2A24",  shape=21,size=4)+
  geom_point(aes(x = pointz[13,1], y = pointz[13,2]), fill="#0E59A2", color="#0E59A2",  shape=21,size=4)+
  geom_point(aes(x = pointz[14,1], y = pointz[14,2]), fill="#CE1B1E", color="#CE1B1E",  shape=21,size=4)+
  geom_point(aes(x = pointz[15,1], y = pointz[15,2]), fill="#084B94", color="#084B94",  shape=21,size=4)+
  geom_point(aes(x = pointz[16,1], y = pointz[16,2]), fill="#AF1117", color="#AF1117",  shape=21,size=4)+
  geom_point(aes(x = pointz[17,1], y = pointz[17,2]), fill="#083D7F", color="#083D7F",  shape=21,size=4)+
  geom_point(aes(x = pointz[18,1], y = pointz[18,2]), fill="#67000D", color="#67000D",  shape=21,size=4)



#POINT SELECTION NEAR 

ggplot(HGOpti, aes(x = beta, y = epsilon)) +
  scale_x_continuous(breaks=seq(2.9e-5,3.4e-5,0.05e-5))+
  geom_point(size=1) +
  stat_smooth(method = "lm", se=FALSE,colour="black",size=1) +
  geom_point(aes(x = Gen_betas[1,1], y = Gen_betas[1,2]), fill="#9e0142",color="#9e0142",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[2,1], y = Gen_betas[2,2]), fill="#d53e4f",color="#d53e4f",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[3,1], y = Gen_betas[3,2]), fill="#f46d43",color="#f46d43",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[4,1], y = Gen_betas[4,2]), fill="#fdae61",color="#fdae61",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[5,1], y = Gen_betas[5,2]), fill="#fee08b",color="#fee08b",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[6,1], y = Gen_betas[6,2]), fill="#e6f598",color="#e6f598",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[7,1], y = Gen_betas[7,2]), fill="#abdda4",color="#abdda4",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[8,1], y = Gen_betas[8,2]), fill="#66c2a5",color="#66c2a5",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[9,1], y = Gen_betas[9,2]), fill="#3288bd",color="#3288bd",   shape=21,size=6)+
  geom_point(aes(x = Gen_betas[10,1], y = Gen_betas[10,2]), fill="#5e4fa2", color="#5e4fa2",shape=21,size=6)+
  geom_point(aes(x = argmin_HG$beta, y = argmin_HG$epsilon), color="red",shape=17, size=6)+
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[1,2]), label=as.character( 1 ), angle=0,size=5, colour="#9e0142")  +
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[2,2]), label=as.character( 2), angle=0,size=5, colour="#d53e4f")  +
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[3,2]), label=as.character( 3), angle=0,size=5, colour="#f46d43")  +
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[4,2]), label=as.character( 4), angle=0,size=5, colour="#fdae61")  +
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[5,2]), label=as.character( 5), angle=0,size=5, colour="#fee08b")  +
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[6,2]), label=as.character( 6), angle=0,size=5, colour="gold")  +
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[7,2]), label=as.character( 7), angle=0,size=5, colour="#abdda4")  +
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[8,2]), label=as.character( 8), angle=0,size=5, colour="#66c2a5"  )+
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[9,2]), label=as.character( 9), angle=0,size=5, colour="#3288bd")  +
  annotate("text",x = as.numeric( Gen_betas[1,1])-0.1e-6, y = as.numeric(Gen_betas[10,2]), label=as.character( 10 ), angle=0,size=5, colour="#5e4fa2")  +
  annotate("text",x = as.numeric( argmin_HG$beta)+0.2e-6, y = argmin_HG$epsilon, label=as.character( round(argmin_HG$Likelihood,2 )), angle=0,size=5,colour="red")  +
  ylab("Importation rate (Ɛ)")+xlab("Transmission rate (β)")+ 
  theme_bw()+ 
  theme(axis.text=element_text(size=32), axis.title=element_text(size=32)   , legend.title = element_text(size=32), legend.text = element_text(size=24))

#SENSITIVITY INNER 

ggplot(Various_Sim_Curves_Opti,aes(x=Date))+
  geom_line(aes(y = Various_Sim_Curves_Opti[,2],color="01"),size=1)+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = Various_Sim_Curves_Opti[,3],color="02"),size=1) +  
  geom_line(aes(y = Various_Sim_Curves_Opti[,4],color="03"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,5],color="04"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,6],color="05"),size=1)  +    
  geom_line(aes(y = Various_Sim_Curves_Opti[,7],color="06"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,8],color="07"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,9],color="08"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,10],color="09"),size=1)+
  geom_line(aes(y = Various_Sim_Curves_Opti[,11],color="10"),size=1)+
  geom_line(aes(y = cumsum(Various_Sim_Curves_Opti[,13])),color="black",size=1.0) +
  geom_line(aes(y = Various_Sim_Curves_Opti[,14],color="MLE"),size=1.1) +
  scale_color_manual(values = Brewed_colours_Opti)+xlab("Date")+ylab("Cumulative weekly incidence")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=120, label="2012 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2013-11-10")),y=100, label="2013 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2022-12-13")),y=100, label="Border vaccination", angle=0,size=10)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=30))  

#SENSITIVITY OUTER 

Sensitivity_analysis_plot<-ggplot(Various_Sim_Curves,aes(x=Date))+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = Various_Sim_Curves[,2],color=   "B1"),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,3],color=  "E1" ),size=1) +  ####this 
  geom_line(aes(y = Various_Sim_Curves[,4],color= "B2"   ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,5],color=   "E2" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,6],color= "B3"  ),size=1)  +    
  geom_line(aes(y = Various_Sim_Curves[,7],color=  "E3"),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,8],color= "B4"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,9],color= "E4"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,10],color="B5" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,11],color="E5" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,12],color= "B6"),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,13],color="E6" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,14],color="B7" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,15],color="E7" ),size=1)+   ## this
  geom_line(aes(y = Various_Sim_Curves[,16],color="B8"),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,17],color="E8"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,18],color="B9"  ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,19],color="E9"  ),size=1)+
  # geom_line(aes(y = Various_Sim_Curves[,20],color= "B10" ),size=1)+
  # geom_line(aes(y = Various_Sim_Curves[,21],color= "E10" ),size=1)+
  geom_line(aes(y = Various_Sim_Curves[,22],color="MLE"),size=1.5) +
  geom_line(aes(y = cumsum(Various_Sim_Curves[,21])),colour="black",size=1.2) +
  scale_color_manual(values = Brewed_colours)+xlab("Date")+ylab("Cumulative weekly incidence")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=120, label="2012 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2013-11-10")),y=100, label="2013 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2022-12-13")),y=100, label="Border vaccination", angle=0,size=10)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=30))  

Sensitivity_analysis_plot 





###########################################################
# UNDER-REPORTING SENSITIVITY ANALYSIS: BINOMIAL THINNING
##########################################################


# Let Pi be the reporting probability- To check sensitivity to UR, we do not estimate the under-reporting rate, but fi it, and see how the estimated parameters vary. 
Pi<-c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
UR_Table<-as.data.frame(Pi)
UR_Table$Beta<-0
UR_Table$Eps1<-0
UR_Table$Eps2<-0
UR_Table$Predicted_Cases_P1<-0
UR_Table$Expected_Cases_P1<-0
UR_Table$Predicted_Cases_Total<-0
UR_Table$Expected_Cases_Total<-0

UR_sim_runs<-NDJ

for (i in 1:dim(UR_Table)[1])
{
  
  #Period 1 estimation
  
  LogL_Pre_UR<-function(ML_pars)
  {
    Data<-NDJ_pre    #Pre_data
    ML_pars<-exp(ML_pars)
    initialPop_LL<-EndemicEQ(Eps_in = ML_pars["epsilon"],Beta_in = ML_pars["beta"],N=35000)[1:4];initialPop_LL<-c(initialPop_LL,0);names(initialPop_LL) <- c("S","E","I","V","C");
    
    LL_pars<-c(ML_pars["epsilon"],ML_pars["beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(initialPop_LL))
    times<-seq(1:dim(Data)[1])
    Cumu_incidence<-SEIV_Incidence_Pre(LL_pars,initialPop_LL,times)$C.C
    Incidence<-diff(Cumu_incidence)
    observations<-Data$TestResult[-1]  
    logL<-sum(dpois(x=observations,lambda = Incidence*UR_Table[i,"Pi"],log=TRUE))
    return(- logL )
  }
  
  MLE_pre_UR<-optim(ML_pars,LogL_Pre_UR);exp(MLE_pre_UR$par)
  UR_Table[i,"Eps1"]<-exp(MLE_pre_UR$par[1]);UR_Table[i,"Beta"]<-exp(MLE_pre_UR$par[2])
  
  
  #period 1 simulation
  initialPop_UR<-EndemicEQ(Eps_in = UR_Table[i,"Eps1"],Beta_in = UR_Table[i,"Beta"],N=35000)[1:4];initialPop_UR<-c(initialPop_UR,0);names(initialPop_UR) <- c("S","E","I","V","C");
  Simulation_UR<-SEIV_Incidence_Pre(parameters<-c(epsilon=UR_Table[i,"Eps1"] ,beta=UR_Table[i,"Beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,delta_1=7,Phi_delta=0.3,N=   sum(initialPop_UR[1:4])     ),initialPop_UR,  time_plus_1   )
  Pop_end_P1_UR<-as.numeric(tail(Simulation_UR,n=1)[c("S.S","E.E","I.I","V.V","C.C")]);names(Pop_end_P1_UR) <- c("S","E","I","V","C");
  
  UR_Table[i,"Predicted_Cases_P1"]<-as.numeric(Pop_end_P1_UR[5])
  UR_Table[i,"Expected_Cases_P1"]<-sum(NDJ_pre$TestResult)/UR_Table[i,"Pi"]
  
  
  
  
  #Period 2 Estimation
  LogL_Intervention_UR<-function(y)
  {
    Data<-NDJ_during
    initialPop_LL<-Pop_end_P1_UR
    LL_pars<-c(epsilon=as.numeric(exp(y)),beta=as.numeric(UR_Table[i,"Beta"]),mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N=sum(initialPop_LL[1:4]))
    times<-seq(1:dim(Data)[1])
    Cumu_incidence<-SEIV_Incidence_Intervention(LL_pars,initialPop_LL,times)$C.C
    Incidence<-diff(Cumu_incidence)
    observations<-Data$TestResult[-1]
    logL<-sum(dpois(x=observations,lambda = Incidence*UR_Table[i,"Pi"],log=TRUE))
    return( - logL )
  }
  
  ymin<-optimise(LogL_Intervention_UR,c(-10,0)) ;  exp(ymin$minimum)
  UR_Table[i,"Eps2"]<-exp(ymin$minimum)
  
  
  #period 2 simulation
  Input_initialPop_UR<-Pop_end_P1_UR
  Input_pars_UR<-c(epsilon=UR_Table[i,"Eps2"],beta=UR_Table[i,"Beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N= sum(Pop_end_P1_UR[1:4])    )
  time_values_UR_P2<-as.numeric(seq(1,1+dim(NDJ_during)[1]))
  Pred_Incidence_Intervention_UR<-SEIV_Incidence_Intervention(Input_pars_UR,Input_initialPop_UR,time_values_UR_P2)
  Pop_end_P2_UR<-as.numeric(tail(Pred_Incidence_Intervention_UR,n=1)[c("S.S","E.E","I.I","V.V","C.C")]);names(Pop_end_P2_UR) <- c("S","E","I","V","C");
  
  
  #period 3 Simulation/ Validation
  Input_pars_UR_post<-c(epsilon=UR_Table[i,"Eps1"],beta=UR_Table[i,"Beta"],mu=6.66e-3,d=6.66e-3,sigma=0.239,delta_0=1.23,N= sum(Pop_end_P2_UR[1:4])    )
  time_values_UR_P3<-as.numeric(seq(1,dim(NDJ_post)[1]))
  Pred_Incidence_Post_UR<-SEIV_Incidence_Post(Input_pars_UR_post,Pop_end_P2_UR,time_values_UR_P3)
  UR_Table[i,"Predicted_Cases_Total"]<-tail(Pred_Incidence_Post_UR,n=1)$C.C
  UR_Table[i,"Expected_Cases_Total"]<-sum(NDJ$TestResult)/UR_Table[i,"Pi"]
  
  
  
  # put runs together per reporting prob i
  run_sims<-cbind(c(Simulation_UR[1:length(time_plus_1)-1,"C.C"],Pred_Incidence_Intervention_UR[1:length(time_values_UR_P2)-1,"C.C"],Pred_Incidence_Post_UR[,"C.C"]))
  run_sims<-c(0,diff(run_sims))
  UR_sim_runs<-cbind(UR_sim_runs,run_sims)
  names(UR_sim_runs)[i+3]<-paste("Rep_prob=",UR_Table[i,"Pi"])
  
  
  
}

colorRampPalette(brewer.pal(9,"Blues"))(100)   


UR_Table
UR_sim_runs<-UR_sim_runs[-1,]

t.test(UR_Table$Eps1,UR_Table$Eps2,alternative = c("greater"),paired=TRUE)
mean(UR_Table$Eps1)/mean(UR_Table$Eps2)



Brewed_colours_UR<-c(
  "0.3"="#92C3DE",
  "0.4"= "#71B1D7",
  "0.5"= "#58A1CE" ,
  "0.6"="#4191C5",
  "0.7"= "#1562A9",
  "0.8"= "#08478E",
  "0.9"="#08306B",
  "1"="black"
)

#2800x1400
True_incidence_plot<-ggplot(UR_sim_runs,aes(x=Date))+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = UR_sim_runs[,4],color=   "0.3"),size=1)+
  geom_line(aes(y = UR_sim_runs[,5],color=  "0.4" ),size=1) +  ####this 
  geom_line(aes(y = UR_sim_runs[,6],color= "0.5"   ),size=1)+
  geom_line(aes(y = UR_sim_runs[,7],color=   "0.6" ),size=1)+
  geom_line(aes(y = UR_sim_runs[,8],color= "0.7"  ),size=1)  +    
  geom_line(aes(y = UR_sim_runs[,9],color=  "0.8"),size=1)+
  geom_line(aes(y = UR_sim_runs[,10],color= "0.9"  ),size=1)+
  geom_line(aes(y = UR_sim_runs[,11],color= "1"  ),size=1)+
  scale_color_manual(values = Brewed_colours_UR)+xlab("Date")+ylab("Simulated true incidence")+
  labs(color = "Reporting
Probability")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=1.8, label="2012 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2013-11-10")),y=1.6, label="2013 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2022-12-13")),y=1.6, label=" Border vaccination", angle=0,size=10)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_text(size=30), legend.text = element_text(size=30))  

True_incidence_plot


Reported_incidence_plot<-ggplot(UR_sim_runs,aes(x=Date))+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = UR_sim_runs[,4]*0.3,color=   "0.3"),size=1)+
  geom_line(aes(y = UR_sim_runs[,5]*0.4,color=  "0.4" ),size=1) +  ####this 
  geom_line(aes(y = UR_sim_runs[,6]*0.5,color= "0.5"   ),size=1)+
  geom_line(aes(y = UR_sim_runs[,7]*0.6,color=   "0.6" ),size=1)+
  geom_line(aes(y = UR_sim_runs[,8]*0.7,color= "0.7"  ),size=1)  +    
  geom_line(aes(y = UR_sim_runs[,9]*0.8,color=  "0.8"),size=1)+
  geom_line(aes(y = UR_sim_runs[,10]*0.9,color= "0.9"  ),size=1)+
  geom_line(aes(y = UR_sim_runs[,11]*1,color= "1"  ),size=1)+
  scale_color_manual(values = Brewed_colours_UR)+xlab("Date")+ylab("Simulated reported incidence")+
  labs(color = "Reporting
Probability")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=0.6, label="2012 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2013-11-10")),y=.6, label="2013 DMV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2022-12-13")),y=0.6, label=" Border vaccination", angle=0,size=10)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_text(size=30), legend.text = element_text(size=30))  

Reported_incidence_plot





###################################
# STOCHASTIC MODEL COMPARAISON
###################################

# INITIAL VALUES AND PARAMETERS
beta_SEIV<-as.numeric(Key_Info$Beta[1])
eps_SEIV_pre<-as.numeric(Key_Info$Epsilon[1])
eps_SEIV_int<-as.numeric(Key_Info$Epsilon[2])

parameters<-c(epsilon=eps_SEIV_pre,beta=beta_SEIV,mu=6.66e-3,d=6.66e-3,sigma=0.239,delta=1.23,N=35000)

initialPop<-EndemicEQ(Eps_in = as.numeric(parameters["epsilon"]),Beta_in =as.numeric(parameters["beta"]) , N_in = as.numeric(parameters["N"]))[1:4]; initialPop<-c(initialPop,0);names(initialPop) <- c("S","E","I","V","C");
alphas2012<-c( 0.02460852 ,0.005369375,  0.011279878, 0.000720511, 0.02438747 , 0.09795482, 0.2789810   , 0.06400862, 0.1644845  , 0.09993705 ,0.01864542,  0.023627322 ,0.01299291)
alphas2013<- c(0.02123361, 0.013955696,  0.009938354, 0.013243540, 0.002048698, 0.03720184, 0.030324127 , 0.1597335 , 0.19858560 ,0.18416825 , 0.3279123, 0.07518223 ,0.04201111 )

#background vaccination rate
initial_alpha<-0
#initial_alpha<-0.0003294947  # 5% 

falpha <- approxfun(x = c(1:(18),(19):(31),(32):(70),(71):(83),(84:136),137), y = as.numeric( c(rep(0,(18)),alphas2012,rep(0,39),alphas2013,rep(initial_alpha,53),initial_alpha) ), method = "constant", rule = 2)

#need to decide on value of future epsilon
fepsilon <- approxfun(x = c(1:522,523:578,579), y = as.numeric(  c(   rep(eps_SEIV_pre,522),     rep(eps_SEIV_int,56),  eps_SEIV_pre   ) ), method = "constant", rule = 2)



#simulation settings
tau<-0.005
T= dim(NDJ)[1] 
iterations<-T/tau
sim_replications<-10001
mu=6.66e-3;d=6.66e-3;sigma=0.239;delta=1.23;w=5.5e-3;N=35000;alpha=0.000;
beta=beta_SEIV;epsilon=eps_SEIV_pre
simulated_incidence<-matrix(0,nrow = iterations+1,ncol=sim_replications+1)
colnames(simulated_incidence)<-c("time",1:sim_replications)
simulated_incidence<-as.data.frame(simulated_incidence)
simulated_incidence[,"time"]<-seq(0,T,by=tau)

library(svMisc)



simulated_incidence_weekly<-matrix(0,nrow = T,ncol=sim_replications+1)
colnames(simulated_incidence_weekly)<-c("time",1:sim_replications)
simulated_incidence_weekly<-as.data.frame(simulated_incidence_weekly)
simulated_incidence_weekly[,"time"]<-seq(1,T,by=1)

####
# Pre-intervention period, mostly to assess the quality of the predictions
####

j=1

for (j in 1:sim_replications )
{
  
  # State initialisation per run
  PopMat<-matrix(0,nrow = iterations+1, ncol = 6)
  colnames(PopMat)<-c("S","E","I","V","time","cumulative_cases")
  PopMat[1,c(1:4)]<-round(initialPop[1:4],digits=0)
  PopMat[1,"time"]<-0
  PopMat[1,"cumulative_cases"]<-0
  
  cases_per_week<-as.data.frame(matrix(0,nrow=T,ncol=3));colnames(cases_per_week)<-c("time","cumulative_incidence","incidence")
  
  #round starting conditions up or down 
  start_E<-runif(1,min=0,max=1)
  start_I<-runif(1,min=0,max=1)
  
  if ( start_E >= as.numeric(initialPop["E"])-trunc(initialPop["E"]) )
  { start_E <- 2 } else { start_E<-3}
  
  if ( start_I >= as.numeric(initialPop["I"]))
  { start_I <- 0 } else { start_I<-1}
  PopMat[1,c("E","I")]<-c(start_E,start_I)
  
  for (i in 1:(iterations))
  {
    
    
    N_t<-sum(PopMat[i,c(1,4)])
    alpha<- falpha(tau*i)
    #epsilon<-fepsilon(tau*i)
    
    # Rates for Poisson distributions per tau
    rate_transmission <-   (beta *PopMat[i,"S"]* PopMat[i,"I"] )*tau  #transmission
    rate_importation<-tau*fepsilon(tau*i) #importation
    rate_birth <- N*mu*tau  #birth
    
    # transition probability over time tau
    p_vax <-   1 - exp(-alpha*tau)    #vaccination
    p_death <- 1 - exp(-d*tau)  #natural death
    p_rabies_death <- 1 - exp(-delta*tau) #disease induced death
    p_wane <- 1 - exp(-w*tau)    #vax waning
    
    p_leave_exposed <-   1 - exp(-(d+sigma)*tau)  
    n_leave_exposed <- rbinom(1,as.numeric(PopMat[i,"E"]), p_leave_exposed)
    n_exposed_to_infectious<-  rbinom(1,n_leave_exposed,sigma/(sigma+d))
    
    # number of events over time tau
    n_vax <- rbinom(1,as.numeric(PopMat[i,"S"]), p_vax) 
    n_death_S <- rbinom(1,as.numeric(PopMat[i,"S"]), p_death)
    n_death_V <- rbinom(1,as.numeric(PopMat[i,"V"]), p_death)
    n_rabies_death <- rbinom(1,as.numeric(PopMat[i,"I"]), p_rabies_death)
    n_vax_wane<-rbinom(1,as.numeric(PopMat[i,"V"]), p_wane)
    
    n_birth <- rpois(1,rate_birth)
    n_transmission <-rpois(1,rate_transmission)
    n_import <- rpois(1,rate_importation) 
    
    # update system
    PopMat[i+1,"S"] <- as.numeric(PopMat[i,"S"]) + n_birth - n_vax - n_transmission - n_death_S + n_vax_wane
    PopMat[i+1,"E"] <- as.numeric(PopMat[i,"E"]) + n_transmission - n_leave_exposed + n_import
    PopMat[i+1,"I"] <- as.numeric(PopMat[i,"I"]) + n_exposed_to_infectious - n_rabies_death 
    PopMat[i+1,"V"] <- as.numeric(PopMat[i,"V"]) + n_vax - n_death_V- n_vax_wane
    PopMat[i+1,"time"] <- as.numeric(PopMat[i,"time"]) + tau
    PopMat[i+1,"cumulative_cases"] <- PopMat[i,"cumulative_cases"]+n_exposed_to_infectious
    
  }
  

  cases_per_week[1,c("time","cumulative_incidence")]<-PopMat[1,c("time","cumulative_cases")] 
  
  for (k in 1:(dim(cases_per_week)[1]-1))
  {
    {cases_per_week[k+1,c("time","cumulative_incidence")]<-PopMat[1+k*(1/tau),c("time","cumulative_cases")] }
    
  }

  
  cases_per_week$incidence<-c(0,diff(cases_per_week$cumulative_incidence))
  simulated_incidence_weekly[,j+1]<-cases_per_week$incidence
  
  progress(j,sim_replications)
  
}

#save(simulated_incidence,file="TauLeap0.005_10001Replications_Tau_stepsize_FINAL")
#save(simulated_incidence_weekly,file="TauLeap0.005_10001Replications_FINAL")


load("TauLeap0.005_10001Replications_Tau_stepsize_FINAL")
load("TauLeap0.005_10001Replications_FINAL")











#MEAN simulation values
averaged_out_simulation<- cbind(simulated_incidence_weekly$time,  rowMeans(simulated_incidence_weekly[,-1] ))  
colnames(averaged_out_simulation)<-c("time","incidence")

ggplot(averaged_out_simulation[-1,], aes(x=time))+
  geom_ma(aes(y=incidence),n=1)
mean(colSums(simulated_incidence_weekly[,-1])) #target: 255 we get 264 



##################

#median simulation values: get all indices hitting the true number of cases 255
sims_max_index<-sim_replications+1
#median_sim_indexes <-  as.numeric( which(colSums(simulated_incidence_weekly[,2:sims_max_index] ) ==     round(median(colSums(simulated_incidence_weekly[,2:sims_max_index]))    )       ) )    +1
median_sim_indexes <-  as.numeric( which(colSums(simulated_incidence_weekly[,2:sims_max_index] ) ==     255   ) )    +1

median_sims<-simulated_incidence_weekly[,c(1,median_sim_indexes)]

#we have 33 simulations hitting the median predicted incidence (subject to change in re-runs)
#we choose four at random without replacement (choose a random seed..)
sampled_runs<-sample(2:dim(median_sims)[2],4,replace=FALSE,set.seed(2026))

median_sims$Date<-All_results$Date
some_stoca_runs<-median_sims[,c(sample(2:dim(median_sims)[2],4,replace=FALSE,set.seed(2026)))]
some_stoca_runs$Date<-median_sims$Date
some_stoca_runs$Data<-All_results$TestResult
colnames(some_stoca_runs)<-c("run1","run2","run3","run4","Date","TestResult")

#individual plots, with MA..
ggplot(median_sims, aes(x=time))+
  geom_ma(aes(y=median_sims[,sampled_runs[1]]),n=10)+ylab("Simulated Incidence")

ggplot(median_sims, aes(x=time))+
  geom_ma(aes(y=median_sims[,sampled_runs[2]]),n=10)+ylab("Simulated Incidence")

ggplot(median_sims, aes(x=Date))+
  geom_ma(aes(y=median_sims[,sampled_runs[3]]),n=10)+ylab("Simulated incidence")
  
ggplot(median_sims, aes(x=time))+
  geom_ma(aes(y=median_sims[,sampled_runs[4]]),n=10)+ylab("Simulated Incidence")

#individual plots, on a monthly basis

library(clock)

Monthly_summary_Stoca<-some_stoca_runs %>%
  mutate(date = date_group(Date, "month")) %>%
  group_by(date) %>%
  summarise(Run_1 = sum(run1, na.rm = TRUE),Run_2 = sum(run2, na.rm = TRUE),Run_3 = sum(run3, na.rm = TRUE),Run_4 = sum(run4, na.rm = TRUE),TestResult = sum(TestResult, na.rm = TRUE), .groups = "drop")


rand_stoca_colours<-c(
  "Run 9054"="#E69F00",   
  "Run 7607" ="#009E73",      #keep
  "Run 345"="#F0E442",  #keep
  "Run 868" = "#56B4E9",  #keep
  "Data"="black")

png(filename = "RandomStoca1.png", width = 2600, height = 1200);
#incidence, Monthly
ggplot(Monthly_summary_Stoca, aes(x=date))+
  ylim(0,max(Monthly_summary_Stoca[,c(2,3,4,5)]))+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y=Run_1,color="Run 9054"),linetype="solid",size=0.8)+
  geom_line(aes(y=Run_2,color="Run 7607"),linetype="solid",size=0.8)+
  geom_line(aes(y=TestResult,color="Data"),linetype="dashed",size=1)+
  scale_color_manual(values = rand_stoca_colours)+  ylab("Monthly incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="longdash", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="longdash", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=11, label="2012 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2013-11-10")),y=10, label="2013 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2022-12-04")),y=10, label=" Border vaccination", angle=0,size=11)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=34),legend.key.size = unit(1, 'cm'))  
dev.off()


png(filename = "RandomStoca2.png", width = 2600, height = 1200);
ggplot(Monthly_summary_Stoca, aes(x=date))+
  ylim(0,max(Monthly_summary_Stoca[,c(2,3,4,5)]))+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y=Run_3,color="Run 345"),linetype="solid",size=0.8)+
  geom_line(aes(y=Run_4,color="Run 868"),linetype="solid",size=0.8)+
  geom_line(aes(y=TestResult,color="Data"),linetype="dashed",size=1)+
  scale_color_manual(values = rand_stoca_colours)+  ylab("Monthly incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="longdash", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="longdash", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=11, label="2012 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2013-11-10")),y=10, label="2013 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2022-12-04")),y=10, label=" Border vaccination", angle=0,size=11)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=34),legend.key.size = unit(1, 'cm'))  
dev.off()


png(filename = "CumuStoca.png", width = 2600, height = 1200);
#cumulative incidence, Monthly
ggplot(Monthly_summary_Stoca, aes(x=date))+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y=cumsum(Run_1),color="Run 9054"),linetype="solid",size=1)+
  geom_line(aes(y=cumsum(Run_2),color="Run 7607"),linetype="solid",size=1)+
  geom_line(aes(y=cumsum(Run_3),color="Run 345"),linetype="solid",size=1)+
  geom_line(aes(y=cumsum(Run_4),color="Run 868"),linetype="solid",size=1)+
  geom_line(aes(y=cumsum(TestResult),color="Data"),linetype="dashed",size=1)+
  scale_color_manual(values = rand_stoca_colours)+  ylab("Monthly cumulative incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="longdash", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="longdash", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=110, label="2012 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2013-11-10")),y=100, label="2013 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2022-12-04")),y=100, label=" Border vaccination", angle=0,size=11)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=34),legend.key.size = unit(1, 'cm'))  
dev.off()


#We can do the same with simulations hitting the mean predicted value (258), instead of median (231)
#Median predicted number seems to suffer of down-ward bias? Zero-inflated simulations
mean_sim_indexes <-  as.numeric( which(colSums(simulated_incidence_weekly[,2:sims_max_index] ) ==     round(mean(colSums(simulated_incidence_weekly[,2:sims_max_index]))    )       ) )    +1
mean_sims<-simulated_incidence_weekly[,c(1,mean_sim_indexes)]

#columns and names...
sample(1:dim(mean_sims)[2],4,replace=FALSE,set.seed(10))
colnames(mean_sims[,c(sample(1:dim(mean_sims)[2],4,replace=FALSE,set.seed(10)))])

mean_sims$Date<-All_results$Date
some_stoca_runs_mean<-mean_sims[,c( sample(1:dim(mean_sims)[2],4,replace=FALSE,set.seed(10)) )]
some_stoca_runs_mean$Date<-mean_sims$Date
some_stoca_runs_mean$Data<-All_results$TestResult
colnames(some_stoca_runs_mean)<-c("run1","run2","run3","run4","Date","TestResult")

#individual plots, with MA..
ggplot(mean_sims, aes(x=time))+
  geom_ma(aes(y=mean_sims[,11]),n=10)+ylab("Simulated Incidence")

ggplot(mean_sims, aes(x=time))+
  geom_ma(aes(y=mean_sims[,9]),n=10)+ylab("Simulated Incidence")

ggplot(mean_sims, aes(x=Date))+
  geom_ma(aes(y=mean_sims[,10]),n=10)+ylab("Simulated incidence")

ggplot(mean_sims, aes(x=time))+
  geom_ma(aes(y=mean_sims[,16]),n=10)+ylab("Simulated Incidence")

#individual plots, on a monthly basis

Monthly_summary_Stoca_mean<-some_stoca_runs_mean %>%
  mutate(date = date_group(Date, "month")) %>%
  group_by(date) %>%
  summarise(Run_1 = sum(run1, na.rm = TRUE),Run_2 = sum(run2, na.rm = TRUE),Run_3 = sum(run3, na.rm = TRUE),Run_4 = sum(run4, na.rm = TRUE),TestResult = sum(TestResult, na.rm = TRUE), .groups = "drop")

rand_stoca_colours_mean<-c(
  "Run 2444"="gold",
  "Run 2346"= "red",
  "Run 2363"= "#08306B",
  "Run 3482"="#009E73",
  "Data"="black")


#incidence, Monthly, mean sim (258)

png(filename = "RandomStocaMean1.png", width = 2600, height = 1200);
ggplot(Monthly_summary_Stoca_mean, aes(x=date))+
  ylim(0,17)+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y=Run_1,color="Run 2444"),linetype="solid",size=0.8)+
  geom_line(aes(y=Run_3,color="Run 2363"),linetype="solid",size=0.8)+
  geom_line(aes(y=TestResult,color="Data"),linetype="dashed",size=1,alpha=0.8)+
  scale_color_manual(values = rand_stoca_colours_mean)+  ylab("Monthly incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="longdash", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="longdash", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=14, label="2012 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2013-11-10")),y=13, label="2013 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2022-12-04")),y=13, label=" Border vaccination", angle=0,size=11)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=34),legend.key.size = unit(1, 'cm'))  
dev.off()


png(filename = "RandomStocaMean2.png", width = 2600, height = 1200);
ggplot(Monthly_summary_Stoca_mean, aes(x=date))+
  ylim(0,17)+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y=Run_2,color="Run 2346"),linetype="solid",size=0.8)+
  geom_line(aes(y=Run_4,color="Run 3482"),linetype="solid",size=0.8)+
  geom_line(aes(y=TestResult,color="Data"),linetype="dashed",size=1,alpha=0.8)+
    scale_color_manual(values = rand_stoca_colours_mean)+  ylab("Monthly incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="longdash", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="longdash", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=14, label="2012 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2013-11-10")),y=13, label="2013 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2022-12-04")),y=13, label=" Border vaccination", angle=0,size=11)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=34),legend.key.size = unit(1, 'cm'))  
dev.off()



#cumulative incidence, Monthly, MEAN sim (258)
png(filename = "RandomStocaMeanCUMU.png", width = 2600, height = 1200);

ggplot(Monthly_summary_Stoca_mean, aes(x=date))+
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y=cumsum(Run_1),color="Run 2444"),linetype="solid",size=0.8)+
  geom_line(aes(y=cumsum(Run_2),color="Run 2346"),linetype="solid",size=0.8)+
  geom_line(aes(y=cumsum(Run_3),color="Run 2363"),linetype="solid",size=0.8)+
  geom_line(aes(y=cumsum(Run_4),color="Run 3482"),linetype="solid",size=0.8)+
  geom_line(aes(y=cumsum(TestResult),color="Data"),linetype="dashed",size=1)+
  scale_color_manual(values = rand_stoca_colours_mean)+  ylab("Monthly cumulative incidence")+xlab("Date")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="longdash", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="longdash", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=110, label="2012 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2013-11-10")),y=100, label="2013 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2022-12-04")),y=100, label=" Border vaccination", angle=0,size=11)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=34),legend.key.size = unit(1, 'cm'))  

dev.off()


# 50% sims-> cut off 25% of each side (Interquartile range) 
a<-colSums(simulated_incidence_weekly[,2:sims_max_index])
a<-as.data.frame(a)
a$index<-1:10001
colnames(a)<-c("colsums","index")
ordered_a<-a %>% arrange(colsums)


q1<-round((sim_replications)/4)
q3<-3*q1

a50percent<-ordered_a[q1:q3,"index"]+1  #because of time columns in first column
IQR_sims<-simulated_incidence_weekly[,c(1,a50percent)]

mean(colSums(IQR_sims[,-1]))

averaged_IQR_sims<- cbind(simulated_incidence_weekly$time,  rowMeans(IQR_sims[,-1] ))  
colnames(averaged_IQR_sims)<-c("time","incidence")

ggplot(averaged_IQR_sims[-1,], aes(x=time))+
  geom_ma(aes(y=incidence),n=1)

sum(averaged_IQR_sims[,2]) #good


dataset_deter_stoca<-as.data.frame( cbind(All_results$Week,trajectories_incidence$incidence,averaged_IQR_sims[,2],averaged_out_simulation[,2]   )        )
dataset_deter_stoca$Date<-All_results$Date
colnames(dataset_deter_stoca)<-c("Week","Deterministic","Tauleap_IQR","Tauleap","Date")
dataset_deter_stoca$TestResult<-All_results$TestResult

dataset_deter_stoca$randomMedian1<-median_sims[,11]
dataset_deter_stoca$randomMedian2<-median_sims[,3]
dataset_deter_stoca$randomMedian3<-median_sims[,7]
dataset_deter_stoca$randomMedian4<-median_sims[,44]

dataset_deter_stoca$randomMean1<-mean_sims[,11]
dataset_deter_stoca$randomMean2<-mean_sims[,9]
dataset_deter_stoca$randomMean3<-mean_sims[,10]
dataset_deter_stoca$randomMean4<-mean_sims[,16]

#save(dataset_deter_stoca,file="BunchOfSimulations")
load("BunchOfSimulations") #update stocas and deter. 


#cases per period...IQR
iqr_p1<-dataset_deter_stoca[which(dataset_deter_stoca$Date<"2022-06-01"),"Tauleap_IQR"];mean(iqr_p1);sum(iqr_p1);
iqr_p2<-dataset_deter_stoca[which(dataset_deter_stoca$Date > "2022-06-01" & dataset_deter_stoca$Date < "2023-07-01"),"Tauleap_IQR"];mean(iqr_p2);sum(iqr_p2);
iqr_p3<-dataset_deter_stoca[which(dataset_deter_stoca$Date >"2023-07-01"),"Tauleap_IQR"];mean(iqr_p3);sum(iqr_p3);

tau_p2<-dataset_deter_stoca[which(dataset_deter_stoca$Date > "2022-06-01" & dataset_deter_stoca$Date < "2023-07-01"),"Tauleap"];mean(tau_p2);sum(tau_p2);



stoca_colours<-c(
  "Data (moving average)"= "black",
  "Deterministic"="gold",
  "τ-leap mean"= "#CE1B1E", 
  "τ-leap IQR mean"= "#08306B" )
#71B1D7


png(filename = "StocaDeterministic.png", width = 2600, height = 1200);

ggplot(dataset_deter_stoca[-1,],aes(x=Date))+
    scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y=Deterministic,color="Deterministic"),size=1.1)+
  geom_ma(aes(y=TestResult,color="Data (moving average)"),n=5)+
  geom_line(aes(y=Tauleap,color="τ-leap mean"),size=0.8)+
  geom_line(aes(y=Tauleap_IQR,color="τ-leap IQR mean"),size=0.8)+
    xlab("Date (in years)")+ylab("Weekly incidence")+
  scale_color_manual(values = stoca_colours)+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="longdash", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="longdash", size=0.8)+
  annotate("text",x=c(as.Date("2012-11-22")),y=1, label="2012 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2013-11-10")),y=0.8, label="2013 DMV", angle=0,size=11)+
  annotate("text",x=c(as.Date("2022-12-04")),y=0.8, label=" Border vaccination", angle=0,size=11)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=34),legend.key.size = unit(1, 'cm'))  

dev.off()




########################################################################
# RAW DATA PLOTS
########################################################################

library(clock)

Monthly_summary<-trajectories_incidence %>%
  mutate(date = date_group(Date, "month")) %>%
  group_by(date) %>%
  summarise(Monthly_Confirmed_Cases = sum(TestResult, na.rm = TRUE),Monthly_Sim_incidence = sum(incidence,na.rm=TRUE), .groups = "drop")

Quick_cols_monthly<-c(
  "Simulated" = "#1966AD",
  "Data (raw)" = "black")

ggplot(Monthly_summary[-145,],aes(x=as.Date(date))) +
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  geom_line(aes(y = Monthly_Confirmed_Cases,color="Data (raw)"),show.legend = TRUE,size=1,linetype="dashed") +
  xlab("Date")+ylab("Monthly confirmed cases")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=30))  +
  annotate("text",x=c(as.Date("2012-11-22")),y=4, label="2012 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date("2013-11-10")),y=3.8, label="2013 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date("2022-12-13")),y=4, label="Border vaccination", angle=0,size=9)+
  scale_color_manual(values = Quick_cols_monthly)



#Full data set from 2005
ggplot(NDJ_full,aes(x=as.Date(Date))) +
  geom_point(aes(y = TestResult),show.legend = TRUE,size=1.5)   + 
  scale_x_date(date_labels="%b %y",date_breaks="24 month")  +
  xlab("Date")+ylab("Weekly confirmed cases")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  geom_vline(xintercept=as.Date(c("2010-01-01","2011-12-31")), color="red", linetype="dashed", size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=30))  +
  annotate("text",x=c(as.Date("2012-11-22")),y=1.5, label="2012 MDV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2013-11-10")),y=1.3, label="2013 MDV", angle=0,size=10)+
  annotate("text",x=c(as.Date("2022-12-13")),y=1.5, label="Border vaccination", angle=0,size=10)+
  annotate("text",x=c(as.Date("2010-12-31")),y=1.5, label="Lab renovations", angle=0,size=10)

#raw data set from 2012
ggplot(NDJ,aes(x=as.Date(Date))) +
  geom_point(aes(y = TestResult),show.legend = TRUE,size=1.5)   + 
  scale_x_date(date_labels="%b %y",date_breaks="12 month")  +
  xlab("Date")+ylab("Weekly confirmed cases")+
  geom_vline(xintercept=as.Date(c("2022-06-01","2023-07-01")), color="purple", linetype="dashed", size=0.8) +
  geom_vline(xintercept=as.Date(c("2012-10-08","2012-12-31","2013-09-30","2013-12-23")), color="blue", linetype="dashed", size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=34), axis.title=element_text(size=34)   , legend.title = element_blank(), legend.text = element_text(size=30))  +
  annotate("text",x=c(as.Date("2012-11-22")),y=1.5, label="2012 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date("2013-11-10")),y=1.3, label="2013 MDV", angle=0,size=9)+
  annotate("text",x=c(as.Date("2022-12-13")),y=1.5, label="Border vaccination", angle=0,size=9)
