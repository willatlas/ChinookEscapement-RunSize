##########################
###    Chinook MARSS models

setwd("C:/data/chinook")

library("MARSS")
library("car")
library("corrplot")
require(ggplot2)
require(ggfortify)
require(dplyr)


##########################################################################
#### Chinook escapement trend analysis - the data

CK_esc<-read.csv("CK_esc_FINAL.csv",header=TRUE)
QMat<-read.csv("corr_mat_CK.csv",header=FALSE)

#################################################################
#### Salish Sea - starting analysis at 1979 (when Skagit Summer and Fall Start)

Salish_CK<-subset(CK_esc,CK_esc$group_1 >= 1 & CK_esc$group_1 <= 3)
Salish_CK<-subset(Salish_CK,Salish_CK$year>1978)

ck_dat1<-matrix(data=NA, nrow=length(seq(1979,2019,1)),ncol=length(unique(Salish_CK$pop_no))+1)
dimnames(ck_dat1)[[2]]<-c("year",unique(Salish_CK$population))

years<-seq(1979,2019,1)
ck_dat1[,1]<-years

pop.nos<-unique(Salish_CK$pop_no)

for(i in 1:length(pop.nos)){
  dat<-subset(Salish_CK,Salish_CK$pop_no==pop.nos[i])
  ck_dat1[which(ck_dat1[,1]==dat$year[1]):length(years),i+1]<-log(dat[,8])
}

dat.ck1 = t(ck_dat1)
years = dat.ck1[1,] # the years
n = nrow(dat.ck1)-1 # the number of time series
dat.ck1 = dat.ck1[2:nrow(dat.ck1),]


### MARSS mod1 - one trend for all populations
#Z.model = factor(c(seq(1,n,1))) # each population is independent
Z.model = factor(c(rep(1,n))) # each population shares a trend - model estimates intercept values for each population relative to X0
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "equal" # populations share a long-term trend
Q.model = "diagonal and equal" # equal process errors

mod_Salish1a = MARSS(dat.ck1, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                     control=list(maxit=200000,conv.test.slope.tol = 0.09))

### extracting parameters
paramCK1a<-MARSSparamCIs(mod_Salish1a)

### MARSS model 1b - unique trends for each population, correlated process errors
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
## make sure I'm calling this u.model properly
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5],QMat[1:n,6],
                   QMat[1:n,7],QMat[1:n,8],QMat[1:n,9],QMat[1:n,10],QMat[1:n,11],QMat[1:n,12],
                   QMat[1:n,13],QMat[1:n,14],QMat[1:n,15],QMat[1:n,16],QMat[1:n,17],QMat[1:n,18],
                   QMat[1:n,19],QMat[1:n,20],QMat[1:n,21],QMat[1:n,22],QMat[1:n,23]),n,n)


mod_Salish1b = MARSS(dat.ck1, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                     control=list(maxit=200000,conv.test.slope.tol = 0.09))

### extracting parameters
paramCK1b<-MARSSparamCIs(mod_Salish1b)

###################################################################################################
#### Washington Coast Chinook
WACoast_CK<-subset(CK_esc,group_1 >= 4 & group_1 < 7)
WACoast_CK<-WACoast_CK[-which(WACoast_CK$population=="hoko"),]## Dropping hoko from analysis

WA_pops<-unique(WACoast_CK$population)#subset(CK_pops,CK_pops$group_1 >= 4 & CK_pops$group_1 < 7)

### starting at 1980
WACoast_CK<-subset(WACoast_CK, WACoast_CK$year>1979)

ck_dat2<-matrix(data=NA, nrow=length(seq(1980,2019,1)),ncol=length(unique(WACoast_CK$pop_no))+1)
dimnames(ck_dat2)[[2]]<-c("year",unique(WACoast_CK$population))

years<-seq(1980,2019,1)
ck_dat2[,1]<-years

pop.nos<-unique(WACoast_CK$pop_no)

for(i in 1:length(pop.nos)){
  dat<-subset(WACoast_CK,WACoast_CK$pop_no==pop.nos[i])
  ck_dat2[which(ck_dat2[,1]==dat$year[1]):length(years),i+1]<-log(dat[,8])
}

## transposing the data for MARSS
dat.ck2 = t(ck_dat2)
years = dat.ck2[1,] # the years
n = nrow(dat.ck2)-1 # the number of time series
dat.ck2 = dat.ck2[2:nrow(dat.ck2),]

### MARSS model 2b - unique trends for both populations (all run timing) since 1980 - correlated process error
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5]
                   ,QMat[1:n,6],QMat[1:n,7],QMat[1:n,8],QMat[1:n,9],QMat[1:n,10]
                   ,QMat[1:n,11]),n,n)

mod_WACoast2b = MARSS(dat.ck2, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model),  
                      control=list(maxit=200000,conv.test.slope.tol = 0.09))

CK_params<-MARSSparamCIs(mod_WACoast2b)

##############################################################################################
#### Oregon Coast Chinook
ORCoast_CK<-subset(CK_esc,group_1 == 13)

### starting at 1975
ORCoast_CK<-subset(ORCoast_CK, ORCoast_CK$year>1974)

ck_dat4<-matrix(data=NA, nrow=length(seq(1975,2019,1)),ncol=length(unique(ORCoast_CK$pop_no))+1)
dimnames(ck_dat4)[[2]]<-c("year",unique(ORCoast_CK$population))

years<-seq(1975,2019,1)
ck_dat4[,1]<-years

pop.nos<-unique(ORCoast_CK$pop_no)
unique(ORCoast_CK$population)

for(i in 1:length(pop.nos)){
  dat<-subset(ORCoast_CK,ORCoast_CK$pop_no==pop.nos[i])
  ck_dat4[which(ck_dat4[,1]==dat$year[1]):length(years),i+1]<-log(dat[,8])
}

## transposing the data for MARSS
dat.ck4 = t(ck_dat4)
years = dat.ck4[1,] # the years
n = nrow(dat.ck4)-1 # the number of time series
dat.ck4 = dat.ck4[2:nrow(dat.ck4),]

## Setting north umpqua observer error to zero (Dam counts = zero OE)
R.model=matrix(list(0),4,4)
R.model[1,1] = "nehalem"
R.model[2,2] = "siletz"
R.model[3,3] = "siuslaw"

### MARSS model 4b - unique trends for both populations (all run timing) since 1980, 
###   correlated process error. N.Umpqua OE fixed to zero
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4]),n,n)

mod_ORCoast4b.2 = MARSS(dat.ck4, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model),  
                        control=list(maxit=200000,conv.test.slope.tol = 0.09))

paramCK1<-MARSSparamCIs(mod_ORCoast4b.2)

###################################################################################################
### Columbia and Snake River populations 

Col_CK<-subset(CK_esc,group_1 >= 7 & group_1 < 13)
Col_pops<-subset(CK_pops, group_1 >=7 & group_1 <13)

unique(Col_CK$population)

Snake_CK1<-subset(Col_CK,group_1==10)
Snake_CK2<-subset(Col_CK,group_1==12)
Snake_CK<-rbind(Snake_CK1,Snake_CK2)

#########################################
####  Just Snake Spring/Summer stocks since 1957

Snake_CK<-subset(Snake_CK, year>=1957)

#### the data matrix
ck_dat3<-matrix(data=NA, nrow=length(seq(1957,2019,1)),ncol=length(unique(Snake_CK$pop_no))+1)
dimnames(ck_dat3)[[2]]<-c("year",unique(Snake_CK$population))

years<-seq(1957,2019,1)
ck_dat3[,1]<-years

pop.nos<-unique(Snake_CK$pop_no)

for(i in 1:length(pop.nos)){
  dat<-subset(Snake_CK,Snake_CK$pop_no==pop.nos[i])
  ck_dat3[which(ck_dat3[,1]==dat$year[1]):length(years),i+1]<-log(dat[,8])
}

## transposing the data for MARSS
dat.ck3 = t(ck_dat3)
years = dat.ck3[1,] # the years
n = nrow(dat.ck3)-1 # the number of time series
dat.ck3 = dat.ck3[2:nrow(dat.ck3),]

### MARSS model 3b - unique trends for each population, correlated process errors
Z.model = factor(c(seq(1,17,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5],QMat[1:n,6],
                   QMat[1:n,7],QMat[1:n,8],QMat[1:n,9],QMat[1:n,10],QMat[1:n,11],QMat[1:n,12],
                   QMat[1:n,13],QMat[1:n,14],QMat[1:n,15],QMat[1:n,16],QMat[1:n,17]),n,n)

mod_Col3b = MARSS(dat.ck3, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                  control=list(maxit=100000,conv.test.slope.tol = 0.09))


paramCK<-MARSSparamCIs(mod_Col3b)

######################################################################################
######   Analyzing subset of Columbia and Snake Chinook trends with data since 1967

pops67<-c(41,45,47,48,61,62,63,64,65,66,70,71,76,77,80)

CK67<-subset(Col_CK,year>1966)

CK67_1<-subset(CK67,CK67$pop_no==pops67[1])
CK67_2<-subset(CK67,CK67$pop_no==pops67[2])
CK67_3<-subset(CK67,CK67$pop_no==pops67[3])
CK67_4<-subset(CK67,CK67$pop_no==pops67[4])
CK67_5<-subset(CK67,CK67$pop_no==pops67[5])
CK67_6<-subset(CK67,CK67$pop_no==pops67[6])
CK67_7<-subset(CK67,CK67$pop_no==pops67[7])
CK67_8<-subset(CK67,CK67$pop_no==pops67[8])
CK67_9<-subset(CK67,CK67$pop_no==pops67[9])
CK67_10<-subset(CK67,CK67$pop_no==pops67[10])
CK67_11<-subset(CK67,CK67$pop_no==pops67[11])
CK67_12<-subset(CK67,CK67$pop_no==pops67[12])
CK67_13<-subset(CK67,CK67$pop_no==pops67[13])
CK67_14<-subset(CK67,CK67$pop_no==pops67[14])
CK67_15<-subset(CK67,CK67$pop_no==pops67[15])


CK67_dat<-rbind(CK67_1,CK67_2,CK67_3,CK67_4,CK67_5,CK67_6,CK67_7,CK67_8,CK67_9,
                CK67_10,CK67_11,CK67_12,CK67_13,CK67_14,CK67_15)
#####################################################################

ck_dat3.1<-matrix(data=NA, nrow=length(seq(1967,2019,1)),ncol=length(unique(CK67_dat$pop_no))+1)
dimnames(ck_dat3.1)[[2]]<-c("year",unique(CK67_dat$population))

years<-seq(1967,2019,1)
ck_dat3.1[,1]<-years

for(i in 1:length(pops67)){
  dat<-subset(CK67_dat,CK67_dat$pop_no==pops67[i])
  ck_dat3.1[which(ck_dat3.1[,1]==dat$year[1]):length(years),i+1]<-log(dat[,8])
}

## transposing the data for MARSS
dat.ck3.1 = t(ck_dat3.1)
years = dat.ck3.1[1,] # the years
n = nrow(dat.ck3.1)-1 # the number of time series
dat.ck3.1 = dat.ck3.1[2:nrow(dat.ck3.1),]

### MARSS model 3b - unique trends for each population, correlated process errors
Z.model = factor(c(seq(1,15,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5],QMat[1:n,6],
                   QMat[1:n,7],QMat[1:n,8],QMat[1:n,9],QMat[1:n,10],QMat[1:n,11],QMat[1:n,12],
                   QMat[1:n,13],QMat[1:n,14],QMat[1:n,15]),n,n)

mod_Col3.1b = MARSS(dat.ck3.1, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                    control=list(maxit=200000,conv.test.slope.tol = 0.09))

paramCK<-MARSSparamCIs(mod_Col3.1b)

######################################################################################
######   Columbia and Snake Chinook trends for subset of stocks with data since 1982
pops82s<-c(43,44,46,52,53,54,57,58)
CK82<-subset(Col_CK,year>1981)

CK82_1<-subset(CK82,CK82$pop_no==pops82s[1])
CK82_2<-subset(CK82,CK82$pop_no==pops82s[2])
CK82_3<-subset(CK82,CK82$pop_no==pops82s[3])
CK82_4<-subset(CK82,CK82$pop_no==pops82s[4])
CK82_5<-subset(CK82,CK82$pop_no==pops82s[5])
CK82_6<-subset(CK82,CK82$pop_no==pops82s[6])
CK82_7<-subset(CK82,CK82$pop_no==pops82s[7])
CK82_8<-subset(CK82,CK82$pop_no==pops82s[8])


CK82_dat<-rbind(CK82_1,CK82_2,CK82_3,CK82_4,CK82_5,CK82_6,CK82_7,CK82_8)

#####################################################################

ck_dat3.2<-matrix(data=NA, nrow=length(seq(1982,2019,1)),ncol=length(unique(CK82_dat$pop_no))+1)
dimnames(ck_dat3.2)[[2]]<-c("year",unique(CK82_dat$population))

years<-seq(1982,2019,1)
ck_dat3.2[,1]<-years

for(i in 1:length(pops82s)){
  dat<-subset(CK82_dat,CK82_dat$pop_no==pops82s[i])
  ck_dat3.2[which(ck_dat3.2[,1]==dat$year[1]):length(years),i+1]<-log(dat[,8])
}

## transposing the data for MARSS
dat.ck3.2 = t(ck_dat3.2)
years = dat.ck3.2[1,] # the years
n = nrow(dat.ck3.2)-1 # the number of time series
dat.ck3.2 = dat.ck3.2[2:nrow(dat.ck3.2),]


######### Setting the R matrix to zero for dam counts

## Setting observer error for dam counts to zero
R.model=matrix(list(0),n,n)
R.model[1,1] = 0
R.model[2,2] = 0
R.model[3,3] = "wenatchee_su"
R.model[4,4] = 0
R.model[5,5] = "sandy"
R.model[6,6] = "mckenzie"
R.model[7,7] = "deschutes_fa"
R.model[8,8] = "deschutes_sp"


### MARSS model 3b - unique trends for each population, correlated process errors
inits=(list(R=0.05))
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5],QMat[1:n,6],
                   QMat[1:n,7],QMat[1:n,8]),n,n)

#mod_Col3.2b = MARSS(dat.ck3.2, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
#                    control=list(maxit=100000,conv.test.slope.tol = 0.2),inits=inits)

mod_Col3.2bs = MARSS(dat.ck3.2, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                     control=list(maxit=100000, conv.test.slope.tol=0.09),inits=inits)

paramCK<-MARSSparamCIs(mod_Col3.2bs)

##############################################################################################
#### Klamath and Rogue Chinook escapement trends
RoKlam_CK<-subset(CK_esc,group_1 == 14)

### starting at
RoKlam_CK<-subset(RoKlam_CK, RoKlam_CK$year>1977)
## dropping shasta, scott, and salmon fall - they're subsets of total Klamath Fall
RoKlam_CK<-RoKlam_CK[-which(RoKlam_CK$population=="shasta_fa"),]
RoKlam_CK<-RoKlam_CK[-which(RoKlam_CK$population=="scott_fa"),]
RoKlam_CK<-RoKlam_CK[-which(RoKlam_CK$population=="salmon_fa"),]

ck_dat5<-matrix(data=NA, nrow=length(seq(1978,2019,1)),ncol=length(unique(RoKlam_CK$pop_no))+1)
dimnames(ck_dat5)[[2]]<-c("year",unique(RoKlam_CK$population))

years<-seq(1978,2019,1)
ck_dat5[,1]<-years

pop.nos<-unique(RoKlam_CK$pop_no)

for(i in 1:length(pop.nos)){
  dat<-subset(RoKlam_CK,RoKlam_CK$pop_no==pop.nos[i])
  ck_dat5[which(ck_dat5[,1]==dat$year[1]):length(years),i+1]<-log(dat[,8])
}

## transposing the data for MARSS
dat.ck5 = t(ck_dat5)
years = dat.ck5[1,] # the years
n = nrow(dat.ck5)-1 # the number of time series
dat.ck5 = dat.ck5[2:nrow(dat.ck5),]

### MARSS model 5b - unique trends for both populations (all run timing) since 1980 - correlated process error
inits=(list(R=0.1))
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5],
                   QMat[1:n,6],QMat[1:n,7]),n,n)

mod_RoKlam5b = MARSS(dat.ck5, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                     control=list(maxit=200000,conv.test.slope.tol = 0.09), inits=inits)

#### exporting parameter estimates & population states
paramCK5b<-MARSSparamCIs(mod_RoKlam5b)

###################################################################################################
#### Sacramento Chinook escapement trends
Sac_CK<-subset(CK_esc,group_1 >= 15)

### starting at 1960
ck_dat6<-matrix(data=NA, nrow=length(seq(1960,2019,1)),ncol=length(unique(Sac_CK$pop_no))+1)
dimnames(ck_dat6)[[2]]<-c("year",unique(Sac_CK$population))

years<-seq(1960,2019,1)
ck_dat6[,1]<-years

pop.nos<-unique(Sac_CK$pop_no)

for(i in 1:length(pop.nos)){
  dat<-subset(Sac_CK,Sac_CK$pop_no==pop.nos[i])
  ck_dat6[which(ck_dat6[,1]==dat$year[1]):length(years),i+1]<-log(dat[,8])
}

## transposing the data for MARSS
dat.ck6= t(ck_dat6)
years = dat.ck6[1,] # the years
n = nrow(dat.ck6)-1 # the number of time series
dat.ck6 = dat.ck6[2:nrow(dat.ck6),]

### MARSS model 6b - unique trends for both populations (all run timing) since 1960 - correlated process error
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5]),n,n)

mod_Sac6b = MARSS(dat.ck6, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                  control=list(maxit=200000,conv.test.slope.tol = 0.09))

#### exporting parameter estimates & population states
paramCK6b<-MARSSparamCIs(mod_Sac6b)


######################################################################################
####          Total Run Size Analysis - 1990 to 2019
####################################################################################

CK_run<-read.table("CK_TotalRun_FINAL.txt", header=TRUE)
CK_run1<-CK_run[which(CK_run$year>1989),]
QMat<-read.csv("corr_mat_CK.csv",header=FALSE)

#######################################################
###   Total run size - regional groups

### SALISH
CK_Salish<-subset(CK_run1,CK_run1$group=="salish")
pop.nos<-unique(CK_Salish$pop_no)

ck_runS<-matrix(data=NA, nrow=length(seq(1990,2019,1)),ncol=length(unique(CK_Salish$population))+1)
dimnames(ck_runS)[[2]]<-c("year",unique(CK_Salish$population))

years<-seq(1990,2019,1)
ck_runS[,1]<-years

for(i in 1:length(pop.nos)){
  dat<-subset(CK_Salish,CK_Salish$pop_no==pop.nos[i])
  ck_runS[which(ck_runS[,1]==dat$year[1]):length(years),i+1]<-log(dat[,9])
}

## transposing the data for MARSS
ck_runS = t(ck_runS)
years = ck_runS[1,] # the years
n = nrow(ck_runS)-1 # the number of time series
ck_runS = ck_runS[2:nrow(ck_runS),]

### MARSS model total runsize - unique trends for each population, correlated process errors
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5],QMat[1:n,6],
                   QMat[1:n,7]),n,n)

mod_CK_Salish.run = MARSS(ck_runS, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                          control=list(maxit=200000,conv.test.slope.tol = 0.09))

paramCK<-MARSSparamCIs(mod_CK_Salish.run)

#################################################
### COLUMBIA
CK_Col<-subset(CK_run1,CK_run1$group=="columbia")
pop.nos<-unique(CK_Col$pop_no)

ck_runC<-matrix(data=NA, nrow=length(seq(1990,2019,1)),ncol=length(unique(CK_Col$population))+1)
dimnames(ck_runC)[[2]]<-c("year",unique(CK_Col$population))

years<-seq(1990,2019,1)
ck_runC[,1]<-years

for(i in 1:length(pop.nos)){
  dat<-subset(CK_Col,CK_Col$pop_no==pop.nos[i])
  ck_runC[which(ck_runC[,1]==dat$year[1]):length(years),i+1]<-log(dat[,9])
}

## transposing the data for MARSS
ck_runC = t(ck_runC)
years = ck_runC[1,] # the years
n = nrow(ck_runC)-1 # the number of time series
ck_runC = ck_runC[2:nrow(ck_runC),]

### MARSS model total runsize - unique trends for each population, correlated process errors
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5],QMat[1:n,6],
                   QMat[1:n,7],QMat[1:n,8]),n,n)

mod_CK_Col.run = MARSS(ck_runC, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                       control=list(maxit=200000,conv.test.slope.tol = 0.09))

paramCK2<-MARSSparamCIs(mod_CK_Col.run)

#######################################
### WA & OR Coastal
CK_Coast<-rbind(subset(CK_run1,CK_run1$group=="wa_coast"),subset(CK_run1,CK_run1$group=="or_coast"))
pop.nos<-unique(CK_Coast$pop_no)

ck_run.Coast<-matrix(data=NA, nrow=length(seq(1990,2019,1)),ncol=length(unique(CK_Coast$population))+1)
dimnames(ck_run.Coast)[[2]]<-c("year",unique(CK_Coast$population))

years<-seq(1990,2019,1)
ck_run.Coast[,1]<-years

for(i in 1:length(pop.nos)){
  dat<-subset(CK_Coast,CK_Coast$pop_no==pop.nos[i])
  ck_run.Coast[which(ck_run.Coast[,1]==dat$year[1]):length(years),i+1]<-log(dat[,9])
}

## transposing the data for MARSS
ck_runC = t(ck_run.Coast)
years = ck_runC[1,] # the years
n = nrow(ck_runC)-1 # the number of time series
ck_runC = ck_runC[2:nrow(ck_runC),]

### MARSS model total runsize - unique trends for each population, correlated process errors
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3],QMat[1:n,4],QMat[1:n,5]),n,n)

mod_CK_Coast.run = MARSS(ck_runC, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                         control=list(maxit=200000,conv.test.slope.tol = 0.09, safe=TRUE))

paramCK3<-MARSSparamCIs(mod_CK_Coast.run)

###################################################
### California
CK_Cali<-subset(CK_run1,CK_run1$group=="california")
pop.nos<-unique(CK_Cali$pop_no)

ck_run.Cali<-matrix(data=NA, nrow=length(seq(1990,2019,1)),ncol=length(unique(CK_Cali$population))+1)
dimnames(ck_run.Cali)[[2]]<-c("year",unique(CK_Cali$population))

years<-seq(1990,2019,1)
ck_run.Cali[,1]<-years

for(i in 1:length(pop.nos)){
  dat<-subset(CK_Cali,CK_Cali$pop_no==pop.nos[i])
  ck_run.Cali[which(ck_run.Cali[,1]==dat$year[1]):length(years),i+1]<-log(dat[,9])
}

## transposing the data for MARSS
ck_runC = t(ck_run.Cali)
years = ck_runC[1,] # the years
n = nrow(ck_runC)-1 # the number of time series
ck_runC = ck_runC[2:nrow(ck_runC),]

### MARSS model total runsize - unique trends for each population, correlated process errors
Z.model = factor(c(seq(1,n,1))) # Each population has it's own trend
R.model = "diagonal and unequal" # observation error is independent and uncorrelated
u.model = "unequal" # each population has it's own long-term trend
Q.model = matrix(c(QMat[1:n,1],QMat[1:n,2],QMat[1:n,3]),n,n)

mod_CK_Cali.run = MARSS(ck_runC, model=list(Z=Z.model, R=R.model, U=u.model, Q=Q.model), 
                        control=list(maxit=200000,conv.test.slope.tol = 0.09, safe=TRUE))

paramCK4<-MARSSparamCIs(mod_CK_Cali.run)
