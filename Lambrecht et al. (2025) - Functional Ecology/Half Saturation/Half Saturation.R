setwd("~/Activité Professionnelle/LIMNO 2019-2023/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Half Saturation")

rm(list=ls())

library(car)
library(cowplot)
library(data.table)
library(DescTools)
library(deSolve)
library(dplyr)
library(foreach)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
library(lattice)
library(lme4)
library(lmtest)
library(magrittr)
library(nlme)
library(plyr)
library(plotly)
library(propagate)
library(qpcR)
library(reshape2)
library(rstatix)
library(scales)

################################################
################################################
##### HALF SATURATION CONSTANT OF ISOLATES #####
################################################
################################################

# Import the dataset
Data=read.table("Data_HS.txt", h=T, dec=".")
Data$Species=gsub("\\_","\\ ", Data$Species)
summary(Data)
names(Data)

# Specify the variables as numeric or factor
Data[,c(6:9,11:12)] %<>% mutate_if(is.character,as.numeric)
Data[,c(6:9,11:12)] %<>% mutate_if(is.integer,as.numeric)
Data[,c(6:9,11:12)] %<>% mutate_if(is.factor,as.numeric)

# Sort the dataset by isolate and nutrient treatments
Data$Iso=factor(Data$Iso, levels=unique(Data$Iso))
Data$Species=factor(Data$Species, levels=unique(Data$Species))
Data$Phos=factor(Data$Phos, levels=unique(Data$Phos))
Data=Data[order(Data$Iso,Data$Phos,Data$Trial),]

# Calculate the densities
Data$Dens=Data$IXM*Data$Volu*Data$Site*Data$Dilu*Data$Cove

# Remove the outlier densities
Data=as.data.frame(Data %>% group_by(Iso,Species,Phos,Day) %>% mutate(Dens=ifelse(Dens < mean(Dens)-2.0*sd(Dens), NA, Dens)))
Data=as.data.frame(Data %>% group_by(Iso,Species,Phos,Day) %>% mutate(Dens=ifelse(Dens > mean(Dens)+2.0*sd(Dens), NA, Dens)))

# Remove the decreasing densities
Data=Data[Data$Day > 0,]
MinDay=as.data.frame(subset(na.omit(Data), Day < 4) %>% group_by(Iso,Species,Phos,Trial) %>% slice(which.min(Dens)))$Day
Data$MinDay=rep(MinDay, each=length(unique(Data$Day)))
Data=Data %>% group_by(Iso,Species,Phos,Trial) %>% filter(Day >= MinDay) %>% ungroup

# Calculate the mean densities
Data2=setDT(Data)[, .(MeanDens=round(mean(Dens),0)), by=list(Iso,Species,Phos,Day)]
Data=as.data.frame(Data)
Data2=as.data.frame(Data2)


##########################
### Growth rate models ###
##########################

### Linear model###

# Function for parameters
SetLi=function(mCall,LHS,data) {
  xy=sortedXyData(mCall[["x"]],LHS,data)
  b=min(xy$y); a=coef(lm(y~x,xy))[2] 
  Out=c(b,a); names(Out)=mCall[c("b","a")]
  return(Out)}
SSline=selfStart(as.formula("~b + a*x"), initial=SetLi, parameters=c("b","a"))

# Function for fitting
FuncLi=function(Model, Day) {
  Params=coef(Model)
  names(Params)=NULL
  b=Params[1]
  a=Params[2]
  
  Parameters=data.frame(b=b,a=a)
  Rates=data.frame(
    Day = Day,
    DensP = b + a * Day,
    AGR = a)
  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Parameters=Parameters, Rates=Rates)
  return(Out)
}

### Logistic model ###

# Function for parameters
SetLo=function(Coef){
  K=Coef[1]; r=1/(Coef[3]); b=K/(1 + exp(Coef[2]/Coef[3]))
  Out=c(b,K,r); names(Out)=c("b","K","r")
  return(Out)}

# Function for fitting
FuncLo=function(Model,Day) {
  Coef=coef(Model)
  Params=SetLo(Coef)
  b=Params[1]
  K=Params[2]
  r=Params[3]
  
  Parameters=data.frame(b=b,K=K,r=r)
  Rates=data.frame(
    Day = Day,
    DensP = (b * K) / (b + (K - b) * exp(-r * Day)),
    AGR = (r * b * K * (K - b) * exp(-r * Day)) / (b + (K - b)*exp(-r * Day))^2)
  
  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Parameters=Parameters, Rates=Rates)
  return(Out)
}


##################################
### Fitting growth rate models ###
##################################

# Split the datasets
SplitData=split(Data, list(Data$Iso, Data$Phos))
SplitData=SplitData[sapply(SplitData, function(x) dim(x)[1]) > 0]

# Clean the datasets
SplitData[[11]]$Dens[4]=NA
SplitData[[11]]$Dens[9]=NA
SplitData[[11]]$Dens[14]=NA
SplitData[[12]]$Dens[9]=NA
SplitData[[14]]$Dens[3]=NA
SplitData[[14]]$Dens[7]=NA
SplitData[[14]]$Dens[11]=NA
SplitData[[17]]$Dens[1]=NA
SplitData[[19]]$Dens[9]=NA
SplitData[[19]]$Dens[14]=NA
SplitData[[38]]$Dens[7]=NA
SplitData[[40]]$Dens[1]=NA
SplitData[[41]]$Dens[1]=NA
SplitData[[41]]$Dens[5]=NA
SplitData[[48]]$Dens[2]=NA
SplitData[[48]]$Dens[11]=NA
SplitData[[49]]$Dens[3]=NA
SplitData[[49]]$Dens[7]=NA
SplitData[[52]]$Dens[1]=NA
SplitData[[60]]$Dens[9]=NA
SplitData[[63]]$Dens[3]=NA
SplitData[[73]]$Dens[4]=NA
SplitData[[73]]$Dens[7]=NA
SplitData[[73]]$Dens[9]=NA
SplitData[[73]]$Dens[14]=NA
SplitData[[82]]$Dens[4]=NA
SplitData[[82]]$Dens[9]=NA
SplitData[[87]]$Dens[9]=NA
SplitData[[89]]$Dens[12]=NA
SplitData[[98]]$Dens[2]=NA
SplitData[[99]]$Dens[5]=NA
SplitData[[99]]$Dens[13]=NA
SplitData[[101]]$Dens[9]=NA
SplitData[[105]]$Dens[1]=NA
SplitData[[105]]$Dens[2]=NA
SplitData[[105]]$Dens[5]=NA
SplitData[[108]]$Dens[1]=NA
SplitData[[109]]$Dens[14]=NA
SplitData[[122]]$Dens[7]=NA
SplitData[[128]]$Dens[5]=NA
SplitData[[128]]$Dens[6]=NA
SplitData[[128]]$Dens[13]=NA
SplitData[[132]]$Dens[6]=NA
SplitData[[138]]$Dens[4]=NA
SplitData[[138]]$Dens[9]=NA
SplitData[[138]]$Dens[14]=NA
SplitData[[147]]$Dens[3]=NA
SplitData[[147]]$Dens[7]=NA
SplitData[[147]]$Dens[12]=NA
SplitData[[148]]$Dens[4]=NA
SplitData[[148]]$Dens[12]=NA
SplitData[[153]]$Dens[1]=NA
SplitData[[154]]$Dens[1]=NA
SplitData[[154]]$Dens[5]=NA
SplitData[[158]]$Dens[4]=NA
SplitData[[158]]$Dens[9]=NA
SplitData[[158]]$Dens[14]=NA
SplitData[[161]]$Dens[4]=NA
SplitData[[161]]$Dens[11]=NA
SplitData[[168]]$Dens[10]=NA
SplitData[[170]]$Dens[1]=NA
SplitData[[170]]$Dens[5]=NA
SplitData[[171]]$Dens[11]=NA
SplitData[[174]]$Dens[2]=NA
SplitData[[185]]$Dens[5]=NA
SplitData[[185]]$Dens[9]=NA
SplitData[[186]]$Dens[1]=NA
SplitData[[187]]$Dens[4]=NA
SplitData[[187]]$Dens[9]=NA
SplitData[[187]]$Dens[13]=NA
SplitData[[192]]$Dens[6]=NA
SplitData[[197]]$Dens[2]=NA
SplitData[[197]]$Dens[3]=NA
SplitData[[197]]$Dens[7]=NA
SplitData[[198]]$Dens[9]=NA
SplitData[[200]]$Dens[6]=NA
SplitData[[206]]$Dens[1]=NA
SplitData[[219]]$Dens[1]=NA
SplitData[[219]]$Dens[6]=NA
SplitData[[227]]$Dens[4]=NA
SplitData[[227]]$Dens[9]=NA
SplitData[[234]]$Dens[4]=NA
SplitData[[234]]$Dens[5]=NA
SplitData[[241]]$Dens[11]=NA
SplitData[[244]]$Dens[7]=NA
SplitData[[245]]$Dens[2]=NA
SplitData[[246]]$Dens[9]=NA
SplitData[[246]]$Dens[11]=NA
SplitData[[251]]$Dens[8]=NA
SplitData[[253]]$Dens[12]=NA
SplitData[[258]]$Dens[4]=NA
SplitData[[258]]$Dens[10]=NA
SplitData[[259]]$Dens[3]=NA
SplitData[[259]]$Dens[7]=NA
SplitData[[265]]$Dens[5]=NA
SplitData[[265]]$Dens[10]=NA
SplitData[[266]]$Dens[1]=NA
SplitData[[266]]$Dens[5]=NA
SplitData[[266]]$Dens[9]=NA
SplitData[[266]]$Dens[13]=NA
SplitData[[270]]$Dens[10]=NA
SplitData[[275]]$Dens[1]=NA
SplitData[[275]]$Dens[10]=NA
SplitData[[279]]$Dens[3]=NA
SplitData[[279]]$Dens[4]=NA
SplitData[[283]]$Dens[3]=NA
SplitData[[284]]$Dens[3]=NA
SplitData[[286]]$Dens[3]=NA
SplitData[[286]]$Dens[11]=NA
SplitData[[293]]$Dens[2]=NA
SplitData[[293]]$Dens[5]=NA
SplitData[[293]]$Dens[6]=NA
SplitData[[297]]$Dens[11]=NA
SplitData[[298]]$Dens[4]=NA
SplitData[[301]]$Dens[5]=NA
SplitData[[304]]$Dens[9]=NA
SplitData[[305]]$Dens[6]=NA
SplitData[[310]]$Dens[15]=NA
SplitData[[311]]$Dens[15]=NA
SplitData[[314]]$Dens[8]=NA
SplitData[[314]]$Dens[9]=NA
SplitData[[317]]$Dens[5]=NA
SplitData[[320]]$Dens[11]=NA
SplitData[[322]]$Dens[14]=NA
SplitData[[325]]$Dens[2]=NA
SplitData[[325]]$Dens[2]=NA
SplitData[[325]]$Dens[7]=NA
SplitData[[328]]$Dens[1]=NA
SplitData[[332]]$Dens[5]=NA
SplitData[[332]]$Dens[6]=NA
SplitData[[333]]$Dens[1]=NA
SplitData[[336]]$Dens[2]=NA
SplitData[[340]]$Dens[11]=NA
SplitData[[340]]$Dens[12]=NA
SplitData[[342]]$Dens[1]=NA

# Separate the datasets
SplitDataLi=SplitData[c(1,9,10,14,34,38,39,40,49,87,88,95,97,98,104,136,137,153,185,186,195,196,209,234,235,236,242,244,245,246,247,248,254,255,283,285,295,296,300,304,334)]
SplitDataLo=SplitData[-c(1,9,10,14,34,38,39,40,49,87,88,95,97,98,104,136,137,153,185,186,195,196,209,234,235,236,242,244,245,246,247,248,254,255,283,285,295,296,300,304,334)]

### Linear model ###

# Extract the combinations of names
Iso=lapply(SplitDataLi,"[",c("Iso"))
Species=lapply(SplitDataLi,"[",c("Species"))
Phos=lapply(SplitDataLi,"[",c("Phos"))
Iso=(cbind(do.call("rbind",Iso)))
Species=(cbind(do.call("rbind",Species)))
Phos=(cbind(do.call("rbind",Phos)))

# Fitting the linear model
ModLi=function(x) {
  FitLi=nls(log(Dens+1) ~ SSline(Day, a, b), data=x)
  OutLi=FuncLi(FitLi, x$Day)
} 
OutLi=lapply(SplitDataLi, ModLi)

# Create a dataset
RateLi=bind_rows(lapply(OutLi, function (x) x[c("Rates")]))
RateLi=round(cbind(do.call("rbind",RateLi)),4)
RateLi=cbind(Iso=Iso,Species=Species,Phos=Phos,Model="Linear",RateLi)
rownames(RateLi)=c()

### Logistic model ###

# Extract combinations of names
Iso=lapply(SplitDataLo,"[",c("Iso"))
Species=lapply(SplitDataLo,"[",c("Species"))
Phos=lapply(SplitDataLo,"[",c("Phos"))
Iso=(cbind(do.call("rbind",Iso)))
Species=(cbind(do.call("rbind",Species)))
Phos=(cbind(do.call("rbind",Phos)))

# Logistic model
ModLo=function(x) {
  FitLo=nls(log(Dens+1) ~ SSlogis(Day, Asym, xmid, scal), data=x)
  OutLo=FuncLo(FitLo, x$Day)
}
lapply(SplitDataLo[54], ModLo)
OutLo=lapply(SplitDataLo, ModLo)

# Create a dataset
RateLo=bind_rows(lapply(OutLo, function (x) x[c("Rates")]))
RateLo=round(cbind(do.call("rbind",RateLo)),4)
RateLo=cbind(Iso=Iso,Species=Species,Phos=Phos,Model="Logistic",RateLo)
rownames(RateLo)=c()


################################
### Calculating growth rates ###
################################

# Combine the datasets
Data3=rbind(RateLi,RateLo)
Data3[,c(6)]=exp(Data3[,c(6)])
Data3[,c(6)]=round(Data3[,c(6)],4)

# Add the initial density
Data3=Data3 %>% group_by(Iso,Species,Phos) %>% do(add_row(., Iso=unique(.$Iso), Species=unique(.$Species), Phos=unique(.$Phos), Model=unique(.$Model), Day=0, DensP=log(5000), AGR=0, RGR=0))
Data3=Data3[order(Data3$Iso, Data3$Phos, Data3$Day),]
Data3=Data3[!duplicated(Data3),]
Data3=as.data.frame(Data3)

# Order the dataset
Data=Data[order(Data$Phos,Data$Iso,Data$Trial),]
Data$Iso=factor(Data$Iso, levels=unique(Data$Iso))
Data$Phos=factor(Data$Phos, levels=unique(Data$Phos))
Data$Trial=factor(Data$Trial, levels=unique(Data$Trial))

# Split the datasets
SplitData3a=split(Data, list(Data$Iso, Data$Phos))
SplitData3b=split(Data, list(Data$Iso, Data$Phos, Data$Trial))
SplitData3a=SplitData3a[sapply(SplitData3a, function(x) dim(x)[1]) > 0]
SplitData3b=SplitData3b[sapply(SplitData3b, function(x) dim(x)[1]) > 0]

# Calculate the maximum growth rates
ModGR=function(x) {summary(lm(log(Dens+1)~Day, data=x))$coefficients[2]}
ModGRSD=function(x) {summary(lm(log(Dens+1)~Day, data=x))$coefficients[4]}
GRInter=lapply(SplitData3a, ModGR)
GRInterSD=lapply(SplitData3a, ModGRSD)
GRIntra=lapply(SplitData3b, ModGR)

# Create an interspecific dataset
Iso=c(unique(Data[,c("Iso", "Species","Phos")])[,1])
Species=c(unique(Data[,c("Iso", "Species","Phos")])[,2])
Phos=c(unique(Data[,c("Iso", "Species","Phos")])[,3])
GRInter=as.data.frame(round(cbind(do.call("rbind",GRInter)),4))[[1]]
GRInterSD=as.data.frame(round(cbind(do.call("rbind",GRInterSD)),4))[[1]]
Data4=data.frame(Iso=Iso, Species=Species, Phos=Phos, GR=GRInter, GRLSD=GRInter-GRInterSD, GRUSD=GRInter+GRInterSD)

# Create an intraspecific dataset
Iso=rep(c(unique(Data[,c("Iso", "Species","Phos")])[,1]),3)
Species=rep(c(unique(Data[,c("Iso", "Species","Phos")])[,2]),3)
Phos=rep(c(unique(Data[,c("Iso", "Species","Phos")])[,3]),3)
Trial=rep(c(unique(Data[,c("Trial")])),each=343)
GRIntra=as.data.frame(round(cbind(do.call("rbind",GRIntra)),4))[[1]]
Data5=data.frame(Iso=Iso, Species=Species, Phos=Phos, Trial=Trial, GR=GRIntra)

# Adding the null growth rate
Data4=as.data.frame(Data4 %>% group_by(Iso,Species) %>% do(add_row(., Iso=unique(.$Iso), Species=unique(.$Species), Phos="0", GR=0, GRLSD=0, GRUSD=0)))
Data5=as.data.frame(Data5 %>% group_by(Iso,Species) %>% do(add_row(., Iso=unique(.$Iso), Species=unique(.$Species), Phos="0", Trial=unique(.$Trial), GR=0)))
Data4$Phos=as.numeric(as.character(Data4$Phos))
Data5$Phos=as.numeric(as.character(Data5$Phos))
Data4=Data4[order(Data4$Iso, Data4$Phos),]
Data5=Data5[order(Data5$Trial, Data5$Iso, Data5$Phos),]


############################
### Fitting Monod models ###
############################

### Growth rate ###

# Remove negative growth rates
Data4$GR[Data4$GR<0]=0
Data5$GR[Data5$GR<0]=0

# Prepare the datasets
Data4$GRM=rep(as.data.frame(Data4 %>% group_by(Iso,Species) %>% summarise(GRM=max(GR)))[,3], each=8)
Data4$GRMLSD=rep(as.data.frame(Data4 %>% group_by(Iso,Species) %>% summarise(GRM=max(GRLSD)))[,3], each=8)
Data4$GRMUSD=rep(as.data.frame(Data4 %>% group_by(Iso,Species) %>% summarise(GRM=max(GRUSD)))[,3], each=8)
Data5$GRM=rep(as.data.frame(Data5 %>% group_by(Trial,Iso,Species) %>% summarise(GRM=max(GR)))[,4], each=8)
GRMInter=subset(Data4, Phos=="0")$GRM
GRMInterLSD=subset(Data4, Phos=="0")$GRMLSD
GRMInterUSD=subset(Data4, Phos=="0")$GRMUSD
GRMIntra=subset(Data5, Phos=="0")$GRM

# Split the datasets
SplitData4a=split(Data4, list(Data4$Iso))
SplitData4b=split(Data5, list(Data5$Iso, Data5$Trial))

# Monod model
FuncH=function(x) {OutH=summary(nls(GR~(GRM*Phos)/(Phos + K), start=c(K=1), data=x))$coefficients[2]}
FuncHSD=function(x) {OutH=summary(nls(GR~(GRM*Phos)/(Phos + K), start=c(K=1), data=x))$coefficients[4]}
HSInter=lapply(SplitData4a, FuncH)
HSInterSD=lapply(SplitData4a, FuncHSD)
HSIntra=lapply(SplitData4b, FuncH)

# Create an interspecific dataset
Iso=c(unique(Data[,c("Iso","Species")])[,1])
Species=c(unique(Data[,c("Iso","Species")])[,2])
HSInter=as.data.frame(round(cbind(do.call("rbind",HSInter)),4))[[1]] 
HSInterSD=as.data.frame(round(cbind(do.call("rbind",HSInterSD)),4))[[1]]
Data6=data.frame(Iso=Iso, Species=Species, HS=HSInter, HSLSD=HSInter-HSInterSD, HSUSD=HSInter+HSInterSD, AF=round(GRMInter/HSInter,4), AFLSD=round(GRMInterLSD/HSInter,4), AFUSD=round(GRMInterUSD/HSInter,4), GRM=GRMInter, GRMLSD=GRMInterLSD, GRMUSD=GRMInterUSD)

# Create an intraspecific dataset
Iso=rep(c(unique(Data[,c("Iso","Species")])[,1]),3)
Species=rep(c(unique(Data[,c("Iso","Species")])[,2]),3)
Trial=rep(c(unique(Data[,c("Trial")])),each=49)
HSIntra=as.data.frame(round(cbind(do.call("rbind",HSIntra)),4))[[1]] 
Data7=data.frame(Iso=Iso, Species=Species, Trial=Trial, HS=HSIntra, AF=round(GRMIntra/HSIntra,4), GRM=GRMIntra)

# Export the datasets
Data6$HS[Data6$HS < 0]=0
Data6$HSLSD[Data6$HSLSD < 0]=0
Data6$HSUSD[Data6$HSUSD < 0]=0
Data7$HS[Data7$HS < 0]=0
Data7=Data7[order(Data7$Iso,Data7$Trial),]
write.table(Data6[,c(1:11)], file="Data_HS_Inter.txt", sep="\t", row.names=F)
write.table(Data7[,c(1:6)], file="Data_HS_Intra.txt", sep="\t", row.names=F)

# Monod model
FuncH=function(x) {OutH=predict(nls(GR~(GRM*Phos)/(Phos + K), start=c(K=1), data=x))}
HSInter=lapply(SplitData4a, FuncH)

# Create a dataset
Iso=c(unique(Data4[,c("Iso", "Species","Phos")])[,1])
Species=c(unique(Data4[,c("Iso", "Species","Phos")])[,2])
Phos=c(unique(Data4[,c("Iso", "Species","Phos")])[,3])
Data8=data.frame(Iso=Iso, Species=Species, Phos=Phos, GRP=unlist(HSInter))

### Error bars ###

# Calculate the error bars for half saturation constant
Data10=as.data.frame(setDT(na.omit(Data6))[, .(HSSpeciesL=mean(HS)-sd(HS), HSSpecies=mean(HS), HSSpeciesU=mean(HS)+sd(HS)), by=list(Species)])
Data10$HSSpecies[Data10$HSSpecies < 0]=0
Data10$HSSpeciesL[Data10$HSSpeciesL < 0]=0
Data10$HSSpeciesU[Data10$HSSpeciesU < 0]=0
Data6$HSSpecies=c(rep(Data10$HSSpecies[1],24), rep(Data10$HSSpecies[2],25))
Data6$HSSpeciesL=c(rep(Data10$HSSpeciesL[1],24), rep(Data10$HSSpeciesL[2],25))
Data6$HSSpeciesU=c(rep(Data10$HSSpeciesU[1],24), rep(Data10$HSSpeciesU[2],25))

# Calculate the error bars for affinity constant
Data11=as.data.frame(setDT(na.omit(Data6))[, .(AFSpeciesL=mean(AF)-sd(AF), AFSpecies=mean(AF), AFSpeciesU=mean(AF)+sd(AF)), by=list(Species)])
Data11$AFSpeciesL[Data11$AFSpeciesL < 0]=0
Data11$AFSpecies[Data11$AFSpecies < 0]=0
Data11$AFSpeciesU[Data11$AFSpeciesU < 0]=0
Data6$AFSpeciesL=c(rep(Data11$AFSpeciesL[1],24), rep(Data11$AFSpeciesL[2],25))
Data6$AFSpecies=c(rep(Data11$AFSpecies[1],24), rep(Data11$AFSpecies[2],25))
Data6$AFSpeciesU=c(rep(Data11$AFSpeciesU[1],24), rep(Data11$AFSpeciesU[2],25))

# Calculate the error bars for maximum growth rate
Data12=as.data.frame(setDT(na.omit(Data6))[, .(GRMSpeciesL=mean(GRM)-sd(GRM), GRMSpecies=mean(GRM), GRMSpeciesU=mean(GRM)+sd(GRM)), by=list(Species)])
Data12$GRMSpeciesL[Data12$GRMSpeciesL < 0]=0
Data12$GRMSpecies[Data12$GRMSpecies < 0]=0
Data12$GRMSpeciesU[Data12$GRMSpeciesU < 0]=0
Data6$GRMSpeciesL=c(rep(Data12$GRMSpeciesL[1],24), rep(Data12$GRMSpeciesL[2],25))
Data6$GRMSpecies=c(rep(Data12$GRMSpecies[1],24), rep(Data12$GRMSpecies[2],25))
Data6$GRMSpeciesU=c(rep(Data12$GRMSpeciesU[1],24), rep(Data12$GRMSpeciesU[2],25))

# Export the dataset
Data6[,c(3:20)]=round(Data6[,c(3:20)],4)
write.table(Data6[,c(1:20)], file="Data_HS_TO.txt", sep="\t", row.names=F)


###########################################
### Plot predicted saturation constants ###
###########################################

tiff('Half Saturation Curves.tiff', units="in", width=30, height=30, res=1000)
ggplot(Data8, aes(Phos, GRP, group=Species)) + 
  geom_line(aes(color=Species),size=1) +
  geom_point(aes(color=Species), alpha=0.8, size=3, pch=16) +
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression('Phosphorus concentration'~'('*µµmol~PO[4]^'3-'~mL^-1*')'))+
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.5,by=0.1), limits=c(0,0.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,3.2,by=0.8), limits=c(0,3.2)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.text.x=element_text(face="plain", color="black", size=18)) +
  theme(strip.background=element_blank()) +
  facet_wrap(~Iso, ncol=7, nrow=7, scales="free") +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="bottom")
dev.off()

tiff('Half Saturation Constants.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data6, aes(Species, HS, group=Species)) +
  geom_pointrange(aes(Species, ymin=HSLSD, ymax=HSUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=0.6+(0.020*0.6), yend=0.6+(0.020*0.6), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=0.6+(0.020*0.6), yend=0.6+(0.020*0.6), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=0.6+(0.100*0.6), yend=0.6+(0.100*0.6), color="black", size=0.8) +
  geom_text(x=1.00, y=0.6+(0.035*0.6), label="**", color="black", size=5) +
  geom_text(x=2.00, y=0.6+(0.035*0.6), label="***", color="black", size=5) +
  geom_text(x=1.50, y=0.6+(0.125*0.6), label="NS", color="black", size=5) +
  ylab(expression('Half-saturation constant'~'('*µmol~PO[4]^{"-"}~mL^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.6,by=0.2), limits=c(-10^-5,0.6+(0.10*0.6))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off() 

tiff('Maximum Growth Rates.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data6, aes(Species, GRM, group=Species)) +
  geom_pointrange(aes(Species, ymin=GRMLSD, ymax=GRMUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=0.6+(0.020*0.6), yend=0.6+(0.020*0.6), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=0.6+(0.020*0.6), yend=0.6+(0.020*0.6), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=0.6+(0.100*0.6), yend=0.6+(0.100*0.6), color="black", size=0.8) +
  geom_text(x=1.00, y=0.6+(0.035*0.6), label="***", color="black", size=5) +
  geom_text(x=2.00, y=0.6+(0.035*0.6), label="***", color="black", size=5) +
  geom_text(x=1.50, y=0.6+(0.125*0.6), label="NS", color="black", size=5) +
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.6,by=0.2), limits=c(-10^-5,0.6+(0.10*0.6))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off() 


######################################
### Statistics interspecific level ###
######################################

### Half saturation ###

# Test for normality
shapiro.test(log(Data6$HS+1))
qqnorm(log(Data6$HS+1))
qqline(log(Data6$HS+1))
hist(log(Data6$HS+1))

# Test for homoscedasticity
bartlett.test(log(Data6$HS+1), Data6$Species)

# Comparisons of means
kruskal.test(log(Data6$HS+1), Data6$Species)

### Maximum growth rate ###

# Test for normality
shapiro.test(log(Data6$GRM+1))
qqnorm(log(Data6$GRM+1))
qqline(log(Data6$GRM+1))
hist(log(Data6$GRM+1))

# Test for homoscedasticity
bartlett.test(log(Data6$GRM+1), Data6$Species)

# Comparisons of means
kruskal.test(log(Data6$GRM+1), Data6$Species)


######################################
### Statistics intraspecific level ###
######################################

# Subset the dataset
DataB=subset(Data7, Species=="Chlamydomonas bilatus")
DataN=subset(Data7, Species=="Chlamydomonas noctigama")

### Half saturation ###

# Test for normality
shapiro.test(log(DataB$HS+1))
qqnorm(log(DataB$HS+1))
qqline(log(DataB$HS+1))
hist(log(DataB$HS+1))

shapiro.test(log(DataN$HS+1))
qqnorm(log(DataN$HS+1))
qqline(log(DataN$HS+1))
hist(log(DataN$HS+1))

# Test for homoscedasticity
bartlett.test(log(DataB$HS+1), DataB$Iso)
bartlett.test(log(DataN$HS+1), DataN$Iso)

# Comparisons of means
kruskal.test(log(DataB$HS+1), DataB$Iso)
kruskal.test(log(DataN$HS+1), DataN$Iso)

### Maximum growth rate ###

# Test for normality
shapiro.test(log(DataB$GRM+1))
qqnorm(log(DataB$GRM+1))
qqline(log(DataB$GRM+1))
hist(log(DataB$GRM+1))

shapiro.test(log(DataN$GRM+1))
qqnorm(log(DataN$GRM+1))
qqline(log(DataN$GRM+1))
hist(log(DataN$GRM+1))

# Test for homoscedasticity
bartlett.test(log(DataB$GRM+1), DataB$Iso)
bartlett.test(log(DataN$GRM+1), DataN$Iso)

# Comparisons of means
kruskal.test(log(DataB$GRM+1), DataB$Iso)
kruskal.test(log(DataN$GRM+1), DataN$Iso)


####################################################
### Interspecific vs intraspecific contributions ###
####################################################

# Mixed effect regressions
Model1=lmer(HS~1 +(1|Iso), REML=T, data=Data)
Model2=lmer(HS~Species +(1|Iso), REML=T, data=Data)
anova(Model1,Model2)

# Mixed effect regressions
Model1=lmer(GRM~1 +(1|Iso), REML=T, data=Data)
Model2=lmer(GRM~Species +(1|Iso), REML=T, data=Data)
anova(Model1,Model2)
