setwd("~/Activité Professionnelle/LIMNO 2019-2023/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Growth Rate")

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

###################################
###################################
##### GROWTH RATE OF ISOLATES #####
###################################
###################################

# Import the dataset
Data=read.table("Data_GR.txt", h=T, dec=",", row.names=NULL)
Data$Species=gsub("\\_","\\ ", Data$Species)
summary(Data)
names(Data)

# Specify the variables as numeric or factor
Data[,c(5:8,11)] %<>% mutate_if(is.character,as.numeric)
Data[,c(5:8,11)] %<>% mutate_if(is.integer,as.numeric)
Data[,c(5:8,11)] %<>% mutate_if(is.factor,as.numeric)

# Sort the dataset by isolate and nutrient
Data$Iso=factor(Data$Iso, levels=unique(Data$Iso))
Data$Species=factor(Data$Species, levels=unique(Data$Species))
Data=Data[order(Data$Iso,Data$Trial),]

# Calculate the densities
Data$Dens=Data$IXM*Data$Volu*Data$Site*Data$Dilu*Data$Cove

# Calculate the mean densities
Data2=setDT(Data)[, .(MeanDens=round(mean(Dens),0)), by=list(Iso,Species,Day)]
Data2=as.data.frame(Data2)


##########################
### Growth rate models ###
##########################

### Linear model ###

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
  
  Out=list(Model=Model, Parameters=Parameters, Rates=Rates)
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
  
  Out=list(Model=Model, Parameters=Parameters, Rates=Rates)
  return(Out)
}


##################################
### Fitting growth rate models ###
##################################

# Separate the datasets
DataLi=subset(Data, Iso=="I24"|Iso=="I37")
DataLo=subset(Data, !c(Iso=="I24"|Iso=="I37"))

### Linear model ###

# Split the dataset by isolate
SplitDataLi=split(DataLi, list(DataLi$Iso))
SplitDataLi=SplitDataLi[sapply(SplitDataLi, function(x) dim(x)[1]) > 0]

# Extract the combinations of names
Names=unique(DataLi[,c("Iso","Species","Day")])
Names=Names[order(Names$Iso, Names$Species),]
Iso=rep(Names$Iso,each=5)
Species=rep(Names$Species,each=5)

# Fitting the linear model
ModLi=function(x) {
  FitLi=nls(log(Dens+1) ~ SSline(Day, a, b), data=x)
  OutLi=FuncLi(FitLi, x$Day)
} 
OutLi=lapply(SplitDataLi, ModLi)

# Create a dataset
RateLi=bind_rows(lapply(OutLi, function (x) x[c("Rates")]))
RateLi=round(cbind(do.call("rbind",RateLi)),4)
RateLi=cbind(Iso,Species,Phos="Low",Model="Linear",RateLi)
rownames(RateLi)=c()

### Logistic model ###

# Split the dataset by isolate
SplitDataLo=split(DataLo, list(DataLo$Iso))
SplitDataLo=SplitDataLo[sapply(SplitDataLo, function(x) dim(x)[1]) > 0]

# Extract the combinations of names
Names=unique(DataLo[,c("Iso","Species","Day")])
Names=Names[order(Names$Iso, Names$Species),]
Iso=rep(Names$Iso,each=5)
Species=rep(Names$Species,each=5)

# Logistic model
ModLo=function(x) {
  FitLo=nls(log(Dens+1) ~ SSlogis(Day, b, K, r), data=x)
  OutLo=FuncLo(FitLo, x$Day)
}
OutLo=lapply(SplitDataLo, ModLo)

# Create a dataset
RateLo=bind_rows(lapply(OutLo, function (x) x[c("Rates")]))
RateLo=round(cbind(do.call("rbind",RateLo)),4)
RateLo=cbind(Iso,Species,Phos="Low",Model="Logistic",RateLo)
rownames(RateLo)=c()


####################################
### Plot predicted growth curves ###
####################################

# Combine the datasets
Data3=rbind(RateLi,RateLo)
Data3[,c(6:8)]=exp(Data3[,c(6:8)])
Data3[,c(6:8)]=round(Data3[,c(6:8)],4)
Data3=Data3[order(Data3$Iso,Data3$Day),]
Data3=Data3[!duplicated(Data3),]

tiff('Growth Curves.tiff', units="in", width=15, height=15, res=1000)
ggplot(Data3, aes(Day, log(DensP+1), group=Species)) + 
  geom_line(aes(color=Species), linetype="solid", size=1) +
  geom_point(data=Data, aes(Day, log(Dens+1), color=Species), alpha=0.4, size=1.5, pch=16, position=position_jitter(0.2)) +
  ylab(expression('Density'~'('*ln~cells~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=16)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=16)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(9,15,by=2), limits=c(9,15)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(2,20,by=6), limits=c(2,20)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas bilatus"="tomato3","Chlamydomonas noctigama"="steelblue2")) +
  theme(strip.text.x=element_text(face="bold", color="black", size=16)) +
  theme(strip.background=element_blank()) +
  facet_wrap(~Iso, scales="free", ncol=7, nrow=7) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()


################################
### Calculating growth rates ###
################################

### Carrying capacity ###

# Calculate carrying capacity
CCIntra=setDT(Data)[, .(CC=round(mean(Dens),0)), by=list(Iso,Species,Trial)]
CCInter=setDT(CCIntra)[, .(CC=round(mean(CC),0), CCLSD=round(mean(CC)-sd(CC),0), CCUSD=round(mean(CC)+sd(CC),0)), by=list(Iso,Species)]
CCIntra=as.data.frame(CCIntra)
CCInter=as.data.frame(CCInter)

### Growth rate ###

# Subset the dataset
Data4=subset(Data, Day < 8)
Data=as.data.frame(Data)

# Split the dataset
SplitData4=split(Data4, list(Data4$Iso,Data4$Trial))

# Extract combinations of names
Names=unique(Data[,c("Iso","Species","Trial")])
Names=Names[order(Names$Iso, Names$Species),]
Iso=Names$Iso; Species=Names$Species; Trial=Names$Trial

# Linear model
ModGR=function(x) {coef(summary(lm(log(Dens+1)~Day, data=x)))[c(2,4)]}
OutGR=lapply(SplitData4, ModGR)
OutGR=round(cbind(do.call("rbind",OutGR)),4)
OutGR=data.frame(Iso=Iso,Species=Species,Trial=Trial,RGR=OutGR[,1],RGRSD=OutGR[,2])
rownames(OutGR)=c()

# Create a dataset
Data5=data.frame(OutGR[,c(1:3)],RGR=OutGR[,4],RGRLSD=OutGR[,4]-OutGR[,5],RGRUSD=OutGR[,4]+OutGR[,5])
Data5=cbind(Data5[,c(1:6)],AGR=exp(Data5[,4]),AGRLSD=exp(Data5[,5]),AGRUSD=exp(Data5[,6]))
Data5[,c(4:9)][Data5[,c(4:9)]<0]=0

# Calculate maximum growth rate
GRIntra=setDT(Data5)[, .(GR=round(max(AGR),4)), by=list(Iso,Species,Trial)]
GRInter=setDT(GRIntra)[, .(GR=round(mean(GR),4), GRLSD=round(mean(GR)-sd(GR),4), GRUSD=round(mean(GR)+sd(GR),4)), by=list(Iso,Species)]
GRIntra=as.data.frame(GRIntra)
GRInter=as.data.frame(GRInter)

# Create the datasets
Data5=data.frame(Iso=GRInter[,1], Species=GRInter[,2], GR=GRInter[,3], GRLSD=GRInter[,4], GRUSD=GRInter[,5], CC=CCInter[,3], CCLSD=CCInter[,4], CCUSD=CCInter[,5])
Data6=data.frame(Iso=GRIntra[,1], Species=GRIntra[,2], Trial=GRIntra[,3], GR=GRIntra[,4], CC=CCIntra[,4])

# Export the datasets
Data5$GR[Data5$GR < 0]=0; Data5$GRLSD[Data5$GRLSD < 0]=0; Data5$GRUSD[Data5$GRUSD < 0]=0
Data6$GR[Data6$GR < 0]=0
Data5$CC[Data5$CC < 0]=0; Data5$CCLSD[Data5$CCLSD < 0]=0; Data5$CCUSD[Data5$CCUSD < 0]=0
Data6$CC[Data6$CC < 0]=0
Data6=Data6[order(Data6$Iso,Data6$Trial),]
write.table(Data5[,c(1:8)], file="Data_GR_Inter.txt", sep="\t", row.names=F)
write.table(Data6[,c(1:5)], file="Data_GR_Intra.txt", sep="\t", row.names=F)

### Error bars ###

# Calculate the error bars for growth rate
Data7=as.data.frame(setDT(na.omit(Data5))[, .(GRSpecies=round(mean(GR),4), GRSpeciesL=round(mean(GR)-sd(GR),4), GRSpeciesU=round(mean(GR)+sd(GR),4)), by=list(Species)])
Data5$GRSpecies=c(rep(Data7$GRSpecies[1],24), rep(Data7$GRSpecies[2],25))
Data5$GRSpeciesL=c(rep(Data7$GRSpeciesL[1],24), rep(Data7$GRSpeciesL[2],25))
Data5$GRSpeciesU=c(rep(Data7$GRSpeciesU[1],24), rep(Data7$GRSpeciesU[2],25))

# Calculate the error bars for carrying capacity
Data8=as.data.frame(setDT(na.omit(Data5))[, .(CCSpecies=round(mean(CC),0), CCSpeciesL=round(mean(CC)-sd(CC),0), CCSpeciesU=round(mean(CC)+sd(CC),0)), by=list(Species)])
Data5$CCSpecies=c(rep(Data8$CCSpecies[1],24), rep(Data8$CCSpecies[2],25))
Data5$CCSpeciesL=c(rep(Data8$CCSpeciesL[1],24), rep(Data8$CCSpeciesL[2],25))
Data5$CCSpeciesU=c(rep(Data8$CCSpeciesU[1],24), rep(Data8$CCSpeciesU[2],25))

# Export the dataset
Data5[,c(3:14)]=round(Data5[,c(3:14)],4)
write.table(Data5[,c(1:14)], file="Data_GR_TO.txt", sep="\t", row.names=F)


###################################
### Plot predicted growth rates ###
###################################

tiff('Maximum Growth Rates.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data5, aes(Species, GR, group=Species)) +
  geom_pointrange(aes(Species, ymin=GRLSD, ymax=GRUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=1.9+(0.020*0.9), yend=1.9+(0.020*0.9), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=1.9+(0.020*0.9), yend=1.9+(0.020*0.9), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=1.9+(0.100*0.9), yend=1.9+(0.100*0.9), color="black", size=0.8) +
  geom_text(x=1.00, y=1.9+(0.035*0.9), label="***", color="black", size=5) +
  geom_text(x=2.00, y=1.9+(0.035*0.9), label="***", color="black", size=5) +
  geom_text(x=1.50, y=1.9+(0.125*0.9), label="NS", color="black", size=5) +
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(1.0,1.9,by=0.3), limits=c(1.0,1.9+(0.10*0.9))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off() 

tiff('Carrying Capacities.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data5, aes(Species, CC/10^5, group=Species)) +
  geom_errorbar(aes(Species, CC/10^5, ymin=CCSpeciesL/10^5, ymax=CCSpeciesU/10^5, color=Species), linetype="solid", alpha=0.3, width=0.2, linewidth=1, position=position_dodge(0.5)) +
  geom_errorbar(aes(Species, CC/10^5, ymin=CCSpecies/10^5, ymax=CCSpecies/10^5, color=Species), linewidth=1, width=0.2, linetype="solid", position=position_dodge(w=0.5)) +
  geom_pointrange(aes(Species, ymin=CCLSD/10^5, ymax=CCUSD/10^5, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=5.4+(0.020*5.4), yend=5.4+(0.020*5.4), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=5.4+(0.020*5.4), yend=5.4+(0.020*5.4), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=5.4+(0.100*5.4), yend=5.4+(0.100*5.4), color="black", size=0.8) +
  geom_text(x=1.00, y=5.4+(0.045*5.4), label="NS", color="black", size=5) +
  geom_text(x=2.00, y=5.4+(0.045*5.4), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=5.4+(0.115*5.4), label="***", color="black", size=5) +
  ylab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,5.4,by=1.8), limits=c(-10^-5,5.4+(0.10*5.4))) +
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

### Growth rate ###

# Test for normality
shapiro.test(log(Data5$GR+1))
qqnorm(log(Data5$GR+1))
qqline(log(Data5$GR+1))
hist(log(Data5$GR+1))

# Test for homoscedasticity
bartlett.test(log(Data5$GR+1), Data5$Species)

# Comparisons of means
kruskal.test(log(Data5$GR+1), Data5$Species)

### Carrying capacity ###

# Test for normality
shapiro.test(log(Data5$CC+1))
qqnorm(log(Data5$CC+1))
qqline(log(Data5$CC+1))
hist(log(Data5$CC+1))

# Test for homoscedasticity
bartlett.test(log(Data5$CC+1), Data5$Species)

# Comparisons of means
kruskal.test(log(Data5$CC+1), Data5$Species)


######################################
### Statistics intraspecific level ###
######################################

# Subset the dataset
DataB=subset(Data6, Species=="Chlamydomonas bilatus")
DataN=subset(Data6, Species=="Chlamydomonas noctigama")

### Growth rate ###

# Test for normality
shapiro.test(log(DataB$GR+1))
qqnorm(log(DataB$GR+1))
qqline(log(DataB$GR+1))
hist(log(DataB$GR+1))

shapiro.test(log(DataN$GR+1))
qqnorm(log(DataN$GR+1))
qqline(log(DataN$GR+1))
hist(log(DataN$GR+1))

# Test for homoscedasticity
bartlett.test(log(DataB$GR+1), DataB$Iso)
bartlett.test(log(DataN$GR+1), DataN$Iso)

# Comparisons of means
kruskal.test(log(DataB$GR+1), DataB$Iso)
kruskal.test(log(DataN$GR+1), DataN$Iso)

### Carrying capacity ###

# Test for normality
shapiro.test(log(DataB$CC+1))
qqnorm(log(DataB$CC+1))
qqline(log(DataB$CC+1))
hist(log(DataB$CC+1))

shapiro.test(log(DataN$CC+1))
qqnorm(log(DataN$CC+1))
qqline(log(DataN$CC+1))
hist(log(DataN$CC+1))

# Test for homoscedasticity
bartlett.test(log(DataB$CC+1), DataB$Iso)
bartlett.test(log(DataN$CC+1), DataN$Iso)

# Comparisons of means
kruskal.test(log(DataB$CC+1), DataB$Iso)
kruskal.test(log(DataN$CC+1), DataN$Iso)


####################################################
### Interspecific vs intraspecific contributions ###
####################################################

# Mixed effect regressions
Model1=lmer(GR~1 +(1|Iso), REML=T, data=Data)
Model2=lmer(GR~Species +(1|Iso), REML=T, data=Data)
anova(Model1,Model2)

# Mixed effect regressions
Model1=lmer(CC~1 +(1|Iso), REML=T, data=Data)
Model2=lmer(CC~Species +(1|Iso), REML=T, data=Data)
anova(Model1,Model2)
