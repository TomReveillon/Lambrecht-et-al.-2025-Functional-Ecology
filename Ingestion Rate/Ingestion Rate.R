setwd("~/Activité Professionnelle/LIMNO 2019-2023/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Ingestion Rate")

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

######################################
######################################
##### INGESTION RATE OF ISOLATES #####
######################################
######################################

# Import the dataset
Data=read.table("Data_IR.txt", h=T, dec=",")
Data$Species=gsub("\\_","\\ ", Data$Species)
summary(Data)
names(Data)

# Specify the variables as numeric numbers
Data[,c(5:8,10:11)] %<>% mutate_if(is.character,as.numeric)
Data[,c(5:8,10:11)] %<>% mutate_if(is.integer,as.numeric)
Data[,c(5:8,10:11)] %<>% mutate_if(is.factor,as.numeric)

# Sort the dataset by isolate
Data$Iso=factor(Data$Iso, levels=unique(Data$Iso))
Data$Species=factor(Data$Species, levels=unique(Data$Species))
Data=Data[order(Data$Iso,Data$Trial),]

# Calculate the densities
Data$IDens=round(Data$IXMI*Data$Volu*Data$Site*Data$Dilu*Data$Cove,0)
Data$FDens=round(Data$IXMF*Data$Volu*Data$Site*Data$Dilu*Data$Cove,0)

### Ingestion rate ###

# Calculate the ingestion rates
Data$Inges=(Data$IDens-Data$FDens)/(8*60*60)

# Correct by rotifer density
Data$Inges=Data$Inges/5

# Correct by volume
Data$Inges=Data$Inges/5

# Calculate the clearance rates
Data$Clear=(Data$Inges/Data$IDens)*8*60*60

# Remove the negative values
Data$Inges[Data$Inges < 0]=0
Data$Clear[Data$Clear < 0]=0

# Calculate the maximum ingestion and clearance rates
Data=setDT(Data)[, .SD[which.max(Inges)], by=list(Iso,Species,Trial)]

# Calculate the mean ingestion and clearance rates
Data2=setDT(Data)[, .(IR=round(mean(Inges),4), IRLSD=round(mean(Inges)-sd(Inges),4), IRUSD=round(mean(Inges)+sd(Inges),4), CR=round(mean(Clear),4), CRLSD=round(mean(Clear)-sd(Clear),4), CRUSD=round(mean(Clear)+sd(Clear),4)), by=list(Iso,Species)]
Data3=setDT(Data)[, .(IR=round(mean(Inges),4), CR=round(mean(Clear),4)), by=list(Iso,Species,Trial)]
Data2=as.data.frame(Data2)
Data3=as.data.frame(Data3)

# Export the datasets
Data2$IR[Data2$IR < 0]=0; Data2$IRLSD[Data2$IRLSD < 0]=0; Data2$IRUSD[Data2$IRUSD < 0]=0
Data3$IR[Data3$IR < 0]=0
Data2$CR[Data2$CR < 0]=0; Data2$CRLSD[Data2$CRLSD < 0]=0; Data2$CRUSD[Data2$CRUSD < 0]=0
Data3$CR[Data3$CR < 0]=0
Data3=Data3[order(Data3$Iso,Data3$Trial),]
write.table(Data2[,c(1:8)], file="Data_IR_Inter.txt", sep="\t", row.names=F)
write.table(Data3[,c(1:5)], file="Data_IR_Intra.txt", sep="\t", row.names=F)

### Error bars ###

# Calculate the error bars for ingestion rate
Data4=as.data.frame(setDT(na.omit(Data2))[, .(IRSpecies=round(mean(IR),4), IRSpeciesL=round(mean(IR)-sd(IR),4), IRSpeciesU=round(mean(IR)+sd(IR),4)), by=list(Species)])
Data2$IRSpecies=c(rep(Data4$IRSpecies[1],24), rep(Data4$IRSpecies[2],25))
Data2$IRSpeciesL=c(rep(Data4$IRSpeciesL[1],24), rep(Data4$IRSpeciesL[2],25))
Data2$IRSpeciesU=c(rep(Data4$IRSpeciesU[1],24), rep(Data4$IRSpeciesU[2],25))

# Calculate the error bars for clearance rate
Data5=as.data.frame(setDT(na.omit(Data2))[, .(CRSpecies=round(mean(CR),4), CRSpeciesL=round(mean(CR)-sd(CR),4), CRSpeciesU=round(mean(CR)+sd(CR),4)), by=list(Species)])
Data2$CRSpecies=c(rep(Data5$CRSpecies[1],24), rep(Data5$CRSpecies[2],25))
Data2$CRSpeciesL=c(rep(Data5$CRSpeciesL[1],24), rep(Data5$CRSpeciesL[2],25))
Data2$CRSpeciesU=c(rep(Data5$CRSpeciesU[1],24), rep(Data5$CRSpeciesU[2],25))

# Export the dataset
Data2[,c(3:14)]=round(Data2[,c(3:14)],4)
write.table(Data2[,c(1:14)], file="Data_IR_TO.txt", sep="\t", row.names=F)


######################################
### Plot predicted ingestion rates ###
######################################

tiff('Ingestion Rates.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data2, aes(Species, IR, group=Species)) +
  geom_pointrange(aes(Species, ymin=IRLSD, ymax=IRUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=0.45+(0.020*0.45), yend=0.45+(0.020*0.45), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=0.45+(0.020*0.45), yend=0.45+(0.020*0.45), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=0.45+(0.100*0.45), yend=0.45+(0.100*0.45), color="black", size=0.8) +
  geom_text(x=1.00, y=0.45+(0.035*0.45), label="***", color="black", size=5) +
  geom_text(x=2.00, y=0.45+(0.045*0.45), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=0.45+(0.115*0.45), label="**", color="black", size=5) +
  ylab(expression('Maximum ingestion rate'~'('*cells~sec^-1~ind^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.2f", x), breaks=seq(0,0.45,by=0.15), limits=c(-10^-5,0.45+(0.10*0.45))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off() 

tiff('Clearance Rates.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data2, aes(Species, CR*10^2, group=Species)) +
  geom_pointrange(aes(Species, ymin=CRLSD*10^2, ymax=CRUSD*100, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=2.1+(0.020*2.1), yend=2.1+(0.020*2.1), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=2.1+(0.020*2.1), yend=2.1+(0.020*2.1), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=2.1+(0.100*2.1), yend=2.1+(0.100*2.1), color="black", size=0.8) +
  geom_text(x=1.00, y=2.1+(0.035*2.1), label="***", color="black", size=5) +
  geom_text(x=2.00, y=2.1+(0.045*2.1), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=2.1+(0.125*2.1), label="NS", color="black", size=5) +
  ylab(expression('Maximum clearance rate'~'('*10^-2~mL~day^-1~ind^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,2.1,by=0.7), limits=c(-10^-5,2.1+(0.10*2.1))) +
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

### Ingestion rate ###

# Test for normality
shapiro.test(log(Data2$IR+1))
qqnorm(log(Data2$IR+1))
qqline(log(Data2$IR+1))
hist(log(Data2$IR+1))

# Test for homoscedasticity
bartlett.test(log(Data2$IR+1), Data2$Species)

# Comparisons of means
kruskal.test(log(Data2$IR+1), Data2$Species)

### Clearance rate ###

# Test for normality
shapiro.test(log(Data2$CR+1))
qqnorm(log(Data2$CR+1))
qqline(log(Data2$CR+1))
hist(log(Data2$CR+1))

# Test for homoscedasticity
bartlett.test(log(Data2$CR+1), Data2$Species)

# Comparisons of means
kruskal.test(log(Data2$CR+1), Data2$Species)


######################################
### Statistics intraspecific level ###
######################################

# Subset the dataset
DataB=subset(Data3, Species=="Chlamydomonas bilatus")
DataN=subset(Data3, Species=="Chlamydomonas noctigama")

### Ingestion rate ###

# Test for normality
shapiro.test(log(DataB$IR+1))
qqnorm(log(DataB$IR+1))
qqline(log(DataB$IR+1))
hist(log(DataB$IR+1))

shapiro.test(log(DataN$IR+1))
qqnorm(log(DataN$IR+1))
qqline(log(DataN$IR+1))
hist(log(DataN$IR+1))

# Test for homoscedasticity
bartlett.test(log(DataB$IR+1), DataB$Iso)
bartlett.test(log(DataN$IR+1), DataN$Iso)

# Comparisons of means
kruskal.test(log(DataB$IR+1), DataB$Iso)
kruskal.test(log(DataN$IR+1), DataN$Iso)

### Clearance rate ###

# Test for normality
shapiro.test(log(DataB$CR+1))
qqnorm(log(DataB$CR+1))
qqline(log(DataB$CR+1))
hist(log(DataB$CR+1))

shapiro.test(log(DataN$CR+1))
qqnorm(log(DataN$CR+1))
qqline(log(DataN$CR+1))
hist(log(DataN$CR+1))

# Test for homoscedasticity
bartlett.test(log(DataB$CR+1), DataB$Iso)
bartlett.test(log(DataN$CR+1), DataN$Iso)

# Comparisons of means
kruskal.test(log(DataB$CR+1), DataB$Iso)
kruskal.test(log(DataN$CR+1), DataN$Iso)


####################################################
### Interspecific vs intraspecific contributions ###
####################################################

# Mixed effect regressions
Model1=lmer(IR~1 +(1|Iso), REML=T, data=Data)
Model2=lmer(IR~Species +(1|Iso), REML=T, data=Data)
anova(Model1,Model2)

# Mixed effect regressions
Model1=lmer(CR~1 +(1|Iso), REML=T, data=Data)
Model2=lmer(CR~Species +(1|Iso), REML=T, data=Data)
anova(Model1,Model2)
