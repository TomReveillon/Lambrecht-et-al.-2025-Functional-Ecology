setwd("~/Activité Professionnelle/LIMNO 2019-2023/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Morphology")

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

##################################
##################################
##### MORPHOLOGY OF ISOLATES #####
##################################
##################################

# Import the dataset
Folder="Data Isolates"
Files=list.files(path=Folder, pattern="*.txt", full.names=T)
Data=ldply(Files, read.table, sep="\t", fill=T, header=T, dec=",")
summary(Data)
names(Data)

# Specify the variables as numeric or factor
Data[,c(1:468)] %<>% mutate_if(is.character,as.numeric)
Data[,c(1:468)] %<>% mutate_if(is.integer,as.numeric)
Data[,c(1:468)] %<>% mutate_if(is.factor,as.numeric)

# Create the combinations of names
Data$Code=rep(c("3A_B1","3A_B5","3A_C1","3A_C4","3A_D5","3B_A4","3B_A6","3B_B6","3B_C4","3B_D4","3B_D6","4A_A3","4A_A5","4A_B2","4A_B5","4A_C3","4A_C4","4B_A2","4B_A5","4B_B2","4B_B3","4B_B4","4B_B6","4B_C2","18A_A1","18A_A2","18A_A3","18A_B4","18A_B5","18A_C3","18A_C5","18A_D1","18A_D4","18A_D6","18B_A1","18B_A3","18B_A5","18B_B1","18B_B2","18B_B4","18B_B6","18B_C3","18B_C4","18B_C5","18B_C6","18B_D2","18B_D3","21A_D6","21B_c6"), each=5000)
Data$Iso=rep(c("I2","I3","I4","I5","I7","I8","I9","I10","I11","I13","I14","I15","I16","I17","I18","I19","I20","I23","I24","I25","I26","I27","I28","I29","I31","I32","I33","I34","I35","I36","I37","I38","I40","I41","I42","I43","I44","I45","I46","I47","I48","I49","I50","I51","I52","I53","I54","I57","I58"), each=5000)
Data$Species=c(rep("Chlamydomonas noctigama",5000*24),rep("Chlamydomonas bilatus",5000*25))
Data$Period=c(rep("Oligotrophic",5000*24), rep("Eutrophic",5000*25))
colnames(Data)=gsub("\\_", "\\.", colnames(Data))

# Sort the dataset by isolate
Data$Iso=factor(Data$Iso, levels=unique(Data$Iso))
Data$Species=factor(Data$Species, levels=unique(Data$Species))
Data=Data[order(Data$Iso),]

# Calculate the mean cell areas
Data2=setDT(Data)[, .(CA=round(mean(Area.M05),4), CALSD=round(mean(Area.M05)-sd(Area.M05),4), CAUSD=round(mean(Area.M05)+sd(Area.M05),4), CD=round(mean(Diameter.M05),4), CDLSD=round(mean(Diameter.M05)-sd(Diameter.M05),4), CDUSD=round(mean(Diameter.M05)+sd(Diameter.M05),4)), by=list(Iso,Species)]
Data3=Data[,c("Iso","Species","Area.M05","Diameter.M05"),]; colnames(Data3)[c(3:4)]=c("CA","CD"); Data3[,c(3:4)]=round(Data3[,c(3:4)],4)
Data2=as.data.frame(Data2)
Data3=as.data.frame(Data3)

# Export the datasets
Data2$CA[Data2$CA < 0]=0; Data2$CALSD[Data2$CALSD < 0]=0; Data2$CAUSD[Data2$CAUSD < 0]=0
Data3$CA[Data3$CA < 0]=0
Data2$CD[Data2$CD < 0]=0; Data2$CDLSD[Data2$CDLSD < 0]=0; Data2$CDUSD[Data2$CDUSD < 0]=0
Data3$CD[Data3$CD < 0]=0
write.table(Data2[,c(1:8)], file="Data_MO_Inter.txt", sep="\t", row.names=F)
write.table(Data3[,c(1:4)], file="Data_MO_Intra.txt", sep="\t", row.names=F)

### Error bars ###

# Calculate the error bars for cell area
Data4=as.data.frame(setDT(na.omit(Data2))[, .(CASpecies=round(mean(CA),4), CASpeciesL=round(mean(CA)-sd(CA),4), CASpeciesU=round(mean(CA)+sd(CA),4)), by=list(Species)])
Data2$CASpecies=c(rep(Data4$CASpecies[1],24), rep(Data4$CASpecies[2],25))
Data2$CASpeciesL=c(rep(Data4$CASpeciesL[1],24), rep(Data4$CASpeciesL[2],25))
Data2$CASpeciesU=c(rep(Data4$CASpeciesU[1],24), rep(Data4$CASpeciesU[2],25))

# Calculate the error bars for cell diameter
Data5=as.data.frame(setDT(na.omit(Data2))[, .(CDSpecies=round(mean(CD),4), CDSpeciesL=round(mean(CD)-sd(CD),4), CDSpeciesU=round(mean(CD)+sd(CD),4)), by=list(Species)])
Data2$CDSpecies=c(rep(Data5$CDSpecies[1],24), rep(Data5$CDSpecies[2],25))
Data2$CDSpeciesL=c(rep(Data5$CDSpeciesL[1],24), rep(Data5$CDSpeciesL[2],25))
Data2$CDSpeciesU=c(rep(Data5$CDSpeciesU[1],24), rep(Data5$CDSpeciesU[2],25))

# Export the dataset
Data2[,c(3:14)]=round(Data2[,c(3:14)],4)
write.table(Data2[,c(1:14)], file="Data_MO_TO.txt", sep="\t", row.names=F)


#################################
### Plot predicted cell areas ###
#################################

tiff('Cell Areas.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data2, aes(Species, CA/10^2, group=Species)) +
  geom_pointrange(aes(Species, ymin=CALSD/10^2, ymax=CAUSD/10^2, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=7.5+(0.020*7.5), yend=7.5+(0.020*7.5), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=7.5+(0.020*7.5), yend=7.5+(0.020*7.5), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=7.5+(0.100*7.5), yend=7.5+(0.100*7.5), color="black", size=0.8) +
  geom_text(x=1.00, y=7.5+(0.045*7.5), label="NS", color="black", size=5) +
  geom_text(x=2.00, y=7.5+(0.045*7.5), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=7.5+(0.115*7.5), label="**", color="black", size=5) +
  ylab(expression('Particle area'~'('*10^2~µm^2*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,7.5,by=2.5), limits=c(-10^-5,7.5+(0.10*7.5))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,16.5),"pt")) +
  theme(legend.position="none")
dev.off() 

tiff('Cell Diameters.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data2, aes(Species, CD/10^1, group=Species)) +
  geom_pointrange(aes(Species, ymin=CDLSD/10^1, ymax=CDUSD/10^1, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=2.8+(0.020*2.1), yend=2.8+(0.020*2.1), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=2.8+(0.020*2.1), yend=2.8+(0.020*2.1), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=2.8+(0.100*2.1), yend=2.8+(0.100*2.1), color="black", size=0.8) +
  geom_text(x=1.00, y=2.8+(0.045*2.1), label="NS", color="black", size=5) +
  geom_text(x=2.00, y=2.8+(0.045*2.1), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=2.8+(0.115*2.1), label="**", color="black", size=5) +
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0.7,2.8,by=0.7), limits=c(0.7,2.8+(0.10*2.1))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,16.5),"pt")) +
  theme(legend.position="none")
dev.off() 


######################################
### Statistics interspecific level ###
######################################

# Subset the dataset
Data3=as.data.frame(Data3 %>% group_by(Iso,Species) %>% filter(row_number()==1000))

### Cell area ###

# Test for normality
shapiro.test(log(Data3$CA+1))
qqnorm(log(Data3$CA+1))
qqline(log(Data3$CA+1))
hist(log(Data3$CA+1))

# Test for homoscedasticity
bartlett.test(log(Data3$CA+1), Data3$Species)

# Comparisons of means
kruskal.test(log(Data3$CA+1), Data3$Species)

### Cell diameter ###

# Test for normality
shapiro.test(log(Data3$CD+1))
qqnorm(log(Data3$CD+1))
qqline(log(Data3$CD+1))
hist(log(Data3$CD+1))

# Test for homoscedasticity
bartlett.test(log(Data3$CD+1), Data3$Species)

# Comparisons of means
kruskal.test(log(Data3$CD+1), Data3$Species)


######################################
### Statistics intraspecific level ###
######################################

# Subset the dataset
DataB=subset(Data3, Species=="Chlamydomonas bilatus")
DataN=subset(Data3, Species=="Chlamydomonas noctigama")

### Cell area ###

# Test for normality
shapiro.test(log(DataB$CA+1))
qqnorm(log(DataB$CA+1))
qqline(log(DataB$CA+1))
hist(log(DataB$CA+1))

shapiro.test(log(DataN$CA+1))
qqnorm(log(DataN$CA+1))
qqline(log(DataN$CA+1))
hist(log(DataN$CA+1))

# Test for homoscedasticity
bartlett.test(log(DataB$CA+1), DataB$Iso)
bartlett.test(log(DataN$CA+1), DataN$Iso)

# Comparisons of means
kruskal.test(log(DataB$CA+1), DataB$Iso)
kruskal.test(log(DataN$CA+1), DataN$Iso)

### Cell diameter ###

# Test for normality
shapiro.test(log(DataB$CD+1))
qqnorm(log(DataB$CD+1))
qqline(log(DataB$CD+1))
hist(log(DataB$CD+1))

shapiro.test(log(DataN$CD+1))
qqnorm(log(DataN$CD+1))
qqline(log(DataN$CD+1))
hist(log(DataN$CD+1))

# Test for homoscedasticity
bartlett.test(log(DataB$CD+1), DataB$Iso)
bartlett.test(log(DataN$CD+1), DataN$Iso)

# Comparisons of means
kruskal.test(log(DataB$CD+1), DataB$Iso)
kruskal.test(log(DataN$CD+1), DataN$Iso)


####################################################
### Interspecific vs intraspecific contributions ###
####################################################

# Mixed effect regressions
Model1=lmer(CA~1 +(1|Iso), REML=T, data=Data)
Model2=lmer(CA~Species +(1|Iso), REML=T, data=Data)
anova(Model1,Model2)

# Mixed effect regressions
Model1=lmer(CD~1 +(1|Iso), REML=T, data=Data)
Model2=lmer(CD~Species +(1|Iso), REML=T, data=Data)
anova(Model1,Model2)
