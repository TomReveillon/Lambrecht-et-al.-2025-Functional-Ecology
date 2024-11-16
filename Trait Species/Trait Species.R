setwd("~/Activité Professionnelle/LIMNO 2019-2023/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Trait Species")

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

###################################################
###################################################
##### TRAITS FITNESS RELATIONSHIP OF ISOLATES #####
###################################################
###################################################

# Import the dataset for ingestion rate of isolates
Data=read.table("Data_IR.txt", h=T, dec=",")
Data$Species=gsub("\\_","\\ ", Data$Species)
summary(Data)
names(Data)

# Import the dataset for traits of isolates
Data2=read.table("Data_Traits.txt", h=T, dec=".")
summary(Data2)
names(Data2)

# Import the dataset for traits of species
Data3=read.table("Data_Species.txt", h=T, dec=".")
summary(Data3)
names(Data3)

# Specify the variables as numeric or factor
Data[,c(5:8,10:11)] %<>% mutate_if(is.character,as.numeric)
Data2[,c(3:29)] %<>% mutate_if(is.character,as.numeric)
Data3[,c(2:13)] %<>% mutate_if(is.character,as.numeric)

# Sort the datasets by species
Data$Species=factor(Data$Species, levels=unique(Data$Species))
Data2$Species=factor(Data2$Species, levels=unique(Data2$Species))
Data3$Species=factor(Data3$Species, levels=unique(Data3$Species))

# Calculate the densities
Data$IDens=round(Data$IXMI*Data$Volu*Data$Site*Data$Dilu*Data$Cove,0)
Data$FDens=round(Data$IXMF*Data$Volu*Data$Site*Data$Dilu*Data$Cove,0)

# Calculate the density differences
Data$DDens=(Data$IDens-Data$FDens)

# Calculate the edibilities
Data$Edib=round(Data$DDens/Data$IDens,4)

# Remove the negative values
Data$Edib[Data$Edib < 0]=0

# Calculate the mean edibilities
Data4=setDT(Data)[, .(ED=round(mean(Edib)/0.8,4), EDLSD=round((mean(Edib)-sd(Edib))/0.8,4), EDUSD=round((mean(Edib)+sd(Edib))/0.8,4)), by=list(Iso,Species)]
Data4=as.data.frame(Data4)

# Combine the datasets
Data2=cbind(Data2[,c(1:29)],Data4[,c(3:5)])
Data2[Data2 < 0]=0


###################################
### Edibility-size relationship ###
###################################

# Include data for strains
ED=c(0.0000,0.0687,0.1659,0.1462,0.5303,0.3056,0.3336,1.0000)
CA=c(1000,323,178,171,108,121,136,0)

# Polynomial regression
ModED=lm(ED~poly(CA,3))
summary(ModED)
anova(ModED)

# Filter edible or non-edible species
Data5=subset(Data3, CA < 1000)
Data6=subset(Data3, CA > 1000)
Data6[,c(2:4)]=0

# Predict edibilities for species
PredED=predict.lm(ModED, data.frame(CA=Data5$CA), se.fit=T)

# Create a dataset
Data5=data.frame(Species=Data5$Species, ED=PredED[[1]], EDLSD=PredED[[1]]-PredED[[2]], EDUSD=PredED[[1]]+PredED[[2]])
Data5=merge(Data5, Data6[,c(1:4)], all=T)
Data5=Data5[order(Data5$Species),]

# Combine the datasets
Data3=cbind(Data5[,c(1:4)],Data3[,c(5:13)])
Data3$ED[Data3$ED < 0]=0
Data3$EDLSD[Data3$EDLSD < 0]=0
Data3$EDUSD[Data3$EDUSD < 0]=0

# Define spaces for isolates
Data7=setDT(Data2)[, .(MinED=min(ED), MaxED=max(ED), MinGR=min(GR), MaxGR=max(GR), MinHS=min(HS), MaxHS=max(HS), MinCA=min(CA), MaxCA=max(CA)), by=list(Species)]
Data7=as.data.frame(Data7)

# Define location for species
Data8=subset(Data3, Species=="Chlamydomonas")[,c(1:2,5,8,11)]


####################################
### Calculating trade-off curves ###
####################################

### Edibility (defense) vs growth rate (competitiveness) ###

# Fitting linear models
ModLi=lm(ED~GR, data=Data3)

# Fitting polynomial models
ModPo=nls(ED~min(ED) - a*(GR-min(GR))^2.0, start=c(a=0.1), data=Data3)

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))

# Calculate the predicted values
Data3$TO.ED.GR=round(predict(ModPo),4)

### Edibility (defense) vs half saturation (competitiveness) ###

# Fitting linear models
ModLi=lm(ED~HS, data=Data3)

# Fitting polynomial models
ModPo=nls(ED~max(ED) - a*(HS+min(HS))^0.2, start=c(a=0.1), data=Data3)

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))

# Calculate the predicted values
Data3$TO.ED.HS=round(predict(ModPo),4)

### Particle area (defense) vs growth rate (competitiveness) ###

# Fitting linear models
ModLi=lm(CA~GR, data=Data3)

# Fitting polynomial models
ModPo=nls(CA~max(CA) + a*(GR+min(GR))^2.0, start=c(a=0.1), data=Data3)

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))

# Calculate the predicted values
Data3$TO.CA.GR=round(predict(ModPo),4)

### Particle area (defense) vs half saturation (competitiveness) ###

# Fitting linear models
ModLi=lm(CA~HS, data=Data3)

# Fitting polynomial models
ModPo=nls(CA~min(CA) + a*(HS+min(HS))^0.5, start=c(a=0.1), data=Data3)

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))

# Calculate the predicted values
Data3$TO.CA.HS=round(predict(ModPo),4)


#################################################
### Plot defense-competitiveness trait spaces ###
#################################################

Plot1=ggplot(Data3) + coord_cartesian(clip="off") +
  geom_smooth(aes(GR, TO.ED.GR), color="black", linetype="solid", size=1.0, se=F) +
  geom_rect(data=Data7[1,], aes(xmin=MinGR, xmax=MaxGR, ymin=MinED, ymax=MaxED), fill=NA, color="steelblue2", linetype="solid", size=1.5) +
  geom_rect(data=Data7[2,], aes(xmin=MinGR, xmax=MaxGR, ymin=MinED, ymax=MaxED), fill=NA, color="tomato3", linetype="solid", size=1.5) +
  geom_point(aes(GR, ED, color=Species), fill="white", size=2, pch=16) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=2.4, ymin=1.0+(1.0-0.0)*0.010, ymax=1.0+(1.0-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=2.4+(2.4+0.0)*0.010, xmax=2.4+(2.4+0.0)*0.010, ymin=0.0, ymax=1.0) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=2.4+(2.4+0.0)*0.040, xmax=2.4+(2.4+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=1.0+(1.0-0.0)*0.040, ymax=1.0+(1.0-0.0)*0.040) + 
  ylab(expression('Edibility')) + xlab(expression('Maximum growth rate'~'('*day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,1.0,by=0.25),fmt="%.1f"), breaks=seq(0,1.0,by=0.25), limits=c(0,1.0)) +
  scale_x_continuous(labels=sprintf(seq(0,2.4,by=0.6),fmt="%.1f"), breaks=seq(0,2.4,by=0.6), limits=c(0,2.4)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(legend.text=element_text(face="plain", colour="black", size=8), legend.key=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot2=ggplot(Data3) + coord_cartesian(clip="off") +
  geom_smooth(aes(HS, TO.ED.HS, color=Species), color="black", linetype="solid", size=1.0, se=F) +
  geom_rect(data=Data7[1,], aes(xmin=MinHS, xmax=MaxHS, ymin=MinED, ymax=MaxED), fill=NA, color="steelblue2", linetype="solid", size=1.5) +
  geom_rect(data=Data7[2,], aes(xmin=MinHS, xmax=MaxHS, ymin=MinED, ymax=MaxED), fill=NA, color="tomato3", linetype="solid", size=1.5) +
  geom_point(aes(HS, ED, color=Species), fill="white", size=2, pch=16) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=2.8, ymin=1.0+(1.0-0.0)*0.010, ymax=1.0+(1.0-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=2.8+(2.8+0.0)*0.010, xmax=2.8+(2.8+0.0)*0.010, ymin=0.0, ymax=1.0) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=2.8+(2.8+0.0)*0.040, xmax=2.8+(2.8+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=1.0+(1.0-0.0)*0.040, ymax=1.0+(1.0-0.0)*0.040) + 
  ylab(expression('Edibility')) + xlab(expression('Half-saturation constant'~'('*µmol~PO[4]^{"-"}~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,1.0,by=0.25),fmt="%.1f"), breaks=seq(0,1.0,by=0.25), limits=c(0,1.0)) +
  scale_x_continuous(labels=sprintf(seq(0,2.8,by=0.7),fmt="%.1f"), breaks=seq(0,2.8,by=0.7), limits=c(0,2.8)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(legend.text=element_text(face="plain", colour="black", size=8), legend.key=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot3=ggplot(Data3) + coord_cartesian(clip="off") +
  geom_smooth(aes(GR, log(TO.CA.GR+1)), color="black", linetype="solid", size=1.0, se=F) +
  geom_rect(data=Data7[1,], aes(xmin=MinGR, xmax=MaxGR, ymin=log(MinCA+1), ymax=log(MaxCA+1)), fill=NA, color="steelblue2", linetype="solid", size=1.5) +
  geom_rect(data=Data7[2,], aes(xmin=MinGR, xmax=MaxGR, ymin=log(MinCA+1), ymax=log(MaxCA+1)), fill=NA, color="tomato3", linetype="solid", size=1.5) +
  geom_point(aes(GR, log(CA+1), color=Species), fill="white", size=2, pch=16) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=2.4, ymin=9.5+(9.5-3.5)*0.010, ymax=9.5+(9.5-3.5)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=2.4+(2.4+0.0)*0.010, xmax=2.4+(2.4+0.0)*0.010, ymin=3.5, ymax=9.5) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=2.4+(2.4+0.0)*0.040, xmax=2.4+(2.4+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=9.5+(9.5-3.5)*0.040, ymax=9.5+(9.5-3.5)*0.040) + 
  ylab(expression('Particle area'~'('*ln~µm^2*')')) + xlab(expression('Maximum growth rate'~'('*day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(3.5,9.5,by=1.5),fmt="%.1f"), breaks=seq(3.5,9.5,by=1.5), limits=c(3.5,9.5)) +
  scale_x_continuous(labels=sprintf(seq(0,2.4,by=0.6),fmt="%.1f"), breaks=seq(0,2.4,by=0.6), limits=c(0,2.4)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(legend.text=element_text(face="plain", colour="black", size=8), legend.key=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,2.5),"pt")) +
  theme(legend.position="none")

Plot4=ggplot(Data3) + coord_cartesian(clip="off") +
  geom_smooth(aes(HS, log(TO.CA.HS+1)), color="black", linetype="solid", size=1.0, se=F) +
  geom_rect(data=Data7[1,], aes(xmin=MinHS, xmax=MaxHS, ymin=log(MinCA+1), ymax=log(MaxCA+1)), fill=NA, color="steelblue2", linetype="solid", size=1.5) +
  geom_rect(data=Data7[2,], aes(xmin=MinHS, xmax=MaxHS, ymin=log(MinCA+1), ymax=log(MaxCA+1)), fill=NA, color="tomato3", linetype="solid", size=1.5) +
  geom_point(aes(HS, log(CA+1), color=Species), fill="white", size=2, pch=16) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=2.8, ymin=9.5+(9.5-3.5)*0.010, ymax=9.5+(9.5-3.5)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=2.8+(2.8+0.0)*0.010, xmax=2.8+(2.8+0.0)*0.010, ymin=3.5, ymax=9.5) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=2.8+(2.8+0.0)*0.040, xmax=2.8+(2.8+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=9.5+(9.5-3.5)*0.040, ymax=9.5+(9.5-3.5)*0.040) + 
  ylab(expression('Particle area'~'('*ln~µm^2*')')) + xlab(expression('Half-saturation constant'~'('*µmol~PO[4]^{"-"}~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(3.5,9.5,by=1.5),fmt="%.1f"), breaks=seq(3.5,9.5,by=1.5), limits=c(3.5,9.5)) +
  scale_x_continuous(labels=sprintf(seq(0,2.8,by=0.7),fmt="%.1f"), breaks=seq(0,2.8,by=0.7), limits=c(0,2.8)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(legend.text=element_text(face="plain", colour="black", size=8), legend.key=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,2.5),"pt")) +
  theme(legend.position="none")

# Panel plot of trait spaces
tiff('Trait Species.tiff', units="in", width=12, height=12, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4)
grid.arrange(grobs=Panel, ncol=2, nrow=2)
dev.off()
