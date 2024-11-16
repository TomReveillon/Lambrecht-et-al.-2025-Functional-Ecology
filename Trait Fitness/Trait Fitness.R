setwd("~/Activité Professionnelle/LIMNO 2019-2023/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Trait Fitness")

rm(list=ls())

library(car)
library(cowplot)
library(data.table)
library(DescTools)
library(deSolve)
library(dplyr)
library(emmeans)
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
library(mgcv)
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

# Import the dataset for traits
DataT=read.table("Data_Traits.txt", h=T, dec=".")
summary(DataT)
names(DataT)

# Import the dataset for fitness
DataF=read.table("Data_Fitness.txt", h=T, dec=".")
summary(DataF)
names(DataF)

# Specify the variables as numeric or factor
DataT[,c(3:29)] %<>% mutate_if(is.character,as.numeric)
DataF[,c(8:9)] %<>% mutate_if(is.character,as.numeric)

# Sort the datasets by isolate
DataT$Iso=factor(DataT$Iso, levels=unique(DataT$Iso))
DataF$Iso=factor(DataF$Iso, levels=unique(DataF$Iso))
DataT$Species=factor(DataT$Species, levels=unique(DataT$Species))
DataF$Species=factor(DataF$Species, levels=unique(DataT$Species))
DataF$Nutri=factor(DataF$Nutri, levels=c("Low","High"))
DataF$Pred=factor(DataF$Pred, levels=c("Without","With"))
DataF=DataF[order(DataF$Iso,DataF$Nutri,DataF$Pred),]

# Calculate the mean areas under the curve
DataF=setDT(DataF)[, .(FitA=mean(GrowA), FitALSD=mean(GrowA)-sd(GrowA), FitAUSD=mean(GrowA)+sd(GrowA)), by=list(Iso,Species,Nutri,Pred)]
DataF=as.data.frame(DataF)

# Extract the combinations of names
Nutri=DataF$Nutri
Pred=DataF$Pred

# Format the datasets
DataF2=dcast(DataF, Iso + Species ~ Nutri + Pred, value.var="FitA")
DataF3=dcast(DataF, Iso + Species ~ Nutri + Pred, value.var="FitALSD")
DataF4=dcast(DataF, Iso + Species ~ Nutri + Pred, value.var="FitAUSD")
colnames(DataF2)[3:6]=c("FitA100","FitA101","FitA1000","FitA1001")
colnames(DataF3)[3:6]=c("FitALSD100","FitALSD101","FitALSD1000","FitALSD1001")
colnames(DataF4)[3:6]=c("FitAUSD100","FitAUSD101","FitAUSD1000","FitAUSD1001")

# Combine the datasets
Data=cbind(DataT[,c(1:29)],DataF2[,c(3:6)],DataF3[,c(3:6)],DataF4[,c(3:6)])


#################################
### Plot trait-fitness spaces ###
#################################

# Melt the dataset
MeltData=melt(Data[,c(1:2,6:8,9:11,12:14,15:17,27:29,30:33)], id.vars=c(colnames(Data[c(1:2,6:8,9:11,12:14,15:17,27:29)]))); colnames(MeltData)[18:19]=c("Treatment","FitA")
MeltDataLSD=melt(Data[,c(1:2,6:8,9:11,12:14,15:17,27:29,34:37)], id.vars=c(colnames(Data[c(1:2,6:8,9:11,12:14,15:17,27:29)]))); colnames(MeltDataLSD)[18:19]=c("Treatment","FitALSD")
MeltDataUSD=melt(Data[,c(1:2,6:8,9:11,12:14,15:17,27:29,38:41)], id.vars=c(colnames(Data[c(1:2,6:8,9:11,12:14,15:17,27:29)]))); colnames(MeltDataUSD)[18:19]=c("Treatment","FitAUSD")

# Include the standard deviations
MeltData$FitALSD=MeltDataLSD$FitALSD
MeltData$FitAUSD=MeltDataUSD$FitAUSD

Plot1=ggplot(MeltData, aes(GR, FitA, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitALSD, ymax=FitAUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=GRLSD, xmax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=1.9-(0.150*1.0), y=0.6-(0.050*1.2), label="NS", color="steelblue2", size=5) +
  geom_text(x=1.9-(0.050*1.0), y=0.6-(0.050*1.2), label="NS", color="tomato3", size=5) +
  ylab(expression('Fitness'~'('*ln~cells~mL^-1*')')) +
  xlab(expression(atop('Maximum growth rate',paste('('*day^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(1.0,1.9,by=0.3), limits=c(1.0,1.9)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitA100"="Low nutrient – No predation","FitA101"="Low nutrient – Predation","FitA1000"="High nutrient – No predation","FitA1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Fitness vs Growth Rate.tiff', units="in", width=5.2, height=18, res=1000)
Plot1
dev.off()

Plot2=ggplot(MeltData, aes(CC/10^5, FitA, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitALSD, ymax=FitAUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CCLSD/10^5, xmax=CCUSD/10^5, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=5.4-(0.150*5.4), y=0.6-(0.050*1.2), label="NS", color="steelblue2", size=5) +
  geom_text(x=5.4-(0.050*5.4), y=0.6-(0.050*1.2), label="NS", color="tomato3", size=5) +
  ylab(expression('Fitness'~'('*ln~cells~mL^-1*')')) +
  xlab(expression(atop('Carrying capacity',paste('('*10^5~cells~mL^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,5.4,by=1.8), limits=c(0,5.4)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitA100"="Low nutrient – No predation","FitA101"="Low nutrient – Predation","FitA1000"="High nutrient – No predation","FitA1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Fitness vs Carrying Capacity.tiff', units="in", width=5.3, height=18, res=1000)
Plot2
dev.off()

Plot3=ggplot(MeltData, aes(HS, FitA, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitALSD, ymax=FitAUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD/10^1, xmax=HSUSD/10^1, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.150*0.3), y=0.6-(0.050*1.2), label="NS", color="steelblue2", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=0.6-(0.050*1.2), label="NS", color="tomato3", size=5) +
  ylab(expression('Fitness'~'('*ln~cells~mL^-1*')')) +
  xlab(expression(atop('Half-saturation constant',paste('('*µmol~PO[4]^{"-"}~mL^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitA100"="Low nutrient – No predation","FitA101"="Low nutrient – Predation","FitA1000"="High nutrient – No predation","FitA1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,3.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Fitness vs Half Saturation.tiff', units="in", width=5.2, height=18, res=1000)
Plot3
dev.off()

Plot4=ggplot(MeltData, aes(CR, FitA), group=Species) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitALSD, ymax=FitAUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CRLSD, xmax=CRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.021-(0.150*0.021), y=0.6-(0.050*1.2), label="NS", color="steelblue2", size=5) +
  geom_text(x=0.021-(0.050*0.021), y=0.6-(0.050*1.2), label="NS", color="tomato3", size=5) +
  ylab(expression('Fitness'~'('*ln~cells~mL^-1*')')) +
  xlab(expression(atop('Maximum clearance rate',paste('('~10^2~mL~sec^-1~ind^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6)) +
  scale_x_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(0,0.021)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.5) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitA100"="Low nutrient – No predation","FitA101"="Low nutrient – Predation","FitA1000"="High nutrient – No predation","FitA1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Fitness vs Clearance Rate.tiff', units="in", width=5.2, height=18, res=1000)
Plot4
dev.off()

Plot5=ggplot(MeltData, aes(CD/10^1, FitA, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitALSD, ymax=FitAUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CDLSD/10^1, xmax=CDUSD/10^1, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=2.8-(0.150*2.1), y=0.6-(0.050*1.2), label="NS", color="steelblue2", size=5) +
  geom_text(x=2.8-(0.050*2.1), y=0.6-(0.050*1.2), label="NS", color="tomato3", size=5) +
  ylab(expression('Fitness'~'('*ln~cells~mL^-1*')')) +
  xlab(expression(atop('Particle diameter',paste('('*10^1~µm*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0.7,2.8,by=0.7), limits=c(0.7,2.8)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitA100"="Low nutrient – No predation","FitA101"="Low nutrient – Predation","FitA1000"="High nutrient – No predation","FitA1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Fitness vs Cell Diameter.tiff', units="in", width=5.2, height=18, res=1000)
Plot5
dev.off()


# Panel plot of fitness spaces
tiff('Fitness Spaces.tiff', units="in", width=20, height=15, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4,Plot5)
Panel[[1]]=Panel[[1]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[2]]=Panel[[2]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[3]]=Panel[[3]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[4]]=Panel[[4]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[5]]=Panel[[5]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Yaxis=textGrob(expression('Fitness'~'('*ln~cells~mL^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
grid.arrange(grobs=Panel, ncol=5, nrow=1, left=Yaxis, right="", layout_matrix=rbind(c(1,1,2,2,3,3,4,4,5,5)))
grid.text("Low nutrient – No predation", x=unit(50,"cm"), y=unit(33.6,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
grid.text("Low nutrient – Predation", x=unit(50,"cm"), y=unit(24.9,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
grid.text("High nutrient – No predation", x=unit(50,"cm"), y=unit(16.1,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
grid.text("High nutrient – Predation", x=unit(50,"cm"), y=unit(7.4,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
dev.off()

tiff('Fitness Spaces.tiff', units="in", width=16.3, height=15, res=1000)
Panel=list(Plot1,Plot2,Plot4,Plot5)
Panel[[1]]=Panel[[1]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[2]]=Panel[[2]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[3]]=Panel[[3]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[4]]=Panel[[4]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Yaxis=textGrob(expression('Fitness'~'('*ln~cells~mL^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
grid.arrange(grobs=Panel, ncol=5, nrow=1, left=Yaxis, right="", layout_matrix=rbind(c(1,1,2,2,3,3,4,4)))
grid.text("Low nutrient – No predation", x=unit(50,"cm"), y=unit(33.6,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
grid.text("Low nutrient – Predation", x=unit(50,"cm"), y=unit(24.9,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
grid.text("High nutrient – No predation", x=unit(50,"cm"), y=unit(16.1,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
grid.text("High nutrient – Predation", x=unit(50,"cm"), y=unit(7.4,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
dev.off()


##################
### Statistics ###
##################

### Clearance rate ###

# Test for normality
shapiro.test(Data$CR)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$CR~Data$Species)

### Cell diameter ###

# Test for normality
shapiro.test(Data$CA)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$CA~Data$Species)

### Maximum growth rate ###

# Test for normality
shapiro.test(Data$GR)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$GR~Data$Species)

### Carrying capacity ###

# Test for normality
shapiro.test(Data$CC)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$CC~Data$Species)

### Half saturation ###

# Test for normality
shapiro.test(Data$HS)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$HS~Data$Species)


##############################################
### Correlation between fitness and traits ###
##############################################

# Create a list of traits
DataB=subset(Data, Species=="Chlamydomonas bilatus")
DataN=subset(Data, Species=="Chlamydomonas noctigama")
ListData1=list(DataN[,c(1,2,30,9,12,15,6,27)],DataB[,c(1,2,30,9,12,15,6,27)])
ListData2=list(DataN[,c(1,2,31,9,12,15,6,27)],DataB[,c(1,2,31,9,12,15,6,27)])
ListData3=list(DataN[,c(1,2,32,9,12,15,6,27)],DataB[,c(1,2,32,9,12,15,6,27)])
ListData4=list(DataN[,c(1,2,33,9,12,15,6,27)],DataB[,c(1,2,33,9,12,15,6,27)])

# Calculate correlation coefficients for species
CorFitA100.GR=lapply(ListData1, function(x) {cor.test(x[,3],x[,4], method="spearman")})
CorFitA100.CC=lapply(ListData1, function(x) {cor.test(x[,3],x[,5], method="spearman")})
CorFitA100.HS=lapply(ListData1, function(x) {cor.test(x[,3],x[,6], method="spearman")})
CorFitA100.CR=lapply(ListData1, function(x) {cor.test(x[,3],x[,7], method="spearman")})
CorFitA100.CD=lapply(ListData1, function(x) {cor.test(x[,3],x[,8], method="spearman")})

CorFitA101.GR=lapply(ListData2, function(x) {cor.test(x[,3],x[,4], method="spearman")})
CorFitA101.CC=lapply(ListData2, function(x) {cor.test(x[,3],x[,5], method="spearman")})
CorFitA101.HS=lapply(ListData2, function(x) {cor.test(x[,3],x[,6], method="spearman")})
CorFitA101.CR=lapply(ListData2, function(x) {cor.test(x[,3],x[,7], method="spearman")})
CorFitA101.CD=lapply(ListData2, function(x) {cor.test(x[,3],x[,8], method="spearman")})

CorFitA1000.GR=lapply(ListData3, function(x) {cor.test(x[,3],x[,4], method="spearman")})
CorFitA1000.CC=lapply(ListData3, function(x) {cor.test(x[,3],x[,5], method="spearman")})
CorFitA1000.HS=lapply(ListData3, function(x) {cor.test(x[,3],x[,6], method="spearman")})
CorFitA1000.CR=lapply(ListData3, function(x) {cor.test(x[,3],x[,7], method="spearman")})
CorFitA1000.CD=lapply(ListData3, function(x) {cor.test(x[,3],x[,8], method="spearman")})

CorFitA1001.GR=lapply(ListData4, function(x) {cor.test(x[,3],x[,4], method="spearman")})
CorFitA1001.CC=lapply(ListData4, function(x) {cor.test(x[,3],x[,5], method="spearman")})
CorFitA1001.HS=lapply(ListData4, function(x) {cor.test(x[,3],x[,6], method="spearman")})
CorFitA1001.CR=lapply(ListData4, function(x) {cor.test(x[,3],x[,7], method="spearman")})
CorFitA1001.CD=lapply(ListData4, function(x) {cor.test(x[,3],x[,8], method="spearman")})


######################################
### Statistics interspecific level ###
######################################

# Melt the dataset
Data2=melt(Data[,-c(34:41)], id.vars=c(names(Data[,c(1:29)])), variable.name="Treatment", value.name="ExpFitA")
Data2$Nutri=rep(rep(unique(Nutri), each=49*2))
Data2$Pred=rep(rep(unique(Pred), each=49),2)

# Calculate variance of inflation
ModA0=glm(ExpFitA~(GR + CC + HS + CR + CD):Species, family=gaussian(link="identity"), data=Data2)
summary(aov(ModA0))
vif(ModA0)

# Calculate correlations
cor(Data2[,c(6,9,12,15,27)])

# Generalized linear regressions
ModA100=glm(FitA100~(CC + HS + CR + CD):Species, family=gaussian(link="identity"), data=Data)
summary(aov(ModA100))

ModA101=glm(FitA101~(CC + HS + CR + CD):Species, family=gaussian(link="identity"), data=Data)
summary(aov(ModA101))

ModA1000=glm(FitA1000~(CC + HS + CR + CD):Species, family=gaussian(link="identity"), data=Data)
summary(aov(ModA1000))

ModA1001=glm(FitA1001~(CC + HS + CR + CD):Species, family=gaussian(link="identity"), data=Data)
summary(aov(ModA1001))

# Generalized linear regression
ModA0=glm(ExpFitA~(CC + HS + CR + CD):Nutri:Pred, family=gaussian(link="identity"), data=Data2)
ModA1=glm(ExpFitA~(CC:CR + CC:CD + HS:CR + HS:CD + CC:HS:CR:CD):Nutri:Pred, family=gaussian(link="identity"), data=Data2)
lrtest(ModA0,ModA1)
summary(aov(ModA1))

# Posthoc pairwise regression
emmeans(ModA1, pairwise~Nutri*Pred)

# Include the predicted fitness
Data2$PredFitA=c(predict(ModA1))


######################################
### Statistics intraspecific level ###
######################################

# Subset the dataset
DataB=subset(Data, Species=="Chlamydomonas bilatus")
DataN=subset(Data, Species=="Chlamydomonas noctigama")

# Melt the dataset
Data3=melt(DataB[,-c(34:41)], id.vars=c(names(DataB[,c(1:29)])), variable.name="Treatment", value.name="ExpFitA")
Data4=melt(DataN[,-c(34:41)], id.vars=c(names(DataN[,c(1:29)])), variable.name="Treatment", value.name="ExpFitA")
Data3$Nutri=rep(rep(unique(Nutri), each=25*2))
Data4$Nutri=rep(rep(unique(Nutri), each=24*2))
Data3$Pred=rep(rep(unique(Pred), each=25),2)
Data4$Pred=rep(rep(unique(Pred), each=24),2)

# Calculate variance of inflation
ModA0B=glm(ExpFitA~(GR + CC + HS + CR + CD), family=gaussian(link="identity"), data=Data3)
summary(aov(ModA0B))
vif(ModA0B)

ModA0N=glm(ExpFitA~(GR + CC + HS + CR + CD), family=gaussian(link="identity"), data=Data4)
summary(aov(ModA0N))
vif(ModA0N)

# Calculate correlations
cor(Data3[,c(6,9,12,15,27)])
cor(Data4[,c(6,9,12,15,27)])

# Generalized linear regressions
ModA100B=glm(FitA100~(CC + HS + CR + CD), family=gaussian(link="identity"), data=DataB)
summary(aov(ModA100B))
ModA100N=glm(FitA100~(CC + HS + CR + CD), family=gaussian(link="identity"), data=DataN)
summary(aov(ModA100N))

ModA101B=glm(FitA101~(CC + HS + CR + CD), family=gaussian(link="identity"), data=DataB)
summary(aov(ModA101B))
ModA101N=glm(FitA101~(CC + HS + CR + CD), family=gaussian(link="identity"), data=DataN)
summary(aov(ModA101N))

ModA1000B=glm(FitA1000~(CC + HS + CR + CD), family=gaussian(link="identity"), data=DataB)
summary(aov(ModA1000B))
ModA1000N=glm(FitA1000~(CC + HS + CR + CD), family=gaussian(link="identity"), data=DataN)
summary(aov(ModA1000N))

ModA1001B=glm(FitA1001~(CC + HS + CR + CD), family=gaussian(link="identity"), data=DataB)
summary(aov(ModA1001B))
ModA1001N=glm(FitA1001~(CC + HS + CR + CD), family=gaussian(link="identity"), data=DataN)
summary(aov(ModA1001N))

# Generalized linear regression
ModA0B=glm(ExpFitA~(CC + HS + CR + CD):Nutri:Pred, family=gaussian(link="identity"), data=Data3)
ModA1B=glm(ExpFitA~(CC:CR + CC:CD + HS:CR + HS:CD + CC:HS:CR:CD):Nutri:Pred, family=gaussian(link="identity"), data=Data3)
lrtest(ModA0B,ModA1B)
summary(aov(ModA1B))

ModA0N=glm(ExpFitA~(CC + HS + CR + CD):Nutri:Pred, family=gaussian(link="identity"), data=Data4)
ModA1N=glm(ExpFitA~(CC:CR + CC:CD + HS:CR + HS:CD + CC:HS:CR:CD):Nutri:Pred, family=gaussian(link="identity"), data=Data4)
lrtest(ModA0N,ModA1N)
summary(aov(ModA1N))

# Posthoc pairwise regression
emmeans(ModA1B, pairwise~Nutri*Pred)
emmeans(ModA1N, pairwise~Nutri*Pred)

# Include the predicted fitness
Data3$PredFitA=c(predict(ModA1B))
Data4$PredFitA=c(predict(ModA1N))


##########################################
### Plot observed vs predicted fitness ###
##########################################

tiff('Observed vs Predicted Fitness.tiff', units="in", width=12, height=12, res=1000)
ggplot(Data2, aes(ExpFitA, PredFitA)) +
  geom_abline(intercept=0, slope=1, color="black", linetype="11", size=1.2) +
  geom_point(aes(color=Species), pch=16, size=4, alpha=0.8) +
  ylab(expression('Observed fitness'~'('*ln~cells~mL^-1~day^-1*')')) +
  xlab(expression('Predicted fitness'~'('*ln~cells~mL^-1~day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(strip.text.x=element_text(face="bold", color="black", size=16)) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=16)) +
  facet_wrap(~Nutri*Pred, ncol=2, nrow=2, labeller=as_labeller(Labels)) +
  facet_wrap(~Treatment, ncol=2, nrow=2, strip.position="right", labeller=as_labeller(c("FitA100"="Low nutrient – No predation","FitA101"="Low nutrient – Predation","FitA1000"="High nutrient – No predation","FitA1001"="High nutrient – Predation"))) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.2) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.2) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()
