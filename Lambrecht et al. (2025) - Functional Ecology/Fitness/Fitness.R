setwd("~/Activité Professionnelle/LIMNO 2019-2023/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Fitness")

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

#########################################################
#########################################################
##### FITNESS OF ISOLATES IN DIFFERENT ENVIRONMENTS #####
#########################################################
#########################################################

# Import the dataset
Data=read.table("Data_F.txt", h=T, dec=",")
Data$Species=gsub("\\_","\\ ", Data$Species)
summary(Data)
names(Data)

# Specify the variables as numeric or factor
Data[,c(7:12,15)] %<>% mutate_if(is.character,as.numeric)
Data[,c(7:12,15)] %<>% mutate_if(is.integer,as.numeric)
Data[,c(7:12,15)] %<>% mutate_if(is.factor,as.numeric)

# Sort the dataset by isolate, nutrient and predator
Data$Iso=factor(Data$Iso, levels=unique(Data$Iso))
Data$Species=factor(Data$Species, levels=unique(Data$Species))
Data$Nutri=factor(Data$Nutri, levels=c("Low","High"))
Data$Pred=factor(Data$Pred, levels=c("Without","With"))
Data=Data[order(Data$Iso,Data$Nutri,Data$Pred,Data$Trial),]

# Calculate the densities
Data$DensA=Data$IXM*Data$Volu*Data$Site*Data$Dilu*Data$Cove
Data$DensR=Data$Rotifer/5

# Calculate the mean densities
Data2=setDT(Data)[, .(MeanDensA=round(mean(DensA),0), MeanDensR=round(mean(DensR),0)), by=list(Iso,Species,Period,Nutri,Pred,Day)]
Data2=as.data.frame(Data2)
Data=as.data.frame(Data)


###############################
### Plot the fitness curves ###
###############################

tiff('Fitness Curves Low.tiff', units="in", width=15, height=15, res=1000)
ggplot(subset(Data2, Nutri=="Low"), aes(Day, log(MeanDensA+1), group=Pred)) +
  geom_smooth(aes(color=Species, linetype=Pred), size=1, se=F) +
  geom_point(data=subset(Data, Nutri=="Low"), aes(Day, log(DensA+1), color=Species), alpha=0.4, size=1.5, pch=16, position=position_jitter(0.2)) +
  ylab(expression('Cell density'~'('*ln~cells~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=16)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=16)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,14,by=3), limits=c(5,14)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=2), limits=c(0,10)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3"))+
  scale_linetype_manual(values=c("Without"="11", "With"="solid")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.text.x=element_text(face="bold", color="black", size=16)) +
  theme(strip.background=element_blank()) +
  facet_wrap(~Iso, scales="free", ncol=7, nrow=7) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()

tiff('Fitness Curves High.tiff', units="in", width=15, height=15, res=1000)
ggplot(subset(Data2, Nutri=="High"), aes(Day, log(MeanDensA+1), group=Pred)) +
  geom_smooth(aes(color=Species, linetype=Pred), size=1, se=F) +
  geom_point(data=subset(Data, Nutri=="High"), aes(Day, log(DensA+1), color=Species), alpha=0.4, size=1.5, pch=16, position=position_jitter(0.2)) +
  ylab(expression('Cell density'~'('*ln~cells~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=16)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=16)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,14,by=3), limits=c(5,14)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=2), limits=c(0,10)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3"))+
  scale_linetype_manual(values=c("Without"="11", "With"="solid")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.text.x=element_text(face="bold", color="black", size=16)) +
  theme(strip.background=element_blank()) +
  facet_wrap(~Iso, scales="free", ncol=7, nrow=7) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()


###############################
### Calculating the fitness ###
###############################

### Fitness ###

# Calculate the areas under the curve
Data3=setDT(subset(Data, Day < 15))[,.(GrowA=as.numeric(round(AUC(Day, log(DensA+1), from=min(Day,na.rm=T), to=max(Day,na.rm=T), method=c("spline"), absolutearea=F, subdivisions=1000, na.rm=T),0)), GrowP=as.numeric(round(AUC(Day, log(DensR+1), from=min(Day,na.rm=T), to=max(Day,na.rm=T), method=c("spline"), absolutearea=F, subdivisions=1000, na.rm=T),0))), by=list(Iso,Species,Period,Nutri,Pred,Phos,Trial)]
Data3=as.data.frame(Data3)

# Calculate the growth rates
Data3=setDT(subset(Data, Day < 15))[,.(GrowA=coef(lm(log(DensA+1)~Day))[2], GrowP=coef(lm(log(DensR+1)~Day))[2]), by=list(Iso,Species,Period,Nutri,Pred,Phos,Trial)]
Data3=as.data.frame(Data3)

# Export the dataset
Data3=Data3[order(Data3$Iso,Data3$Trial),]
write.table(Data3[,c(1:9)], file="Data_Fitness.txt", sep="\t", row.names=F)

### Error bars ###

# Calculate the error bars for growth rate
Data4=setDT(Data3)[, .(GrowA=round(mean(GrowA),4), GrowALSD=round(mean(GrowA)-sd(GrowA),4), GrowAUSD=round(mean(GrowA)+sd(GrowA),4), GrowP=round(mean(GrowP),4), GrowPLSD=round(mean(GrowP)-sd(GrowP),4), GrowPUSD=round(mean(GrowP)+sd(GrowP),4)), by=list(Iso,Species,Nutri,Pred)]
Data4=as.data.frame(Data4)

# Calculate the error bars for growth rate
Data5=as.data.frame(setDT(na.omit(Data4))[, .(GrowASpeciesL=round(mean(GrowA)-sd(GrowA),4), GrowASpecies=round(mean(GrowA),4), GrowASpeciesU=round(mean(GrowA)+sd(GrowA),4)), by=list(Species,Nutri,Pred)])
Data4$GrowASpeciesL=c(rep(Data5$GrowASpeciesL[c(1:4)],24), rep(Data5$GrowASpeciesL[c(5:8)],25))
Data4$GrowASpecies=c(rep(Data5$GrowASpecies[c(1:4)],24), rep(Data5$GrowASpecies[c(5:8)],25))
Data4$GrowASpeciesU=c(rep(Data5$GrowASpeciesU[c(1:4)],24), rep(Data5$GrowASpeciesU[c(5:8)],25))


######################################
### Plot the areas under the curve ###
######################################

# Combine the treatments
Data3$Treat=paste(Data3$Nutri, Data3$Pred, sep="-")
Data4$Treat=paste(Data4$Nutri, Data4$Pred, sep="-")
Data4$Treat=paste(Data4$Treat, Data4$Species, sep="-")
Data4$Treat=gsub("Chlamydomonas bilatus", "CB", Data4$Treat)
Data4$Treat=gsub("Chlamydomonas noctigama", "CN", Data4$Treat)
Data3$Treat=factor(Data3$Treat, levels=c("Low-Without","Low-With","High-Without","High-With"))
Data4$Treat=factor(Data4$Treat, levels=c("Low-Without-CN","Low-Without-CB","Low-With-CN","Low-With-CB","High-Without-CN","High-Without-CB","High-With-CN","High-With-CB"))

# Calculate the mean areas under the curve
Data6=setDT(Data4)[, .(MeanGrowA=round(mean(GrowA),4), MeanGrowP=round(mean(GrowP),4)), by=list(Species,Nutri,Pred,Treat)]
Data6=as.data.frame(Data6)

# Calculate the significance between species
Stats1=Data3 %>% group_by(Treat) %>% wilcox_test(GrowA~Species) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats1=cbind(Stats1[,c(1,9:10)], x=c(1.00,3.00,5.00,7.00), xend=c(2.00,4.00,6.00,8.00), y=rep(0.6+(0.020*0.6),4), yend=rep(0.6+(0.020*0.6),4))
colnames(Stats1)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between predator
Stats2=Data3 %>% group_by(Nutri) %>% wilcox_test(GrowA~Pred) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats2=cbind(Stats2[,c(1,9:10)], x=c(1.50,5.50), xend=c(3.50,7.50), y=rep(0.6+(0.145*0.6),2), yend=rep(0.6+(0.145*0.6),2))
colnames(Stats2)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between nutrient
Stats3=Data3 %>% group_by(Pred) %>% wilcox_test(GrowA~Nutri) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats3=cbind(Stats3[,c(1,9:10)], x=c(2.50), xend=c(6.50), y=rep(0.6+(0.270*0.6),1), yend=rep(0.6+(0.270*0.6),1))
colnames(Stats3)[c(1:3)]=c("Trait","P","S")

# Calculate coefficients of variation
DataCV=setDT(Data3)[, .(GrowA=abs(round((sd(GrowA)/mean(GrowA))*100,1))), by=list(Species,Nutri,Pred)]
DataCV=as.data.frame(DataCV)
Data3=as.data.frame(Data3)

tiff('Fitness.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data4, aes(Treat, GrowA)) + coord_cartesian(clip="off") +
  geom_pointrange(aes(Treat, ymin=GrowALSD, ymax=GrowAUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_text(data=DataCV[1,], mapping=aes(y=-0.6-(-0.6-0.6)*0.040, x=1.00, label=paste(format(GrowA,nsmall=1))), color="steelblue2", size=5, hjust=0.5) +
  geom_text(data=DataCV[2,], mapping=aes(y=-0.6-(-0.6-0.6)*0.040, x=2.00, label=paste(format(GrowA,nsmall=1))), color="tomato3", size=5, hjust=0.5) +
  geom_text(data=DataCV[3,], mapping=aes(y=-0.6-(-0.6-0.6)*0.040, x=3.00, label=paste(format(GrowA,nsmall=1))), color="steelblue2", size=5, hjust=0.5) +
  geom_text(data=DataCV[4,], mapping=aes(y=-0.6-(-0.6-0.6)*0.040, x=4.00, label=paste(format(GrowA,nsmall=1))), color="tomato3", size=5, hjust=0.5) +
  geom_text(data=DataCV[5,], mapping=aes(y=-0.6-(-0.6-0.6)*0.040, x=5.00, label=paste(format(GrowA,nsmall=1))), color="steelblue2", size=5, hjust=0.5) +
  geom_text(data=DataCV[6,], mapping=aes(y=-0.6-(-0.6-0.6)*0.040, x=6.00, label=paste(format(GrowA,nsmall=1))), color="tomato3", size=5, hjust=0.5) +
  geom_text(data=DataCV[7,], mapping=aes(y=-0.6-(-0.6-0.6)*0.040, x=7.00, label=paste(format(GrowA,nsmall=1))), color="steelblue2", size=5, hjust=0.5) +
  geom_text(data=DataCV[8,], mapping=aes(y=-0.6-(-0.6-0.6)*0.040, x=8.00, label=paste(format(GrowA,nsmall=1))), color="tomato3", size=5, hjust=0.5) +
  geom_segment(data=Stats1, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats2, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats3, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_text(data=Stats1, aes(x=x+0.50, y=0.6+(0.050*0.6), label=S), color="black", size=5) +
  geom_text(data=Stats2, aes(x=x+1.00, y=0.6+(0.175*0.6), label=S), color="black", size=5) +
  geom_text(data=Stats3, aes(x=x+2.00, y=0.6+(0.300*0.6), label=S), color="black", size=5) +
  ylab(expression('Fitness'~'('*day^-1*')')) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=1.5, xmax=1.5, ymin=-0.70, ymax=-0.70) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=3.5, xmax=3.5, ymin=-0.70, ymax=-0.70) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=5.5, xmax=5.5, ymin=-0.70, ymax=-0.70) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=7.5, xmax=7.5, ymin=-0.70, ymax=-0.70) +
  annotation_custom(grob=textGrob(label="Low nutrient", gp=gpar(cex=1.5, color="black")), xmin=2.5, xmax=2.5, ymin=-0.77, ymax=-0.77) +
  annotation_custom(grob=textGrob(label="High nutrient", gp=gpar(cex=1.5, color="black")), xmin=6.5, xmax=6.5, ymin=-0.77, ymax=-0.77) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour=NA, size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour=NA, size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6+(0.3*0.6))) +
  scale_x_discrete(labels=c("Low-Without"="No predation", "Low-With"="Predation", "High-Without"="No predation", "High-With"="Predation")) +
  scale_fill_manual(values=alpha(c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3"),0.3)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()


##################
### Statistics ###
##################

# Test for normality
shapiro.test(Data3$GrowA)
qqnorm(Data3$GrowA)
qqline(Data3$GrowA)
hist(Data3$GrowA)

# Test for homoscedasticity
bartlett.test(Data3$GrowA, Data3$Species)
bartlett.test(Data3$GrowA, Data3$Nutri)
bartlett.test(Data3$GrowA, Data3$Pred)

# Treatment differences
summary(aov(glm((GrowA+1)~Species*Nutri*Pred, family=Gamma(link="identity"), data=Data3)))
summary(aov(glm((GrowA+1)~Iso*Nutri*Pred, family=Gamma(link="identity"), data=subset(Data3, Species=="Chlamydomonas bilatus"))))
summary(aov(glm((GrowA+1)~Iso*Nutri*Pred, family=Gamma(link="identity"), data=subset(Data3, Species=="Chlamydomonas noctigama"))))
