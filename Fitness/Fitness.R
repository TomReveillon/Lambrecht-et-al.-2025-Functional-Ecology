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

# Calculate the areas under curve
Data3=setDT(subset(Data, Day < 15))[,.(AucA=as.numeric(round(AUC(Day, log(DensA+1), from=min(Day,na.rm=T), to=max(Day,na.rm=T), method=c("spline"), absolutearea=F, subdivisions=1000, na.rm=T),0)), AucP=as.numeric(round(AUC(Day, log(DensR+1), from=min(Day,na.rm=T), to=max(Day,na.rm=T), method=c("spline"), absolutearea=F, subdivisions=1000, na.rm=T),0))), by=list(Iso,Species,Period,Nutri,Pred,Phos,Trial)]
Data3=as.data.frame(Data3)

# Calculate the growth rates
Data4=setDT(subset(Data, Day < 15))[,.(GrowA=round(coef(lm(log(DensA+1)~Day))[2],6), GrowP=round(coef(lm(log(DensR+1)~Day))[2],6)), by=list(Iso,Species,Period,Nutri,Pred,Phos,Trial)]
Data4=as.data.frame(Data4)

# Export the dataset
Data4=cbind(Data3[,c(1:9)],Data4[,c(8:9)])
Data4=Data4[order(Data4$Iso,Data3$Trial),]
write.table(Data4[,c(1:11)], file="Data_Fitness.txt", sep="\t", row.names=F)

### Error bars ###

# Calculate the error bars
Data5=setDT(Data4)[, .(AucA=round(mean(AucA),0), AucALSD=round(mean(AucA)-sd(AucA),0), AucAUSD=round(mean(AucA)+sd(AucA),0), GrowA=round(mean(GrowA),6), GrowALSD=round(mean(GrowA)-sd(GrowA),6), GrowAUSD=round(mean(GrowA)+sd(GrowA),6)), by=list(Iso,Species,Nutri,Pred)]
Data5=as.data.frame(Data5)

# Calculate the error bars for area under curve
Data6=as.data.frame(setDT(na.omit(Data4))[, .(AucASpeciesL=round(mean(AucA)-sd(AucA),0), AucASpecies=round(mean(AucA),0), AucASpeciesU=round(mean(AucA)+sd(AucA),0)), by=list(Species,Nutri,Pred)])
Data5$AucASpeciesL=c(rep(Data6$AucASpeciesL[c(1:4)],24), rep(Data6$AucASpeciesL[c(5:8)],25))
Data5$AucASpecies=c(rep(Data6$AucASpecies[c(1:4)],24), rep(Data6$AucASpecies[c(5:8)],25))
Data5$AucASpeciesU=c(rep(Data6$AucASpeciesU[c(1:4)],24), rep(Data6$AucASpeciesU[c(5:8)],25))

# Calculate the error bars for growth rate
Data6=as.data.frame(setDT(na.omit(Data4))[, .(GrowASpeciesL=round(mean(GrowA)-sd(GrowA),6), GrowASpecies=round(mean(GrowA),6), GrowASpeciesU=round(mean(GrowA)+sd(GrowA),6)), by=list(Species,Nutri,Pred)])
Data5$GrowASpeciesL=c(rep(Data6$GrowASpeciesL[c(1:4)],24), rep(Data6$GrowASpeciesL[c(5:8)],25))
Data5$GrowASpecies=c(rep(Data6$GrowASpecies[c(1:4)],24), rep(Data6$GrowASpecies[c(5:8)],25))
Data5$GrowASpeciesU=c(rep(Data6$GrowASpeciesU[c(1:4)],24), rep(Data6$GrowASpeciesU[c(5:8)],25))


########################
### Plot the fitness ###
########################

# Combine the treatments
Data4$Treat=paste(Data4$Nutri, Data4$Pred, sep="-")
Data5$Treat=paste(Data5$Nutri, Data5$Pred, sep="-")
Data5$Treat=paste(Data5$Treat, Data5$Species, sep="-")
Data5$Treat=gsub("Chlamydomonas bilatus", "CB", Data5$Treat)
Data5$Treat=gsub("Chlamydomonas noctigama", "CN", Data5$Treat)
Data4$Treat=factor(Data4$Treat, levels=c("Low-Without","Low-With","High-Without","High-With"))
Data5$Treat=factor(Data5$Treat, levels=c("Low-Without-CN","Low-Without-CB","Low-With-CN","Low-With-CB","High-Without-CN","High-Without-CB","High-With-CN","High-With-CB"))

# Calculate the mean areas and growth rates
Data7=setDT(Data5)[, .(MeanAucA=round(mean(AucA),4), MeanGrowA=round(mean(GrowA),4)), by=list(Species,Nutri,Pred,Treat)]
Data7=as.data.frame(Data7)

# Calculate coefficients of variation
DataCV=setDT(Data5)[, .(AucA=abs(round((sd(AucA)/mean(AucA))*100,1)), GrowA=abs(round((sd(GrowA)/mean(GrowA))*100,1))), by=list(Species,Nutri,Pred)]
DataCV=as.data.frame(DataCV)
Data5=as.data.frame(Data5)

### Growth rate ###

# Calculate the significance between species
Stats1=Data4 %>% group_by(Treat) %>% wilcox_test(GrowA~Species) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats1=cbind(Stats1[,c(1,9:10)], x=c(1.00,3.00,5.00,7.00), xend=c(2.00,4.00,6.00,8.00), y=rep(0.6+(0.020*0.6),4), yend=rep(0.6+(0.020*0.6),4))
colnames(Stats1)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between predator
Stats2=Data4 %>% group_by(Nutri) %>% wilcox_test(GrowA~Pred) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats2=cbind(Stats2[,c(1,9:10)], x=c(1.50,5.50), xend=c(3.50,7.50), y=rep(0.6+(0.145*0.6),2), yend=rep(0.6+(0.145*0.6),2))
colnames(Stats2)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between nutrient
Stats3=Data4 %>% group_by(Pred) %>% wilcox_test(GrowA~Nutri) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats3=cbind(Stats3[,c(1,9:10)], x=c(2.50), xend=c(6.50), y=rep(0.6+(0.270*0.6),1), yend=rep(0.6+(0.270*0.6),1))
colnames(Stats3)[c(1:3)]=c("Trait","P","S")

tiff('Fitness.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data5, aes(Treat, GrowA)) + coord_cartesian(clip="off") +
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
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=1.5, xmax=1.5, ymin=-0.703, ymax=-0.703) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=3.5, xmax=3.5, ymin=-0.703, ymax=-0.703) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=5.5, xmax=5.5, ymin=-0.703, ymax=-0.703) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=7.5, xmax=7.5, ymin=-0.703, ymax=-0.703) +
  annotation_custom(grob=textGrob(label="Low nutrient", gp=gpar(cex=1.5, color="black")), xmin=2.5, xmax=2.5, ymin=-0.768, ymax=-0.768) +
  annotation_custom(grob=textGrob(label="High nutrient", gp=gpar(cex=1.5, color="black")), xmin=6.5, xmax=6.5, ymin=-0.768, ymax=-0.768) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour=NA, size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="red", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.6,0.6,by=0.3), limits=c(-0.6,0.6+(0.3*0.6))) +
  scale_x_discrete(labels=c("Low-Without"="No predation", "Low-With"="Predation", "High-Without"="No predation", "High-With"="Predation")) +
  scale_fill_manual(values=alpha(c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3"),0.3)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()

### Areas under curve ###

# Calculate the significance between species
Stats1=Data4 %>% group_by(Treat) %>% wilcox_test(GrowA~Species) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats1=cbind(Stats1[,c(1,9:10)], x=c(1.00,3.00,5.00,7.00), xend=c(2.00,4.00,6.00,8.00), y=rep(125+(0.020*37.5),4), yend=rep(125+(0.020*37.5),4))
colnames(Stats1)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between predator
Stats2=Data4 %>% group_by(Nutri) %>% wilcox_test(GrowA~Pred) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats2=cbind(Stats2[,c(1,9:10)], x=c(1.50,5.50), xend=c(3.50,7.50), y=rep(125+(0.145*37.5),2), yend=rep(125+(0.145*37.5),2))
colnames(Stats2)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between nutrient
Stats3=Data4 %>% group_by(Pred) %>% wilcox_test(GrowA~Nutri) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats3=cbind(Stats3[,c(1,9:10)], x=c(2.50), xend=c(6.50), y=rep(125+(0.270*37.5),1), yend=rep(125+(0.270*37.5),1))
colnames(Stats3)[c(1:3)]=c("Trait","P","S")

tiff('Fitness.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data5, aes(Treat, AucA)) + coord_cartesian(clip="off") +
  geom_pointrange(aes(Treat, ymin=AucALSD, ymax=AucAUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_text(data=DataCV[1,], mapping=aes(y=50+(125-50)*0.040, x=1.00, label=paste(format(AucA,nsmall=1))), color="steelblue2", size=5, hjust=0.5) +
  geom_text(data=DataCV[2,], mapping=aes(y=50+(125-50)*0.040, x=2.00, label=paste(format(AucA,nsmall=1))), color="tomato3", size=5, hjust=0.5) +
  geom_text(data=DataCV[3,], mapping=aes(y=50+(125-50)*0.040, x=3.00, label=paste(format(AucA,nsmall=1))), color="steelblue2", size=5, hjust=0.5) +
  geom_text(data=DataCV[4,], mapping=aes(y=50+(125-50)*0.040, x=4.00, label=paste(format(AucA,nsmall=1))), color="tomato3", size=5, hjust=0.5) +
  geom_text(data=DataCV[5,], mapping=aes(y=50+(125-50)*0.040, x=5.00, label=paste(format(AucA,nsmall=1))), color="steelblue2", size=5, hjust=0.5) +
  geom_text(data=DataCV[6,], mapping=aes(y=50+(125-50)*0.040, x=6.00, label=paste(format(AucA,nsmall=1))), color="tomato3", size=5, hjust=0.5) +
  geom_text(data=DataCV[7,], mapping=aes(y=50+(125-50)*0.040, x=7.00, label=paste(format(AucA,nsmall=1))), color="steelblue2", size=5, hjust=0.5) +
  geom_text(data=DataCV[8,], mapping=aes(y=50+(125-50)*0.040, x=8.00, label=paste(format(AucA,nsmall=1))), color="tomato3", size=5, hjust=0.5) +
  geom_segment(data=Stats1, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats2, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats3, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_text(data=Stats1, aes(x=x+0.50, y=125+(0.050*37.5), label=S), color="black", size=5) +
  geom_text(data=Stats2, aes(x=x+1.00, y=125+(0.175*37.5), label=S), color="black", size=5) +
  geom_text(data=Stats3, aes(x=x+2.00, y=125+(0.300*37.5), label=S), color="black", size=5) +
  ylab(expression('Fitness'~'('*ln~cells~day^-1*')')) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=1.5, xmax=1.5, ymin=43.5, ymax=43.5) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=3.5, xmax=3.5, ymin=43.5, ymax=43.5) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=5.5, xmax=5.5, ymin=43.5, ymax=43.5) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=7.5, xmax=7.5, ymin=43.5, ymax=43.5) +
  annotation_custom(grob=textGrob(label="Low nutrient", gp=gpar(cex=1.5, color="black")), xmin=2.5, xmax=2.5, ymin=39.5, ymax=39.5) +
  annotation_custom(grob=textGrob(label="High nutrient", gp=gpar(cex=1.5, color="black")), xmin=6.5, xmax=6.5, ymin=39.5, ymax=39.5) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour=NA, size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour=NA, size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(50,125,by=25), limits=c(50,125+(0.3*37.5))) +
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
