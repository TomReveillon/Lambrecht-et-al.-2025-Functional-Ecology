setwd("~/Activité Professionnelle/LIMNO 2019-2023/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Trait Space")

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

####################################
####################################
##### TRAIT SPACES OF ISOLATES #####
####################################
####################################

# Import the dataset for ingestion rate 
DataIR=read.table("Data_IR_TO.txt", h=T, dec=".")
summary(DataIR)
names(DataIR)

# Import the dataset for growth rate
DataGR=read.table("Data_GR_TO.txt", h=T, dec=".")
summary(DataGR)
names(DataGR)

# Import the dataset for half saturation constant
DataHS=read.table("Data_HS_TO.txt", h=T, dec=".")
summary(DataHS)
names(DataHS)

# Import the dataset for morphology
DataMO=read.table("Data_MO_TO.txt", h=T, dec=".")
summary(DataMO)
names(DataMO)

# Specify the variables as numeric or factor
DataIR[,c(3:14)] %<>% mutate_if(is.character,as.numeric)
DataGR[,c(3:14)] %<>% mutate_if(is.character,as.numeric)
DataHS[,c(3:20)] %<>% mutate_if(is.character,as.numeric)
DataMO[,c(3:14)] %<>% mutate_if(is.character,as.numeric)

# Sort the datasets by isolate
DataIR$Iso=factor(DataIR$Iso, levels=unique(DataIR$Iso))
DataGR$Iso=factor(DataGR$Iso, levels=unique(DataGR$Iso))
DataHS$Iso=factor(DataHS$Iso, levels=unique(DataHS$Iso))
DataMO$Iso=factor(DataMO$Iso, levels=unique(DataMO$Iso))

# Combine the datasets
Data=cbind(DataIR[,c(1:14)],DataGR[,c(3:14)],DataHS[,c(3:20)],DataMO[,c(3:14)])
Data=Data[,c(1:8,15:20,27:35,45:50,9:14,21:26,36:44,51:56)]

# Calculate coefficients of variation
DataCV=setDT(Data)[, .(CR=round((sd(CR)/mean(CR))*100,1), GR=round((sd(GR)/mean(GR))*100,1), CC=round((sd(CC)/mean(CC))*100,1), HS=round((sd(HS)/mean(HS))*100,1), CD=round((sd(CD)/mean(CD))*100,1)), by=list(Species)]
DataCV=as.data.frame(DataCV)
Data=as.data.frame(Data)


###############################################
### Plot defense and competitiveness traits ###
###############################################

Plot1=ggplot(Data, aes(Species, GR, group=Species)) +
  geom_pointrange(aes(Species, ymin=GRLSD, ymax=GRUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_text(data=DataCV[2,], mapping=aes(y=1.9-(1.9-1.0)*0.040, x=0.50, color=Species, label=paste(format(GR,nsmall=1))), size=5, hjust=0) +
  geom_text(data=DataCV[1,], mapping=aes(y=1.9-(1.9-1.0)*0.040, x=2.50, color=Species, label=paste(format(GR,nsmall=1))), size=5, hjust=1) +
  geom_segment(x=0.70, xend=1.30, y=1.9+(0.020*0.9), yend=1.9+(0.020*0.9), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=1.9+(0.020*0.9), yend=1.9+(0.020*0.9), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=1.9+(0.100*0.9), yend=1.9+(0.100*0.9), color="black", size=0.8) +
  geom_text(x=1.00, y=1.9+(0.035*0.9), label="*", color="black", size=5) +
  geom_text(x=2.00, y=1.9+(0.035*0.9), label="*", color="black", size=5) +
  geom_text(x=1.50, y=1.9+(0.125*0.9), label="NS", color="black", size=5) +
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(1.0,1.9,by=0.3),fmt="%.1f"), breaks=seq(1.0,1.9,by=0.3), limits=c(1.0,1.9+(0.10*0.9))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot2=ggplot(Data, aes(Species, CC, group=Species)) +
  geom_pointrange(aes(Species, ymin=CCLSD, ymax=CCUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_text(data=DataCV[2,], mapping=aes(y=5.4*10^5-(5.4*10^5-0)*0.040, x=0.50, color=Species, label=paste(format(CC,nsmall=1))), size=5, hjust=0) +
  geom_text(data=DataCV[1,], mapping=aes(y=5.4*10^5-(5.4*10^5-0)*0.040, x=2.50, color=Species, label=paste(format(CC,nsmall=1))), size=5, hjust=1) +
  geom_segment(x=0.70, xend=1.30, y=5.4*10^5+(0.020*5.4*10^5), yend=5.4*10^5+(0.020*5.4*10^5), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=5.4*10^5+(0.020*5.4*10^5), yend=5.4*10^5+(0.020*5.4*10^5), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=5.4*10^5+(0.100*5.4*10^5), yend=5.4*10^5+(0.100*5.4*10^5), color="black", size=0.8) +
  geom_text(x=1.00, y=5.4*10^5+(0.035*5.4*10^5), label="***", color="black", size=5) +
  geom_text(x=2.00, y=5.4*10^5+(0.035*5.4*10^5), label="***", color="black", size=5) +
  geom_text(x=1.50, y=5.4*10^5+(0.115*5.4*10^5), label="***", color="black", size=5) +
  ylab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5+(0.10*5.4*10^5))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot3=ggplot(Data, aes(Species, HS, group=Species)) +
  geom_pointrange(aes(Species, ymin=HSLSD, ymax=HSUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_text(data=DataCV[2,], mapping=aes(y=0.6-(0.6-0)*0.040, x=0.50, color=Species, label=paste(format(HS,nsmall=1))), size=5, hjust=0) +
  geom_text(data=DataCV[1,], mapping=aes(y=0.6-(0.6-0)*0.040, x=2.50, color=Species, label=paste(format(HS,nsmall=1))), size=5, hjust=1) +
  geom_segment(x=0.70, xend=1.30, y=0.6+(0.020*0.6), yend=0.6+(0.020*0.6), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=0.6+(0.020*0.6), yend=0.6+(0.020*0.6), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=0.6+(0.100*0.6), yend=0.6+(0.100*0.6), color="black", size=0.8) +
  geom_text(x=1.00, y=0.6+(0.045*0.6), label="**", color="black", size=5) +
  geom_text(x=2.00, y=0.6+(0.035*0.6), label="***", color="black", size=5) +
  geom_text(x=1.50, y=0.6+(0.125*0.6), label="NS", color="black", size=5) +
  ylab(expression('Half-saturation constant'~'('*µmol~PO[4]^{"-"}~mL^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(0,0.6,by=0.2),fmt="%.1f"), breaks=seq(0,0.6,by=0.2), limits=c(-10^-5,0.6+(0.10*0.6))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot4=ggplot(Data, aes(Species, CR, group=Species)) +
  geom_pointrange(aes(Species, ymin=CRLSD, ymax=CRUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_text(data=DataCV[2,], mapping=aes(y=0.021-(0.021-0)*0.040, x=0.50, color=Species, label=paste(format(CR,nsmall=1))), size=5, hjust=0) +
  geom_text(data=DataCV[1,], mapping=aes(y=0.021-(0.021-0)*0.040, x=2.50, color=Species, label=paste(format(CR,nsmall=1))), size=5, hjust=1) +
  geom_segment(x=0.70, xend=1.30, y=0.021+(0.020*0.021), yend=0.021+(0.020*0.021), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=0.021+(0.020*0.021), yend=0.021+(0.020*0.021), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=0.021+(0.100*0.021), yend=0.021+(0.100*0.021), color="black", size=0.8) +
  geom_text(x=1.00, y=0.021+(0.035*0.021), label="***", color="black", size=5) +
  geom_text(x=2.00, y=0.021+(0.045*0.021), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=0.021+(0.125*0.021), label="NS", color="black", size=5) +
  ylab(expression('Maximum clearance rate'~'('*10^-2~mL~day^-1~ind^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(-10^-5,0.021+(0.10*0.021))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot5=ggplot(Data, aes(Species, CD, group=Species)) +
  geom_pointrange(aes(Species, ymin=CDLSD, ymax=CDUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_text(data=DataCV[2,], mapping=aes(y=2.8*10^1-(2.8*10^1-0.7*10^1)*0.040, x=0.50, color=Species, label=paste(format(CD,nsmall=1))), size=5, hjust=0) +
  geom_text(data=DataCV[1,], mapping=aes(y=2.8*10^1-(2.8*10^1-0.7*10^1)*0.040, x=2.50, color=Species, label=paste(format(CD,nsmall=1))), size=5, hjust=1) +
  geom_segment(x=0.70, xend=1.30, y=2.8*10^1+(0.020*2.1*10^1), yend=2.8*10^1+(0.020*2.1*10^1), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=2.8*10^1+(0.020*2.1*10^1), yend=2.8*10^1+(0.020*2.1*10^1), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=2.8*10^1+(0.100*2.1*10^1), yend=2.8*10^1+(0.100*2.1*10^1), color="black", size=0.8) +
  geom_text(x=1.00, y=2.8*10^1+(0.045*2.1*10^1), label="NS", color="black", size=5) +
  geom_text(x=2.00, y=2.8*10^1+(0.045*2.1*10^1), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=2.8*10^1+(0.115*2.1*10^1), label="**", color="black", size=5) +
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression(italic('Chlamydomonas')~'spp. isolates')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(0.7,2.8,by=0.7),fmt="%.1f"), breaks=seq(0.7*10^1,2.8*10^1,by=0.7*10^1), limits=c(0.7*10^1,2.8*10^1+(0.10*2.1*10^1))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,16.5),"pt")) +
  theme(legend.position="none")

# Panel plot of traits
tiff('Traits.tiff', units="in", width=12, height=17, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4,Plot5)
grid.arrange(grobs=Panel, ncol=2, nrow=3, layout_matrix=rbind(c(1,1,1,1,2,2,2,2),c(3,3,3,3,4,4,4,4),c(NA,NA,5,5,5,5,NA,NA)))
dev.off()


#####################################
### Calculating correlation lines ###
#####################################

# Export the dataset
Data=Data[,c(1:29)]
write.table(Data[,c(1:29)], file="Data_Traits.txt", sep="\t", row.names=F)

# Rescale the dataset
Data[,c(12:14)]=round(Data[,c(12:14)]/10^5,4)

### Clearance rate (defense) vs growth rate (competitiveness) ###

# Fitting linear models
ModLi=lm(CR~GR, data=Data)
ModLiN=lm(CR~GR, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(CR~GR, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(CR~max(CR) + a*(min(GR)+GR)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(CR~max(CR) + a*(min(GR)+GR)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(CR~max(CR) + a*(min(GR)+GR)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.CR.GR=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.CR.GR.All=round(c(predict(ModPo)),4)

### Clearance rate (defense) vs carrying capacity (competitiveness) ###

# Fitting linear models
ModLi=lm(CR~CC, data=Data)
ModLiN=lm(CR~CC, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(CR~CC, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(CR~max(CR) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(CR~max(CR) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(CR~max(CR) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.CR.CC=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.CR.CC.All=round(c(predict(ModPo)),4)

### Clearance rate (defense) vs half saturation (competitiveness) ###

# Fitting linear models
ModLi=lm(CR~HS, data=Data)
ModLiN=lm(CR~HS, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(CR~HS, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(CR~max(CR) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(CR~max(CR) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(CR~max(CR) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.CR.HS=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.CR.HS.All=round(c(predict(ModPo)),4)

### Cell diameter (defense) vs growth rate (competitiveness) ###

# Fitting linear models
ModLi=lm(CD~GR, data=Data)
ModLiN=lm(CD~GR, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(CD~GR, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(CD~max(CD) + a*(min(GR)+GR)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(CD~max(CD) + a*(min(GR)+GR)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(CD~max(CD) + a*(min(GR)+GR)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.CD.GR=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.CD.GR.All=round(c(predict(ModPo)),4)

### Cell area (defense) vs carrying capacity (competitiveness) ###

# Fitting linear models
ModLi=lm(CD~CC, data=Data)
ModLiN=lm(CD~CC, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(CD~CC, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(CD~max(CD) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(CD~max(CD) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(CD~max(CD) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.CD.CC=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.CD.CC.All=round(c(predict(ModPo)),4)

### Cell area (defense) vs half saturation (competitiveness) ###

# Fitting linear models
ModLi=lm(CD~HS, data=Data)
ModLiN=lm(CD~HS, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(CD~HS, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(CD~max(CD) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(CD~max(CD) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(CD~max(CD) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.CD.HS=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.CD.HS.All=round(c(predict(ModPo)),4)

### Growth rate (competitiveness) vs carrying capacity (competitiveness) ###

# Fitting linear models
ModLi=lm(GR~CC, data=Data)
ModLiN=lm(GR~CC, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(GR~CC, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(GR~max(GR) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(GR~max(GR) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(GR~max(GR) + a*(min(CC)+CC)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.GR.CC=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.GR.CC.All=round(c(predict(ModPo)),4)

### Growth rate (competitiveness) vs half saturation (competitiveness) ###

# Fitting linear models
ModLi=lm(GR~HS, data=Data)
ModLiN=lm(GR~HS, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(GR~HS, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(GR~max(GR) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(GR~max(GR) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(GR~max(GR) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.GR.HS=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.GR.HS.All=round(c(predict(ModPo)),4)

### Carrying capacity (competitiveness) vs half saturation (competitiveness) ###

# Fitting linear models
ModLi=lm(CC~HS, data=Data)
ModLiN=lm(CC~HS, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModLiB=lm(CC~HS, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModPo=nls(CC~max(CC) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=Data)
ModPoN=nls(CC~max(CC) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModPoB=nls(CC~max(CC) + a*(min(HS)+HS)^b, start=c(a=1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Compare linear and polynomial models
akaike.weights(c(AIC(ModLi),AIC(ModPo)))
akaike.weights(c(AIC(ModLiN),AIC(ModPoN)))
akaike.weights(c(AIC(ModLiB),AIC(ModPoB)))

# Calculate the predicted values
Data$TO.CC.HS=round(c(predict(ModLiN),predict(ModLiB)),4)
Data$TO.CC.HS.All=round(c(predict(ModPo)),4)


##################################
### Correlation between traits ###
##################################

# Rescale the dataset
Data[,c(12:14,46:47)]=round(Data[,c(12:14,46:47)]*10^5,4)

# Test for normality
shapiro.test(Data$CR)
shapiro.test(Data$GR)
shapiro.test(Data$CC)
shapiro.test(Data$HS)
shapiro.test(Data$CD)

# Create a list of traits
DataB=subset(Data, Species=="Chlamydomonas bilatus")
DataN=subset(Data, Species=="Chlamydomonas noctigama")
ListData1=list(DataN[,c(1,2,6,9,12,15)],DataB[,c(1,2,6,9,12,15)])
ListData2=list(DataN[,c(1,2,27,9,12,15)],DataB[,c(1,2,27,9,12,15)])
ListData3=list(DataN[,c(1,2,9,12,15)],DataB[,c(1,2,9,12,15)])

# Calculate correlation coefficients for species
CorCR.GR=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorCR.GR=c(rep(CorCR.GR[1],24),rep(CorCR.GR[2],25))
CorCR.CC=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorCR.CC=c(rep(CorCR.CC[1],24),rep(CorCR.CC[2],25))
CorCR.HS=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,6], method="spearman")[[4]]})),4)
Data$CorCR.HS=c(rep(CorCR.HS[1],24),rep(CorCR.HS[2],25))

CorCD.GR=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorCD.GR=c(rep(CorCD.GR[1],24),rep(CorCD.GR[2],25))
CorCD.CC=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorCD.CC=c(rep(CorCD.CC[1],24),rep(CorCD.CC[2],25))
CorCD.HS=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,6], method="spearman")[[4]]})),4)
Data$CorCD.HS=c(rep(CorCD.HS[1],24),rep(CorCD.HS[2],25))

CorGR.CC=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorGR.CC=c(rep(CorGR.CC[1],24),rep(CorGR.CC[2],25))
CorGR.HS=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorGR.HS=c(rep(CorGR.HS[1],24),rep(CorGR.HS[2],25))
CorCC.HS=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,4],x[,5], method="spearman")[[4]]})),4)
Data$CorCC.HS=c(rep(CorCC.HS[1],24),rep(CorCC.HS[2],25))

# Calculate correlation coefficients for genus
CorCR.GR.All=round(do.call("rbind",lapply(list(bind_rows(ListData1)), function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorCR.GR.All=c(rep(CorCR.GR.All[1],49))
CorCR.CC.All=round(do.call("rbind",lapply(list(bind_rows(ListData1)), function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorCR.CC.All=c(rep(CorCR.CC.All[1],49))
CorCR.HS.All=round(do.call("rbind",lapply(list(bind_rows(ListData1)), function(x) {cor.test(x[,3],x[,6], method="spearman")[[4]]})),4)
Data$CorCR.HS.All=c(rep(CorCR.HS.All[1],49))

CorCD.GR.All=round(do.call("rbind",lapply(list(bind_rows(ListData2)), function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorCD.GR.All=c(rep(CorCD.GR.All[1],49))
CorCD.CC.All=round(do.call("rbind",lapply(list(bind_rows(ListData2)), function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorCD.CC.All=c(rep(CorCD.CC.All[1],49))
CorCD.HS.All=round(do.call("rbind",lapply(list(bind_rows(ListData2)), function(x) {cor.test(x[,3],x[,6], method="spearman")[[4]]})),4)
Data$CorCD.HS.All=c(rep(CorCD.HS.All[1],49))

CorGR.CC.All=round(do.call("rbind",lapply(list(bind_rows(ListData3)), function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorGR.CC.All=c(rep(CorGR.CC.All[1],49))
CorGR.HS.All=round(do.call("rbind",lapply(list(bind_rows(ListData3)), function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorGR.HS.All=c(rep(CorGR.HS.All[1],49))
CorCC.HS.All=round(do.call("rbind",lapply(list(bind_rows(ListData3)), function(x) {cor.test(x[,4],x[,5], method="spearman")[[4]]})),4)
Data$CorCC.HS.All=c(rep(CorCC.HS.All[1],49))

# Calculate significance values for species
SigCR.GR=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigCR.GR=c(rep(SigCR.GR[1],24),rep(SigCR.GR[2],25))
SigCR.CC=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigCR.CC=c(rep(SigCR.CC[1],24),rep(SigCR.CC[2],25))
SigCR.HS=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,6], method="spearman")[[3]]})),4)
Data$SigCR.HS=c(rep(SigCR.HS[1],24),rep(SigCR.HS[2],25))

SigCD.GR=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigCD.GR=c(rep(SigCD.GR[1],24),rep(SigCD.GR[2],25))
SigCD.CC=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigCD.CC=c(rep(SigCD.CC[1],24),rep(SigCD.CC[2],25))
SigCD.HS=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,6], method="spearman")[[3]]})),4)
Data$SigCD.HS=c(rep(SigCD.HS[1],24),rep(SigCD.HS[2],25))

SigGR.CC=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigGR.CC=c(rep(SigGR.CC[1],24),rep(SigGR.CC[2],25))
SigGR.HS=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigGR.HS=c(rep(SigGR.HS[1],24),rep(SigGR.HS[2],25))
SigCC.HS=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,4],x[,5], method="spearman")[[3]]})),4)
Data$SigCC.HS=c(rep(SigCC.HS[1],24),rep(SigCC.HS[2],25))

# Calculate significance values for genus
SigCR.GR.All=round(do.call("rbind",lapply(list(bind_rows(ListData1)), function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigCR.GR.All=c(rep(SigCR.GR.All[1],49))
SigCR.CC.All=round(do.call("rbind",lapply(list(bind_rows(ListData1)), function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigCR.CC.All=c(rep(SigCR.CC.All[1],49))
SigCR.HS.All=round(do.call("rbind",lapply(list(bind_rows(ListData1)), function(x) {cor.test(x[,3],x[,6], method="spearman")[[3]]})),4)
Data$SigCR.HS.All=c(rep(SigCR.HS.All[1],49))

SigCD.GR.All=round(do.call("rbind",lapply(list(bind_rows(ListData2)), function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigCD.GR.All=c(rep(SigCD.GR.All[1],49))
SigCD.CC.All=round(do.call("rbind",lapply(list(bind_rows(ListData2)), function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigCD.CC.All=c(rep(SigCD.CC.All[1],49))
SigCD.HS.All=round(do.call("rbind",lapply(list(bind_rows(ListData2)), function(x) {cor.test(x[,3],x[,6], method="spearman")[[3]]})),4)
Data$SigCD.HS.All=c(rep(SigCD.HS.All[1],49))

SigGR.CC.All=round(do.call("rbind",lapply(list(bind_rows(ListData3)), function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigGR.CC.All=c(rep(SigGR.CC.All[1],49))
SigGR.HS.All=round(do.call("rbind",lapply(list(bind_rows(ListData3)), function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigGR.HS.All=c(rep(SigGR.HS.All[1],49))
SigCC.HS.All=round(do.call("rbind",lapply(list(bind_rows(ListData3)), function(x) {cor.test(x[,4],x[,5], method="spearman")[[3]]})),4)
Data$SigCC.HS.All=c(rep(SigCC.HS.All[1],49))

# Identify non-significant correlations for species
Data$SigCR.GR=ifelse(Data$SigCR.GR > 0.05, "No", "Yes")
Data$SigCR.CC=ifelse(Data$SigCR.CC > 0.05, "No", "Yes")
Data$SigCR.HS=ifelse(Data$SigCR.HS > 0.05, "No", "Yes")
Data$SigCD.GR=ifelse(Data$SigCD.GR > 0.05, "No", "Yes")
Data$SigCD.CC=ifelse(Data$SigCD.CC > 0.05, "No", "Yes")
Data$SigCD.HS=ifelse(Data$SigCD.HS > 0.05, "No", "Yes")
Data$SigGR.CC=ifelse(Data$SigGR.CC > 0.05, "No", "Yes")
Data$SigGR.HS=ifelse(Data$SigGR.HS > 0.05, "No", "Yes")
Data$SigCC.HS=ifelse(Data$SigCC.HS > 0.05, "No", "Yes")

# Identify non-significant correlations for genus
Data$SigCR.GR.All=ifelse(Data$SigCR.GR.All > 0.05, "No", "Yes")
Data$SigCR.CC.All=ifelse(Data$SigCR.CC.All > 0.05, "No", "Yes")
Data$SigCR.HS.All=ifelse(Data$SigCR.HS.All > 0.05, "No", "Yes")
Data$SigCD.GR.All=ifelse(Data$SigCD.GR.All > 0.05, "No", "Yes")
Data$SigCD.CC.All=ifelse(Data$SigCD.CC.All > 0.05, "No", "Yes")
Data$SigCD.HS.All=ifelse(Data$SigCD.HS.All > 0.05, "No", "Yes")
Data$SigGR.CC.All=ifelse(Data$SigGR.CC.All > 0.05, "No", "Yes")
Data$SigGR.HS.All=ifelse(Data$SigGR.HS.All > 0.05, "No", "Yes")
Data$SigCC.HS.All=ifelse(Data$SigCC.HS.All > 0.05, "No", "Yes")


#################################################
### Plot defense-competitiveness trait spaces ###
#################################################

Plot1=ggplot(Data, aes(GR, CR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(GR, TO.CR.GR, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CRLSD, ymax=CRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=GRLSD, xmax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=1.9-(0.120*0.9), y=0.021-(0.050*0.021), label="NS", color="steelblue2", size=5) +
  geom_text(x=1.9-(0.050*0.9), y=0.021-(0.050*0.021), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=1.0, xmax=1.9, ymin=0.021+(0.021-0.0)*0.010, ymax=0.021+(0.021-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=1.9+(1.9+0.0)*0.010, xmax=1.9+(1.9+0.0)*0.010, ymin=0.0, ymax=0.021) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=1.9+(1.9+1.0)*0.040, xmax=1.9+(1.9+1.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=0.021+(0.021-0.0)*0.040, ymax=0.021+(0.021-0.0)*0.040) + 
  ylab(expression('Maximum clearance rate'~'('~10^2~mL~sec^-1~ind^-1*')')) +
  xlab(expression('Maximum growth rate'~'('*day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(0,0.021)) +
  scale_x_continuous(labels=sprintf(seq(1.0,1.9,by=0.3),fmt="%.1f"), breaks=seq(1.0,1.9,by=0.3), limits=c(1.0,1.9)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Clearance Rate vs Growth Rate.tiff', units="in", width=8, height=8, res=1000)
Plot1
dev.off()

Plot2=ggplot(Data, aes(CC, CR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(CC, TO.CR.CC, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CRLSD, ymax=CRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CCLSD, xmax=CCUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=5.4*10^5-(0.120*5.4*10^5), y=0.021-(0.050*0.021), label="NS", color="steelblue2", size=5) +
  geom_text(x=5.4*10^5-(0.050*5.4*10^5), y=0.021-(0.050*0.021), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=5.4*10^5, ymin=0.021+(0.021-0.0)*0.010, ymax=0.021+(0.021-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=5.4*10^5+(5.4*10^5+0.0)*0.010, xmax=5.4*10^5+(5.4*10^5+0.0)*0.010, ymin=0.0, ymax=0.021) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=5.4*10^5+(5.4*10^5+0.0)*0.040, xmax=5.4*10^5+(5.4*10^5+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=0.021+(0.021-0.0)*0.040, ymax=0.021+(0.021-0.0)*0.040) + 
  ylab(expression('Maximum clearance rate'~'('~10^2~mL~sec^-1~ind^-1*')')) +
  xlab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(0,0.021)) +
  scale_x_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Clearance Rate vs Carrying Capacity.tiff', units="in", width=8, height=8, res=1000)
Plot2
dev.off()

Plot3=ggplot(Data, aes(HS, CR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(HS, TO.CR.HS, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CRLSD, ymax=CRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.120*0.3), y=0.021-(0.050*0.021), label="NS", color="steelblue2", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=0.021-(0.050*0.021), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=0.3, ymin=0.021+(0.021-0.0)*0.010, ymax=0.021+(0.021-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.3+(0.3+0.0)*0.010, xmax=0.3+(0.3+0.0)*0.010, ymin=0.0, ymax=0.021) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=0.3+(0.3+0.0)*0.040, xmax=0.3+(0.3+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=0.021+(0.021-0.0)*0.040, ymax=0.021+(0.021-0.0)*0.040) + 
  ylab(expression('Maximum clearance rate'~'('~10^2~mL~sec^-1~ind^-1*')')) +
  xlab(expression('Half-saturation constant'~'('*µmol~PO[4]^{"-"}~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(0,0.021)) +
  scale_x_continuous(labels=sprintf(seq(0,0.3,by=0.1),fmt="%.1f"), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Clearance Rate vs Half Saturation.tiff', units="in", width=8, height=8, res=1000)
Plot3
dev.off()

Plot4=ggplot(Data, aes(GR, CD, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(GR, TO.CD.GR, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CDLSD, ymax=CDUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=GRLSD, xmax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=1.9-(0.120*0.9), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="steelblue2", size=5) +
  geom_text(x=1.9-(0.050*0.9), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=1.0, xmax=1.9, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=1.9+(1.0+0.0)*0.010, xmax=1.9+(1.0+0)*0.010, ymin=0.7*10^1, ymax=2.8*10^1) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=1.9+(1.0+0.0)*0.040, xmax=1.9+(1.0+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040) + 
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression('Maximum growth rate'~'('*day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0.7,2.8,by=0.7),fmt="%.1f"), breaks=seq(0.7*10^1,2.8*10^1,by=0.7*10^1), limits=c(0.7*10^1,2.8*10^1)) +
  scale_x_continuous(labels=sprintf(seq(1.0,1.9,by=0.3),fmt="%.1f"), breaks=seq(1.0,1.9,by=0.3), limits=c(1.0,1.9)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Cell Area vs Growth Rate.tiff', units="in", width=8, height=8, res=1000)
Plot4
dev.off()

Plot5=ggplot(Data, aes(CC, CD, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(CC, TO.CD.CC, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CDLSD, ymax=CDUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CCLSD, xmax=CCUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=5.4*10^5-(0.120*5.4*10^5), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="steelblue2", size=5) +
  geom_text(x=5.4*10^5-(0.050*5.4*10^5), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=5.4*10^5, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=5.4*10^5+(5.4*10^5+0.0)*0.010, xmax=5.4*10^5+(5.4*10^5+0)*0.010, ymin=0.7*10^1, ymax=2.8*10^1) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=5.4*10^5+(5.4*10^5+0.0)*0.040, xmax=5.4*10^5+(5.4*10^5+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040) + 
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0.7,2.8,by=0.7),fmt="%.1f"), breaks=seq(0.7*10^1,2.8*10^1,by=0.7*10^1), limits=c(0.7*10^1,2.8*10^1)) +
  scale_x_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Cell Area vs Carrying Capacity.tiff', units="in", width=8, height=8, res=1000)
Plot5
dev.off()

Plot6=ggplot(Data, aes(HS, CD, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(HS, TO.CD.HS, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CDLSD, ymax=CDUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.120*0.3), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="steelblue2", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=0.3, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.3+(0.3+0.0)*0.010, xmax=0.3+(0.3+0)*0.010, ymin=0.7*10^1, ymax=2.8*10^1) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=0.3+(0.3+0.0)*0.040, xmax=0.3+(0.3+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040) + 
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression('Half-saturation constant'~'('*µmol~PO[4]^{"-"}~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0.7,2.8,by=0.7),fmt="%.1f"), breaks=seq(0.7*10^1,2.8*10^1,by=0.7*10^1), limits=c(0.7*10^1,2.8*10^1)) +
  scale_x_continuous(labels=sprintf(seq(0,0.3,by=0.1),fmt="%.1f"), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Cell Area vs Half Saturation.tiff', units="in", width=8, height=8, res=1000)
Plot6
dev.off()

# Panel plot of trait spaces
tiff('Defense vs Competitiveness Trait Spaces.tiff', units="in", width=18, height=12, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4,Plot5,Plot6)
grid.arrange(grobs=Panel, ncol=3, nrow=2, layout_matrix=rbind(c(1,1,2,2,3,3),c(4,4,5,5,6,6)))
dev.off()


#########################################################
### Plot competitiveness-competitiveness trait spaces ###
#########################################################

Plot7=ggplot(Data, aes(CC, GR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(CC, TO.GR.CC, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=GRLSD, ymax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CCLSD, xmax=CCUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=5.4*10^5-(0.120*5.4*10^5), y=1.9-(0.050*0.9), label="NS", color="steelblue2", size=5) +
  geom_text(x=5.4*10^5-(0.050*5.4*10^5), y=1.9-(0.050*0.9), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=5.4*10^5, ymin=1.9+(1.0-0.0)*0.010, ymax=1.9+(1.0-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=5.4*10^5+(5.4*10^5+0.0)*0.010, xmax=5.4*10^5+(5.4*10^5+0.0)*0.010, ymin=1.0, ymax=1.9) +
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=5.4*10^5+(5.4*10^5+0.0)*0.040, xmax=5.4*10^5+(5.4*10^5+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=1.9+(1.0-0.0)*0.040, ymax=1.9+(1.0-0.0)*0.040) + 
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(1.0,1.9,by=0.3),fmt="%.1f"), breaks=seq(1.0,1.9,by=0.3), limits=c(1.0,1.9)) +
  scale_x_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Growth Rate vs Carrying Capacity.tiff', units="in", width=8, height=8, res=1000)
Plot7
dev.off()

Plot8=ggplot(Data, aes(HS, GR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(HS, TO.GR.HS, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=GRLSD, ymax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.120*0.3), y=1.9-(0.050*0.9), label="NS", color="steelblue2", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=1.9-(0.050*0.9), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=0.3, ymin=1.9+(1.0-0.0)*0.010, ymax=1.9+(1.0-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.3+(0.3+0.0)*0.010, xmax=0.3+(0.3+0.0)*0.010, ymin=1.0, ymax=1.9) +
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=0.3+(0.3+0.0)*0.040, xmax=0.3+(0.3+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=1.9+(1.0-0.0)*0.040, ymax=1.9+(1.0-0.0)*0.040) + 
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression('Half-saturation constant'~'('*µmol~PO[4]^{"-"}~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(1.0,1.9,by=0.3),fmt="%.1f"), breaks=seq(1.0,1.9,by=0.3), limits=c(1.0,1.9)) +
  scale_x_continuous(labels=sprintf(seq(0,0.3,by=0.1),fmt="%.1f"), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Growth Rate vs Half Saturation.tiff', units="in", width=8, height=8, res=1000)
Plot8
dev.off()

Plot9=ggplot(Data, aes(HS, CC, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(HS, TO.CC.HS, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CCLSD, ymax=CCUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.120*0.3), y=5.4*10^5-(0.050*5.4*10^5), label="NS", color="steelblue2", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=5.4*10^5-(0.050*5.4*10^5), label="NS", color="tomato3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=0.3, ymin=5.4*10^5+(5.4*10^5-0.0)*0.010, ymax=5.4*10^5+(5.4*10^5-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.3+(0.3+0.0)*0.010, xmax=0.3+(0.3+0.0)*0.010, ymin=0.0, ymax=5.4*10^5) +
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=90), xmin=0.3+(0.3+0.0)*0.040, xmax=0.3+(0.3+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col="black", fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=5.4*10^5+(5.4*10^5-0.0)*0.040, ymax=5.4*10^5+(5.4*10^5-0.0)*0.040) + 
  ylab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Half-saturation constant'~'('*µmol~PO[4]^{"-"}~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5)) +
  scale_x_continuous(labels=sprintf(seq(0,0.3,by=0.1),fmt="%.1f"), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="steelblue2", "Chlamydomonas bilatus"="tomato3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Carrying Capacity vs Half Saturation.tiff', units="in", width=8, height=8, res=1000)
Plot9
dev.off()

# Panel plot of trait spaces
tiff('Competitiveness vs Competitiveness Trait Spaces.tiff', units="in", width=18, height=6, res=1000)
Panel=list(Plot7,Plot8,Plot9)
grid.arrange(grobs=Panel, ncol=3, nrow=1, layout_matrix=rbind(c(1,1,2,2,3,3)))
dev.off()

# Panel plot of trait spaces
tiff('Trait Spaces.tiff', units="in", width=18, height=18, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4,Plot5,Plot6,Plot7,Plot8,Plot9)
grid.arrange(grobs=Panel, ncol=3, nrow=3, layout_matrix=rbind(c(1,1,2,2,3,3),c(4,4,5,5,6,6),c(7,7,8,8,9,9)))
dev.off()


##################
### Statistics ###
##################

# Fixed effect regressions
Model1=lmer(GR~1 +(1|Species), REML=T, data=Data)
Model2=lmer(CC~1 +(1|Species), REML=T, data=Data)
Model3=lmer(HS~1 +(1|Species), REML=T, data=Data)
Model4=lmer(CR~1 +(1|Species), REML=T, data=Data)
Model5=lmer(CD~1 +(1|Species), REML=T, data=Data)

# Calculate the interspecific vs intraspecific contributions
summary(Model1)
InterGR=0.0599/(0.1273+0.0599)*100
IntraGR=0.1273/(0.1273+0.0599)*100

summary(Model2)
InterCC=0.6350/(0.7601+0.6350)*100
IntraCC=0.7601/(0.7601+0.6350)*100

summary(Model3)
InterHS=0.1204/(0.0001+0.1204)*100
IntraHS=0.0001/(0.0001+0.1204)*100

summary(Model4)
InterCR=0.0008/(0.0025+0.0008)*100
IntraCR=0.0025/(0.0025+0.0008)*100

summary(Model5)
InterCA=0.7969/(1.5964+0.7969)*100
IntraCA=1.5964/(1.5964+0.7969)*100
