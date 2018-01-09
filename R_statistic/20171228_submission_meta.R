###### Import data sheets ------------------------------------------------------
# This version is close to the final version, however, the dataset clean querteria is reads has to be greater than >0.001% is not sufficient for this dataset.
# so we employed reads <10.
# OTU_ is for data analyse

#set the work routine
setwd("/Users/JiayiQin/Dropbox/2nd_Miseq")

## Pretreat working excel #####
OTU<-read.csv("Enchy_262_1f_OTU.filtered.table.csv",header=TRUE,sep=",")
Local_blast<- read.csv("/Users/JiayiQin/Dropbox/2nd_Miseq/RECHECK/Enchy_262_1f_OTU.filtered.local.10output.262length.filtered.top.blast",header=TRUE,sep="\t")


Ench<-merge(OTU,Local_blast, by="amplicon",all.y=TRUE,all.x = TRUE)

Ench$Species<-as.character(Ench$Species)
for ( i in 1:nrow(Ench)) 
{
  if (is.na(Ench$Species[i]))
  {Ench$Species[i]="Unannoated"}
}


## Organise dataset ######

Ench_c_t <- cbind(Ench$Species,Ench[,14:45])
colnames(Ench_c_t)<-c("Taxa",colnames(Ench)[14:45])
Ench_c_t<-Ench_c_t[,-8]

Ged<-aggregate(Ench_c_t[,2:32], by=list(Taxa=Ench_c_t$Taxa), FUN=sum)

Ged.colnames<-substr(colnames(Ench_c_t[,2:32]),1,unlist(gregexpr("_",colnames(Ench_c_t[,2:32])))-1)
colnames(Ged)<-c("Taxa",Ged.colnames)

Ged.t<-t(Ged[,2:32])
Ged.t<-cbind(Ged.t,row.names(Ged.t))
colnames(Ged.t)<-c(as.character(Ged$Taxa),"Sample")
  
# Treatment setting
Treat<-read.csv("Treat.csv",header=TRUE)

Ged_Treatment<-merge(Ged.t,Treat,by="Sample")

# Define function to change the variable type
IntToNum <- function(x) {if(is.integer(x)) 
{as.numeric(paste(x))} else {x}} 
FacToNum <- function(x) {if(is.factor(x)) 
{as.numeric(paste(x))} else {x}} 

# Change variable property in OTU table
for (i in 2:7){
  Ged_Treatment[,i]<-FacToNum(Ged_Treatment[,i])}

# Get rid of the irreleastic reads, rltv abundance < 0.001%
for ( i in 1:31){
  for (j in 2:6){
    if (Ged_Treatment[i,j]/sum(Ged_Treatment[i,2:7])<0.00001){
      Ged_Treatment[i,j]<-0
    }
  }
}


t_Ged<-t(Ged[,2:32])
t_Ged<-as.data.frame( t_Ged )

Taxa<-Ged[,1]

Ged_10<-Ged_Treatment[,1:7]
t_Ged_10<-t(Ged_10[,2:7])
colnames(t_Ged_10)<-Ged_10[,1]

Ged_10_rafy<-cbind(row.names(t_Ged_10),t_Ged_10)
library(vegan)
Ged_10_rafy[,2:32]<-decostand(t_Ged_10,"pa",2)

# remove un assigned data

Ged_10_rafy<-Ged_10_rafy[-6,]

########## No. of species #####################
Ged_10_rafy<-as.data.frame(Ged_10_rafy)
for (i in 2:32){
  Ged_10_rafy[,i]<-FacToNum(Ged_10_rafy[,i])}

SpNum<-colSums(Ged_10_rafy[,2:32])

SpNum<-as.data.frame(SpNum)
SpNum[,2]<-row.names(SpNum)
colnames(SpNum)<-c("SpNum","Sample")
Ged_Treatment_sum<-merge(Ged_Treatment,SpNum,by.Ged_Treatment="Sample",by.SpNum="Sample")
# SpNum is not normal distributed, but the residual of fitted modle is normal distributed, which reaches the requirement of lmer function
shapiro.test(Ged_Treatment_sum$SpNum)
hist(Ged_Treatment_sum$SpNum)

kruskal.test(SpNum~Ged_Treatment_sum$Treat,data=Ged_Treatment_sum) #fast test to look at the Treatment effect

Ged_Treatment_sum$Treat<-factor(Ged_Treatment_sum$Treat)
SpNum_lm<-lmer(SpNum~Treat+(1|Block/Subsample), data=Ged_Treatment_sum)
SpNum_lm_2<-lmer(SpNum~(1|Block/Subsample), data=Ged_Treatment_sum)
anova(SpNum_lm,SpNum_lm_2)
shapiro.test(residuals(SpNum_lm))

########## shanon_alpha diversity ##########
shannon<-diversity(Ged_10[,2:6], index = "shannon", MARGIN = 1, base = exp(1))

# Organise the shannon output, associate the sample name to the result
shannon<-as.data.frame(shannon)
shannon$Sample<-Ged_10$Sample
# Merge the Shannon output to the big table
Ged_Treatment_sum<-merge(Ged_Treatment_sum,shannon,by.Ged_Treatment_sum="Sample",by.shannon="Sample")


kruskal.test(shannon~Ged_Treatment_sum$Treat,data=Ged_Treatment_sum)

shannon_lm<-lmer(shannon~Treat+(1|Block/Subsample), data=Ged_Treatment_sum)
shannon_lm_2<-lmer(shannon~(1|Block/Subsample), data=Ged_Treatment_sum)
anova(shannon_lm,shannon_lm_2)
residuals(shannon_lm)
shapiro.test(residuals(shannon_lm))

################################################################
fisher<-fisher.alpha(Ged_10[,2:6], MARGIN = 1, se = FALSE)
fisher<-as.data.frame(fisher)
fisher$Sample<-Ged_10$Sample
Ged_Treatment_sum<-merge(Ged_Treatment_sum,fisher,by.Ged_Treatment_sum="Sample",by.fisher="Sample")

kruskal.test(fisher~Ged_Treatment_sum$Treat,data=Ged_Treatment_sum)

fisher_lm<-lmer(fisher~Treat+(1|Block/Subsample), data=Ged_Treatment_sum)
fisher_lm_2<-lmer(fisher~(1|Block/Subsample), data=Ged_Treatment_sum)
anova(fisher_lm,fisher_lm_2)
residuals(fisher_lm)
shapiro.test(residuals(fisher_lm))
################################################################

########### Pielou’s evenness ########## 
# Pielou's evenness: J = H0/ log(S)
J <- shannon$shannon/log(specnumber(Ged_10[,2:6]))

J<-as.data.frame(J)
J$Sample<-Ged_10$Sample
Ged_Treatment_sum<-merge(Ged_Treatment_sum,J,by.Ged_Treatment_sum="Sample",by.J="Sample")


kruskal.test(J~Treat,data=Ged_Treatment_sum)
J_lm<-lmer(J~Treat+(1|Block/Subsample), data=Ged_Treatment_sum)
J_lm_2<-lmer(J~(1|Block/Subsample), data=Ged_Treatment_sum)
anova(J_lm,J_lm_2)
residuals(J_lm)
shapiro.test(residuals(J_lm))


######### PERMANOVA incidence dataset ###################
t_Ged_10_rafy <- t(Ged_10_rafy[,2:32])
t_Ged_10_rafy<-as.data.frame(t_Ged_10_rafy)
t_Ged_10_rafy$Sample<-row.names(t_Ged_10_rafy)
t_Ged_10_rafy<-merge(t_Ged_10_rafy,Ged_Treatment_sum,by="Sample")

betad<-betadiver(t_Ged_10_rafy[,2:6] ,"z")

adonis(betad~factor(t_Ged_10_rafy$Treat),method="Jaccard",
       strata=t_Ged_10_rafy$Block, 
       perm=1000) # p=0.968 with Block as Strata, p=0.971 without strata


mod<-betadisper(betad,t_Ged_10_rafy$Treat)
# calculate the percentage, PCoA plot could explain in each axis
eig<-eigenvals(mod)
axes_percent<-eig/sum(eig)
axes_percent # display the result
plot(mod,
     col= c("blue","green","red","orange"),
     main="",
     pch = c(15,22,17,24),
     xlab="",
     ylab="",
     ylim=c(-0.2,0.25),
     xlim=c(-0.2,0.4),
     sub="",
     segments = FALSE, # control grey line
     hull = FALSE, ellipse = FALSE,label = F,lwd=0.5)
mtext("PCoA1 - Percent variance explained 79.69%", side = 1, line = 2)
mtext("PCoA2 - Percent variance explained 32.87%", side = 2, line = 2)
mtext("PCoA - PC1 vs PC2", side = 3, line = 1)
ordiellipse(mod, t_Ged_10_rafy$Treat, 
            kind="se", 
            conf=0.95, 
            lwd=1, 
            col= c("blue","green","red","orange")
)
legend("bottomleft", c("0 t/ha", "3 t/ha", "4.5 t/ha", "6 t/ha"), 
       col = c("blue","green","red","orange"),
       title="incidence dataset",
       text.col = "black", 
       # lty = c(1, 1, 1,1),
       pch = c(15,22,17,24),
       bty="n")

boxplot(mod)
anova(mod)

plot(mod)
boxplot(mod)

######### PERMANOVA relative abundance dataset ###################
Ged_10_rltv <- decostand(t(Ged_10[,2:7]),"total",2)
t_Ged_10_rltv <- t(Ged_10_rltv)
t_Ged_10_rltv<-as.data.frame(t_Ged_10_rltv)
t_Ged_10_rltv$Sample<-Ged_10$Sample
t_Ged_10_rltv<-merge(t_Ged_10_rltv,Ged_Treatment_sum,by="Sample")

betad<-betadiver(t_Ged_10_rltv[,2:6] ,"z")
adonis(betad~factor(t_Ged_10_rltv$Treat),method="bray",
       strata=t_Ged_10_rltv$Block, 
       perm=1000) 

mod<-betadisper(betad,t_Ged_10_rltv$Treat)
plot(mod,
     col= c("blue","green","red","orange"),
     main="Multivariate Dispersion",
     pch = c(15,22,17,24),
     segments = FALSE, # control grey line
     hull = FALSE, ellipse = FALSE,label = F,lwd=0.5)
ordiellipse(mod, t_Ged_10_rltv$Treat, 
            kind="se", 
            conf=0.95, 
            lwd=1, 
            col= c("blue","green","red","orange")
)
legend("bottomleft", c("0 t/ha", "3 t/ha", "4.5 t/ha", "6 t/ha"), 
       col = c("blue","green","red","orange"),
       title="relative abundance dataset",
       text.col = "black", 
       # lty = c(1, 1, 1,1),
       pch = c(15,22,17,24),
       bty="n")

#### Ratio ##########################
#no significant difference
V2_V1<-Ged_10[,3]/Ged_10[,2]
V4_V1<-Ged_10[,5]/Ged_10[,2] # Ch. sphagnetorum and C. chlorophilus 
V5_V1<-Ged_10[,6]/Ged_10[,2] # M.  pelicensis and C. chlorophilus
V3_V1<-Ged_10[,4]/Ged_10[,2]

ratio<-rbind(V2_V1,V3_V1,V4_V1,V5_V1)
ratio<-t(ratio)
ratio<-cbind(Ged_Treatment[,9:11],ratio)


ratio$Treat<-factor(ratio$Treat)
#V2_V1
shapiro.test(sqrt(ratio$V2_V1))
aov_V2_V1<- lmer(ratio$V2_V1~ratio$Treat+(1|Block/Subsample),data=ratio,REML=FALSE) 
shapiro.test(residuals(aov_V2_V1)) # The residuals are not normal distributed

summary(aov_V2_V1)
fixef(aov_V2_V1)
vcov(aov_V2_V1)
plot(aov_V2_V1)

aov_V2_V1_2<- lmer(ratio$V2_V1~(1|Block/Subsample),data=ratio,REML=FALSE)
shapiro.test(residuals(aov_V2_V1_2))
anova(aov_V2_V1,aov_V2_V1_2)

lsmeans(aov_V2_V1,pairwise~Treat,data=ratio,adjust="tukey")

# V4_V1 # the residuals are not normal distributed
aov_V4_V1<- lmer(V4_V1~Treat+(1|Block/Subsample),data=ratio,REML=FALSE)
aov_V4_V1_2<- lmer(V4_V1~(1|Block/Subsample),data=ratio,REML=FALSE)
shapiro.test(residuals(aov_V4_V1_2))
anova(aov_V4_V1,aov_V4_V1_2)

summary(aov_V4_V1)
plot(aov_V4_V1)
lsmeans(aov_V4_V1,pairwise~Treat,data=ratio,adjust="tukey")

# V5_V1 not normal distributed yet
aov_V5_V1<- lmer(V5_V1~Treat+(1|Block/Subsample),data=ratio,REML=FALSE)
shapiro.test(residuals(aov_V5_V1))
aov_V5_V1_2<- lmer(V5_V1~(1|Block/Subsample),data=ratio,REML=FALSE)
shapiro.test(residuals(aov_V5_V1_2))
anova(aov_V5_V1,aov_V5_V1_2)

summary(aov_V5_V1)
plot(aov_V5_V1)
lsmeans(aov_V5_V1,pairwise~Treat,data=ratio,adjust="tukey")

# V3_V1 not normal distributed yet
aov_V3_V1<- lmer(V3_V1~Treat+(1|Block/Subsample),data=ratio,REML=FALSE)
shapiro.test(residuals(aov_V3_V1))
aov_V3_V1_2<- lmer(V3_V1~(1|Block/Subsample),data=ratio,REML=FALSE)
shapiro.test(residuals(aov_V3_V1_2))
anova(aov_V3_V1,aov_V3_V1_2)

summary(aov_V3_V1)
plot(aov_V3_V1)
lsmeans(aov_V3_V1,pairwise~Treat,data=ratio,adjust="tukey")


aov_V4_V1<-  kruskal.test(V4_V1~factor(ratio$Treat),data=ratio)
summary(aov_V4_V1)
aov_V5_V1<-  kruskal.test(V5_V1~factor(Treat),data=ratio)
summary(aov_V5_V1)
aov_V3_V1<-  kruskal.test(V3_V1~factor(Treat),data=ratio)
summary(aov_V3_V1)

# Visualize the ratio with error bars
library(Rmisc)
V2_V1_sum<- summarySE(data = ratio, 
                      measurevar = "V2_V1", 
                      groupvars = c("Treat"))
V4_V1_sum<- summarySE(data = ratio, 
                      measurevar = "V4_V1", 
                      groupvars = c("Treat"))
V3_V1_sum<- summarySE(data = ratio, 
                      measurevar = "V3_V1", 
                      groupvars = c("Treat"))
V5_V1_sum<- summarySE(data = ratio, 
                      measurevar = "V5_V1", 
                      groupvars = c("Treat"))
library(gplots)
plotCI(x = V2_V1_sum$V2_V1, 
       uiw = V2_V1_sum$se, 
       xaxt ="n", 
       las = 1,
       #xlim = c(0.5,2.5), 
       #ylim = c(0,100), 
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       #type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "blue")

plotCI(x = V4_V1_sum$V4_V1, 
       uiw = V4_V1_sum$se, 
       xaxt ="n", 
       las = 1,
       #xlim = c(0.5,2.5), 
       #ylim = c(0,100), 
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       #type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "blue")
plotCI(x = V3_V1_sum$V3_V1, 
       uiw = V3_V1_sum$se, 
       xaxt ="n", 
       las = 1,
       #xlim = c(0.5,2.5), 
       #ylim = c(0,100), 
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "blue")
plotCI(x = V5_V1_sum$V5_V1, 
       uiw = V5_V1_sum$se, 
       xaxt ="n", 
       las = 1,
       #xlim = c(0.5,2.5), 
       #ylim = c(0,100), 
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       #type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "blue")







library("vegan")


############### use dominante species as indicator ##########
Chamaedrilus_chlorophilus<-t_Ged_10_rltv[,2]
x<-log(sqrt(Chamaedrilus_chlorophilus)+0.3)
NormalDistribution(x)
hist(x)
aov_Chamaedrilus_chlorophilus<-aov(log(sqrt(Chamaedrilus_chlorophilus)+0.3)~factor(Ged_Treatment_sum$Treat))
summary(aov_Chamaedrilus_chlorophilus)
anova(aov_Chamaedrilus_chlorophilus)

plot(Chamaedrilus_chlorophilus~factor(Ged_Treatment_sum$Treat))



### Barplot of OTUs ---------------------------------------------------------
Ged_barplot<-cbind(Ged_Treatment_sum$Treat,Ged_Treatment_sum[,2:7])

Ged_barplot<-t(Ged_barplot)
Ged_barplot<-as.data.frame(Ged_barplot)
names(Ged_barplot)<-Ged_barplot[1,]
Ged_barplot<-Ged_barplot[-1,]
for (i in 1:31){
  Ged_barplot[,i]<-FacToNum( Ged_barplot[,i])
}
Ged_barplot<-decostand(Ged_barplot,"total",2)

library(reshape2)
Ged_barplot$row <- row.names(Ged_barplot)

Ged_species_SE<-t(Ged_barplot[1:31])
Ged_species_SE<-as.data.frame(Ged_species_SE)
Ged_species_SE$Treat<-Ged_Treatment_sum$Treat

Chamaedrilus_chlorophilus_se<-summarySE(data = Ged_species_SE, 
                                        measurevar = "Chamaedrilus chlorophilus", 
                                        groupvars = c("Treat"))
Chamaedrilus_sphagnetorum_se<-summarySE(data = Ged_species_SE, 
                                            measurevar = "Chamaedrilus sphagnetorum", 
                                            groupvars = c("Treat"))
Chamaedrilus_glandulosus_se<-summarySE(data = Ged_species_SE, 
                                              measurevar = "Chamaedrilus glandulosus", 
                                              groupvars = c("Treat"))
Chamaedrilus_lapponicus_se<-summarySE(data = Ged_species_SE, 
                                             measurevar = "Chamaedrilus lapponicus", 
                                             groupvars = c("Treat"))
Mesenchytraeus_pelicensis_se<-summarySE(data = Ged_species_SE, 
                                        measurevar = "Mesenchytraeus pelicensis", 
                                        groupvars = c("Treat"))

Ged_barplot_2<-aggregate(Ged_species_SE[,1:6],by=list(Ged_species_SE$Treat),mean) # calculate average

colnames(Ged_barplot_2)[1]<-"Treat"
Ged_barplot_2 <- melt(Ged_barplot_2, id.vars = "Treat") # reorganise the table

library(ggplot2)
require(scales)
library(plyr)
### without break #####
# https://groups.google.com/forum/#!topic/ggplot2/jSrL_FnS8kc
# It is impossible to produce stacked barplot with break in y-axis
charts.data <- ddply(Ged_barplot_2, .(variable),
                     transform, pos = cumsum(value) - (0.5 * value)) #calculate the error bar

Ged_barplot_2$Treat<-factor(Ged_barplot_2$Treat)

require("RColorBrewer")
fill<-c("#B2DF8A","#33A02C","#1F78B4","#A6CEE3","#E31A1C","#FDBF6F")

stackbarplot<- ggplot(Ged_barplot_2, aes(x=Treat, y=value,fill=variable) )+ 
  geom_bar(stat="identity") +
  xlab("\nTreatment") +
  ylab("Percent\n") +
  scale_y_continuous(labels = percent_format(),
                     #limits=c(0.7,1),
                     oob = rescale_none) +
  scale_fill_manual(values=fill) +
  theme_bw()

stackbarplot

stackbarplot+geom_text(data=charts.data, aes(x = variable, y = value, label = paste0(value,"%")),
                       size=4)

