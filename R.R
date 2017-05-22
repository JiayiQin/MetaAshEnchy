###### Import data sheets ------------------------------------------------------
# 
# OTU_ is for data analyse

#set the work routine
setwd("/Users/JiayiQin/Dropbox/2nd_Miseq")

## Pretreat working excel #####
OTU<-read.csv("Enchy_262_1f_OTU.filtered.table.csv",header=TRUE,sep=",")
Local_blast<-read.csv("Enchy_262_1f_OTU.filtered.local.blast.csv",header=TRUE)
Blast<-read.csv("Enchy_262_1f_OTU.filtered.blast.csv",header=TRUE)

Ench<-merge(OTU,Local_blast, by.OTU="amplicon", by.Local_blast="amplicon")
names(Blast)<-c("amplicon",1,2,3,4,5,6,7,8,9,10,11,12,13)
Ench<-merge(Ench,Blast, by.Ench="amplicon", by.Blast="amplicon")

write.csv(Ench,"Enchy_merge.csv")

Ench_c <- read.csv("Enchy_merge_correction_2.csv",header=TRUE)

## Import dataset ######
Ench_c_t <- read.csv("Enchy_merge_correction_tail.csv",header=TRUE)

Ged<-aggregate(Ench_c_t[,2:32], by=list(Taxa=Ench_c_t$Taxa), FUN=sum)
write.csv(Ged,"Enchyreads_to_present.csv")

# Treatment setting
Treat<-read.csv("Treat.csv",header=TRUE)
#
names(Ged)<-c("Taxa",Treat$Treat)
IntToNum <- function(x) {if(is.integer(x)) 
{as.numeric(paste(x))} else {x}} 

for (i in 2:32){
Ged[,i]<-IntToNum(Ged[,i])}

t_Ged<-t(Ged[,2:32])
t_Ged<-as.data.frame( t_Ged )


Taxa<-Ged[,1]

#### ***USEFUL: delete reads lower than 1o and rafy   ######
Ged_10<-Ged

for ( i in 1:8){
  for (j in 2:32){
    if (Ged[i,j]>10 )
    {Ged_10[i,j]<-Ged[i,j]}
    else
    {Ged_10[i,j]<-0}
  }
}

Ged_10_rafy<-Ged_10
Ged_10_rafy[,2:32]<-decostand(Ged_10[2:32],"pa",2)
# remove un assigned data
 # Ged_10_rafy<-Ged_10_rafy[-5,]

# No. of species
SpNum<-colSums(Ged_10_rafy[,2:32])
 # SpNum is not normal distributed, but the residual of fitted modle is normal distributed, which reaches the requirement of lmer function
shapiro.test(log(sqrt(SpNum)))
hist(log(sqrt(SpNum)))

kruskal.test(SpNum~Treat$Treat)

Treat$Treat<-factor(Treat$Treat)
SpNum_lm<-lmer(SpNum~Treat+(1|Block/Subsample), data=Treat)
SpNum_lm_2<-lmer(SpNum~(1|Block/Subsample), data=Treat)
anova(SpNum_lm,SpNum_lm_2)
shapiro.test(residuals(SpNum_lm))

# shanon_alpha diversity
shannon<-diversity(t(Ged_10[-5,2:32]), index = "shannon", MARGIN = 1, base = exp(1))
kruskal.test(shanon~Treat$Treat)

shannon_lm<-lmer(shannon~Treat+(1|Block/Subsample), data=Treat)
shannon_lm_2<-lmer(shannon~(1|Block/Subsample), data=Treat)
anova(shannon_lm,shannon_lm_2)
residuals(shannon_lm)
shapiro.test(residuals(shannon_lm))

fisher<-fisher.alpha(t(Ged_10[2:32]), MARGIN = 1, se = FALSE)
kruskal.test(fisher~Treat$Treat)

# Pielou’s evenness J = H0/ log(S) is easily found as:
J <- shannon/log(specnumber(t(Ged_10[-5,2:32])))
kruskal.test(J~Treat$Treat)
J_lm<-lmer(J~Treat+(1|Block/Subsample), data=Treat)
J_lm_2<-lmer(J~(1|Block/Subsample), data=Treat)
anova(J_lm,J_lm_2)
residuals(J_lm)
shapiro.test(residuals(J_lm))


t_Ged_10_rafy <- t(Ged_10_rafy[,2:32])
betad<-betadiver(t_Ged_10_rafy ,"z")

adonis(betad~factor(Treat$Treat),method="bray",
       strata=Treat$Block, 
       perm=1000) 


mod<-betadisper(betad,Treat$Treat)
eig<-eigenvals(mod)
axes_percent<-eig/sum(eig)
axes_percent
plot(mod,
     col= c("blue","green","red","orange"),
     main="",
     pch = c(15,22,17,24),
     xlab="",
     ylab="",
     sub="",
     segments = FALSE, # control grey line
     hull = FALSE, ellipse = FALSE,label = F,lwd=0.5)
mtext("PC1 - Percent variance explained 63.21%", side = 1, line = 2)
mtext("PC2 - Percent variance explained 29.90%", side = 2, line = 2)
mtext("PCoA - PC1 vs PC2", side = 3, line = 1)
ordiellipse(mod, Treat$Treat, 
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

Ged_10_rltv <- decostand(Ged_10[,2:32],"total",2)
t_Ged_10_rltv <- t(Ged_10_rltv)

betad<-betadiver(t_Ged_10_rltv ,"z")
adonis(betad~factor(Treat$Treat),method="bray",
       strata=Treat$Block, 
       perm=1000) 

mod<-betadisper(betad,Treat$Treat)
plot(mod,
     col= c("blue","green","red","orange"),
     main="Multivariate Dispersion",
     pch = c(15,22,17,24),
     segments = FALSE, # control grey line
     hull = FALSE, ellipse = FALSE,label = F,lwd=0.5)
ordiellipse(mod, Treat$Treat, 
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
#### ***USEFUL: Ratio ##########################
#no significant difference
V1_V3<-Ged[1,2:32]/Ged[3,2:32]
V2_V3<-Ged[2,2:32]/Ged[3,2:32]
V4_V3<-Ged[4,2:32]/Ged[3,2:32]
V6_V3<-Ged[6,2:32]/Ged[3,2:32]
V7_V3<-Ged[7,2:32]/Ged[3,2:32]
V2_V7<-Ged[2,2:32]/Ged[7,2:32]

ratio<-rbind(V1_V3,V2_V3,V4_V3,V6_V3,V7_V3,V2_V7)
ratio<-t(ratio)
ratio<-cbind(Treat,ratio)
names(ratio)<-c("Sample","No.","Block","Treat","Subsample","V1_V3","V2_V3","V4_V3","V6_V3","V7_V3","V2_V7")

ratio$Treat<-factor(ratio$Treat)
#V1_V3
NormalDistribution(sqrt(ratio$V1_V3))
aov_V1_V3<- lmer(sqrt(ratio$V1_V3)~Treat+(1|Block/Subsample),data=ratio,REML=FALSE)
NormalDistribution(residuals(aov_V1_V3))

summary(aov_V1_V3)
fixef(aov_V1_V3)
vcov(aov_V1_V3)
plot(aov_V1_V3)

aov_V1_V3_2<- lmer(sqrt(ratio$V1_V3)~(1|Block/Subsample),data=ratio,REML=FALSE)
NormalDistribution(residuals(aov_V1_V3_2))
anova(aov_V1_V3,aov_V1_V3_2)

lsmeans(aov_V1_V3,pairwise~Treat,data=ratio,adjust="tukey")

# V2_V3 # the residuals are not normal distributed
aov_V2_V3<- lmer(V2_V3~Treat+(1|Block/Subsample),data=ratio,REML=FALSE)
NormalDistribution(residuals(aov_V2_V3))
aov_V2_V3_2<- lmer(V2_V3~(1|Block/Subsample),data=ratio,REML=FALSE)
shapiro.test(residuals(aov_V2_V3_2))
anova(aov_V2_V3,aov_V2_V3_2)

summary(aov_V2_V3)
plot(aov_V2_V3)
lsmeans(aov_V2_V3,pairwise~Treat,data=ratio,adjust="tukey")

# V4_V3 not normal distributed yet
aov_V4_V3<- lmer(V4_V3~Treat+(1|Block/Subsample),data=ratio,REML=FALSE)
NormalDistribution(residuals(aov_V4_V3))
aov_V4_V3_2<- lmer(V4_V3~(1|Block/Subsample),data=ratio,REML=FALSE)
NormalDistribution(residuals(aov_V4_V3_2))
anova(aov_V4_V3,aov_V4_V3_2)

summary(aov_V4_V3)
plot(aov_V4_V3)
lsmeans(aov_V4_V3,pairwise~Treat,data=ratio,adjust="tukey")

# V6_V3 not normal distributed yet
aov_V6_V3<- lmer(V6_V3~Treat+(1|Block/Subsample),data=ratio,REML=FALSE)
NormalDistribution(residuals(aov_V6_V3))
aov_V6_V3_2<- lmer(V6_V3~(1|Block/Subsample),data=ratio,REML=FALSE)
NormalDistribution(residuals(aov_V6_V3_2))
anova(aov_V6_V3,aov_V6_V3_2)

summary(aov_V6_V3)
plot(aov_V6_V3)
lsmeans(aov_V6_V3,pairwise~Treat,data=ratio,adjust="tukey")

# V7_V3 not normal distributed yet
aov_V7_V3<- lmer(V7_V3~Treat+(1|Block/Subsample),data=ratio,REML=FALSE)
NormalDistribution(residuals(aov_V7_V3))
aov_V7_V3_2<- lmer(V7_V3~(1|Block/Subsample),data=ratio,REML=FALSE)
shapiro.test(residuals(aov_V7_V3_2))
anova(aov_V7_V3,aov_V7_V3_2)
summary(aov_V7_V3)
plot(aov_V7_V3)
lsmeans(aov_V7_V3,pairwise~Treat,data=ratio,adjust="tukey")


aov_V2_V3<-  kruskal.test(V2_V3~factor(Treat),data=ratio)
summary(aov_V2_V3)
aov_V4_V3<-  kruskal.test(V4_V3~factor(Treat),data=ratio)
summary(aov_V4_V3)
aov_V6_V3<-  kruskal.test(V6_V3~factor(Treat),data=ratio)
summary(aov_V6_V3)
aov_V7_V3<-  kruskal.test(V7_V3~factor(Treat),data=ratio)
summary(aov_V7_V3)

V1_V3_sum<- summarySE(data = ratio, 
            measurevar = "V1_V3", 
            groupvars = c("Treat"))
V2_V3_sum<- summarySE(data = ratio, 
                      measurevar = "V2_V3", 
                      groupvars = c("Treat"))
V6_V3_sum<- summarySE(data = ratio, 
                      measurevar = "V6_V3", 
                      groupvars = c("Treat"))
V7_V3_sum<- summarySE(data = ratio, 
                      measurevar = "V7_V3", 
                      groupvars = c("Treat"))
library(gplots)
plotCI(x = V1_V3_sum$V1_V3, 
       uiw = V1_V3_sum$se, 
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

plotCI(x = V2_V3_sum$V2_V3, 
       uiw = V2_V3_sum$se, 
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
plotCI(x = V6_V3_sum$V6_V3, 
       uiw = V6_V3_sum$se, 
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
plotCI(x = V7_V3_sum$V8_V3, 
       uiw = V8_V3_sum$se, 
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

################### rltv ######################
Ged_rltv<-Ged
Ged_rltv[,2:32]<-decostand(Ged[,2:32],"total",2)

t_Ged_rltv <- t(Ged_rltv[,2:32])
t_Ged_rltv<-t_Ged_rltv[sort(row.names(t_Ged_rltv)),]

betad<-betadiver(t_Ged_rltv ,"z")

TREAT<-sort(Treat$Treat)
adonis(betad~Treat$Treat+Treat$Block,method="bray",
       #strata=Extract, 
       perm=1000) 

mod<-betadisper(betad,Treat$Treat)

plot(mod)

boxplot(mod)

############# rafy 
Ged_rafy<-Ged
Ged_rafy[,2:32]<-decostand(Ged[,2:32],"pa",2)
t_Ged_rafy <- t(Ged_rafy[,2:32])

betad<-betadiver(t_Ged_rafy ,"z")

TREAT<-sort(Treat$Treat)
adonis(betad~Treat$Treat,method="bray",
       strata=Treat$Block, 
       perm=1000) 
mod<-betadisper(betad,Treat$Treat)

plot(mod)

boxplot(mod)
 


############### use dominante species as indicator ##########
Cognettia_sphagnetorum<-t_Ged_rltv[,2]
x<-log(sqrt(Cognettia_sphagnetorum)+0.3)
NormalDistribution(x)
hist(x)
aov_Cognettia_sphagnetorum<-aov(log(sqrt(Cognettia_sphagnetorum)+0.3)~factor(Treat$Treat))
summary(aov_Cognettia_sphagnetorum)
anova(aov_Cognettia_sphagnetorum)

plot(Cognettia_sphagnetorum~factor(Treat$Treat))

############# CAT_1234 ############
Ged_sort_rltv <-decostand(Ged[2:32],"total",2)
Ged_sort_rltv_cat<-Ged_sort_rltv

for ( i in 1:8){
  for (j in 1:31){
    if (Ged_sort_rltv[i,j]>=0.75)
    {Ged_sort_rltv_cat[i,j]=4}
    else if (Ged_sort_rltv[i,j]<0.75 & Ged_sort_rltv[i,j]>=0.5)
    {Ged_sort_rltv_cat[i,j]=3}
    else if (Ged_sort_rltv[i,j]<0.5 & Ged_sort_rltv[i,j]>=0.25)
    {Ged_sort_rltv_cat[i,j]=2}
    else if (Ged_sort_rltv[i,j]<0.25 & Ged_sort_rltv[i,j]>0)
    {Ged_sort_rltv_cat[i,j]=1}
    else
    {Ged_sort_rltv_cat[i,j]=0}
  }
}

betad<-betadiver(t(Ged_sort_rltv_cat) ,"z")

TREAT<-sort(Treat$Treat)
adonis(betad~TREAT,method="bray",
       #strata=Treat$Block, 
       perm=1000) 
mod_cat<-betadisper(betad,TREAT)

plot(mod_cat,
     col= c("blue","green","red","orange"),
     main="Multivariate Dispersion",
     pch = c(15,22,17,24),
     segments = FALSE, # control grey line
     hull = FALSE, ellipse = FALSE,label = F,lwd=0.5)
ordiellipse(mod_cat, Treat$Treat, 
            kind="ehull", 
            conf=0.95, 
            lwd=1, 
            col= c("blue","green","red","orange")
)
boxplot(mod)


library("ade4")
names(OTU_merge_minrep)
Ged_sort_rltv_cat

OTU_pca <-Ged_sort_rltv_cat
OTU_pca <-t(OTU_pca)

boxplot(OTU_pca,las=2)
boxplot(data.frame(scale(OTU_pca)),las=2)

#OTU_transpca<-log(OTU_pca+1)

#Perform PCA
row1<-c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(6,8)) #row
pca1<-dudi.pca(OTU_pca,
               row.w=TREAT,
               scale = TRUE,
               scannf = F,
               nf=3)
#How many axes shall we keep? This is used to judge the important components.

barplot(pca1$eig)
explain<-pca1$eig/sum(pca1$eig)
barplot(explain)
#Correlation circle
s.corcircle(pca1$co)

#Interpretation
#plot()
scatter(pca1)
gcol = c("green","orange","purple","grey")
Frow1<-as.factor(c(rep("0",9),
                   rep("3",6),
                   rep("4.5",9),
                   rep("6",7)))
par(mfrow=c(2,2))
s.class(pca1$li,
        fac=Frow1,
        col = gcol, 
        xax = 1, 
        yax = 2,
        addaxes=T,
        grid=F,
        include.origin = T,
        origin=c(0,0),# the fixed point for the origin axes
        clabel = 1)
s.class(pca1$li,
        fac=Frow1,
        col = gcol, 
        xax = 2, 
        yax = 3,
        grid=F,
        addaxes=T,
        origin=c(0,0),# the fixed point for the origin axes
        clabel = 1)

s.class(pca1$li,
        fac=Frow1,
        col = gcol, 
        xax = 1, 
        yax = 3,
        addaxes=T,
        origin=c(0,0),# the fixed point for the origin axes
        grid=F,
        clabel = 1)

s.label(pca1$li,clab=0.5,add.p=T)
s.arrow(pca1$c1,add.p=T)

plot3d(pca1$li,col=gcol,type="p",radius=0.05)


### Barplot of OTUs ---------------------------------------------------------
Ged_barplot<-read.csv("Ged_barplot.csv",header=T) # import the organised dataset

Ged_barplot<-t(Ged_barplot)
Ged_barplot<-as.data.frame(Ged_barplot)
names(Ged_barplot)<-Ged_barplot[1,]
Ged_barplot<-Ged_barplot[-1,]
Ged_barplot<-decostand(Ged_barplot,"total",2)

library(reshape2)
Ged_barplot$row <- row.names(Ged_barplot)

Ged_species_SE<-t(Ged_barplot[1:31])
Ged_species_SE<-as.data.frame(Ged_species_SE)
Ged_species_SE$Treat<-sort(Treat$Treat)

Chamaedrilus_chlorophilus_se<-summarySE(data = Ged_species_SE, 
                                       measurevar = "Chamaedrilus_chlorophilus", 
                                       groupvars = c("Treat"))
Chamaedrilus_aff._sphagnetorum_A<-summarySE(data = Ged_species_SE, 
                                       measurevar = "Chamaedrilus_aff._sphagnetorum_A", 
                                       groupvars = c("Treat"))
Chamaedrilus_aff._glandulosus_B_se<-summarySE(data = Ged_species_SE, 
                                       measurevar = "Chamaedrilus_aff._glandulosus_B", 
                                       groupvars = c("Treat"))
Mesenchytraeus_cf._glandulosus_se<-summarySE(data = Ged_species_SE, 
                                       measurevar = "Mesenchytraeus_cf._glandulosus", 
                                       groupvars = c("Treat"))
Mesenchytraeus_cf._flavus_se<-summarySE(data = Ged_species_SE, 
                                       measurevar = "Mesenchytraeus_cf._flavus", 
                                       groupvars = c("Treat"))

Ged_barplot_2<-aggregate(Ged_species_SE,by=list(Ged_species_SE$Treat),mean) # calculate average
Ged_barplot_2<-Ged_barplot_2[,2:length(Ged_barplot_2)] # remove unnecessray columns
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
"#FB9A99"
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
  #theme(legend.position = "bottom", 
  #      legend.box = "vertical",
  #     legend.title = element_blank())+ # Marks the tick in y axis 
  # rescale ylim,limits=c(0.95,1),oob = rescale_none
  #guides(fill=TRUE) + # delete legend by delete the guides(fill=FALSE) line
  



barplot(Ged_barplot_2,legend=rownames(c))
