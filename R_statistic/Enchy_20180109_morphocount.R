setwd("/Users/JiayiQin/Dropbox/2nd_Miseq/Morpho_count")

Enchy<-read.csv("Enchytraeid_sum.csv",header = T, sep =",",check.names=FALSE)

Enchy$Treatment<-factor(Enchy$Treatment) # adjust the variable tye
Enchy$Number<-as.numeric(Enchy$Number) 
Enchy$Total<-as.numeric(Enchy$Total)

# convert the counts into enchytraeid density
r <- 0.058 / 2 # soil core diameter is 5.8 cm 


Enchy$Number <- Enchy$Number/(3.1415926 *r*r)
Enchy$Total <- Enchy$Total/(3.1415926 *r*r)
Treatment<-factor(Enchy$Treatment)

Enchy[,2]<-Treatment
Time<-Enchy$Time

# reference
# http://stats.stackexchange.com/questions/159947/using-lmer-for-repeated-measures-and-rcbd
# http://stats.stackexchange.com/questions/58745/using-lmer-for-repeated-measures-linear-mixed-effect-model
# Plots are nested within Blocks (random). 
# Time is not interact with Treatment
# Repeated measurement indicated by Time has 
# https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q4/004530.html
# warning message produced, therefore run the code to check what happened
# reference: https://github.com/lme4/lme4/issues/120
# https://stats.stackexchange.com/questions/64535/structuring-a-linear-mixed-model-in-r-with-nesting

########## stat Top& Bottom together #######
### entire data set ####
library(lme4)
library(lsmeans)
# lsmeans seems to be a correct way to perform pairwase interractions.
# But don't forget to correct your P-values accoring to the interractions you get.
# Also, if the FACTOR1*FACTOR2 is not significant, you shoud put a bar | instead of a star * like this: 
# lsmeans(model, pairwise ~ FACTOR1 | FACTOR2, adjust = "tukey").
# ResearchGate. Available from: https://www.researchgate.net/post/Post_hoc_test_in_linear_mixed_models_how_to_do [accessed Jan 3, 2017].

plot(Enchy$Number)
Enchy$Treatment<-factor(Enchy$Treatment)

#### Test the effect of Treatment*Time*Position #######
Number_1 <- lmer(log(Enchy$Number+2000)~Treatment*Time*Position + (Time|Block/Subsamples),data=Enchy,REML=FALSE)
shapiro.test(residuals(Number_1))
Number_2<- lmer(log(Enchy$Number+2000)~Time*Treatment+Position+(Time|Block/Subsamples),data=Enchy,REML=FALSE)
shapiro.test(residuals(Number_2))
anova(Number_1,Number_2) # test interation with position, significant

Number_3<- lmer(log(Enchy$Number+2000)~Time*Position+(Time|Block/Subsamples),
                data=Enchy,REML=FALSE)
shapiro.test(residuals(Number_3))
anova(Number_3,Number_1) # interaction with treatment, significant

Number_4<- lmer(log(Enchy$Number+2000)~Time*Treatment+(Time|Block/Subsamples),data=Enchy,REML=FALSE)
shapiro.test(residuals(Number_4))
anova(Number_1,Number_4)# position, significant

Number_5<- lmer(log(Enchy$Number+2000) ~ Treatment*Position+Time+(1|Block/Subsamples),data=Enchy,REML=FALSE)
shapiro.test(residuals(Number_5))
anova(Number_5,Number_1) # the time is not interacted with others.

Number_6<- lmer(log(Enchy$Number+2000)~Treatment*Position+(1|Block/Subsamples),data=Enchy,REML=FALSE)
shapiro.test(residuals(Number_6))
anova(Number_1,Number_6) # time factor significant


Enchy_total<-Enchy[1:144,]


# Subset_Jun2014 ----------------------------------------------------------

Enchy_Jun2014<-subset(Enchy,Time=="Jun-14")

shapiro.test(log(Enchy_Jun2014$Number+1000))
Number_Jun2014_1 <- lmer(log(Enchy_Jun2014$Number+1000) ~ Treatment*Position+(1|Block/Subsamples),data=Enchy_Jun2014,REML=FALSE)
shapiro.test(residuals(Number_Jun2014_1))

Number_Jun2014_2 <- lmer(log(Enchy_Jun2014$Number+1000) ~ Position + (1|Block/Subsamples),data=Enchy_Jun2014,REML=FALSE)
shapiro.test(residuals(Number_Jun2014_2))

anova(Number_Jun2014_1,Number_Jun2014_2) # not significant

# Subset_Oct2014 ----------------------------------------------------------

Enchy_Oct2014<-subset(Enchy,Time=="Oct-14")

Number_Oct2014_1 <- lmer(log(Enchy_Oct2014$Number+1000)~Treatment*Position+(1|Block/Subsamples),data=Enchy_Oct2014,REML=FALSE)
shapiro.test(residuals(Number_Oct2014_1))
Number_Oct2014_2 <- lmer(log(Enchy_Oct2014$Number+1000)~Treatment+Position+(1|Block/Subsamples),data=Enchy_Oct2014,REML=FALSE)
shapiro.test(residuals(Number_Oct2014_2))
Number_Oct2014_3 <- lmer(log(Enchy_Oct2014$Number+1000)~Position+(1|Block/Subsamples),data=Enchy_Oct2014,REML=FALSE)
shapiro.test(residuals(Number_Oct2014_3))

anova(Number_Oct2014_1,Number_Oct2014_2) #tedency of significant
anova(Number_Oct2014_2,Number_Oct2014_3) #significant
anova(Number_Oct2014_1,Number_Oct2014_3) #significant

ls_Number_Oct2014<-lsmeans(Number_Oct2014_1,pairwise~Treatment*Position,at = list(Position = "Top"), adjust="tukey")
ls_Number_Oct2014

ls_Number_Oct2014<-lsmeans(Number_Oct2014_1,pairwise~Treatment*Position,at = list(Position = "Bottom"), adjust="tukey")
ls_Number_Oct2014

##### EC50 Oct2014 #########
library (drc)
UPPER<- max(Enchy_Oct2014$Number[1:36])
Enchy_Oct2014_drc<-Enchy_Oct2014[1:36,]
Enchy_Oct2014_drc$Treatment<-as.numeric(paste(Enchy_Oct2014_drc$Treatment))
Enchy_Oct2014_drc$Total<-as.numeric(paste(Enchy_Oct2014_drc$Total))
Enchy_Oct2014_drc$Number<-as.numeric(paste(Enchy_Oct2014_drc$Number))

drc_Oct2014<- drm(Total~Treatment,
                  data=Enchy_Oct2014_drc,
                  fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "ED50")),
                  control=drmc(),
                  lowerl = c(-Inf,0, 0, -Inf), #Limitation for parameters
                  upperl = c(Inf, Inf, UPPER, Inf)
)

summary(drc_Oct2014)
ED(drc_Oct2014,c(10,50,90),interval="delta")
plot(drc_Oct2014,xlim=c(0,10))
plot(Total~Treatment, data=Enchy_Oct2014_drc)
# Subset_Jun2015 ----------------------------------------------------------

Enchy_Jun2015<-subset(Enchy,Time=="Jun-15")

Number_Jun2015_1 <- lmer(log(Enchy_Jun2015$Number+3000)~Treatment*Position+(1|Block/Subsamples),data=Enchy_Jun2015,REML=FALSE)
shapiro.test(residuals(Number_Jun2015_1))
Number_Jun2015_2 <- lmer(log(Enchy_Jun2015$Number+3000)~Position+(1|Block/Subsamples),data=Enchy_Jun2015,REML=FALSE)
shapiro.test(residuals(Number_Jun2015_2))
anova(Number_Jun2015_1,Number_Jun2015_2) # not significant


# Subset_Jun2016 ----------------------------------------------------------

Enchy_Jun2016<-subset(Enchy,Time=="Jun-16")

Number_Jun2016_1 <- lmer(log(Enchy_Jun2016$Number+1000)~Treatment*Position+(1|Block/Subsamples),data=Enchy_Jun2016,REML=FALSE)
shapiro.test(residuals(Number_Jun2016_1))
Number_Jun2016_2 <- lmer(log(Enchy_Jun2016$Number+1000)~Position+(1|Block/Subsamples),data=Enchy_Jun2016,REML=FALSE)
shapiro.test(residuals(Number_Jun2016_2))
anova(Number_Jun2016_1,Number_Jun2016_2) #significant

ls_Number_Jun2016<-lsmeans(Number_Jun2016_1,pairwise~Treatment*Position,at = list(Position = "Top"), adjust="tukey")
ls_Number_Jun2016
ls_Number_Jun2016<-lsmeans(Number_Jun2016_1,pairwise~Treatment*Position,at = list(Position = "Bottom"), adjust="tukey")
ls_Number_Jun2016

##Jun2016 EC 50######################
library (drc)
Enchy_Jun2016<-subset(Enchy,Time=="Jun-16")
UPPER<- max(Enchy_Jun2016$Total[1:36])
Enchy_Jun2016_drc<-Enchy_Jun2016[1:36,]
Enchy_Jun2016_drc$Treatment<-as.numeric(paste(Enchy_Jun2016_drc$Treatment))
Enchy_Jun2016_drc$Total<-as.numeric(paste(Enchy_Jun2016_drc$Total))

drc_Jun2016<- drm(Total~Treatment,
                   data=Enchy_Jun2016_drc,
                   fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "ED50")),
                   control=drmc(),
                   lowerl = c(-Inf,0, -Inf, -Inf), #Limitation for parameters
                   upperl = c(Inf, Inf, UPPER, Inf)
)

summary(drc_Jun2016)
ED(drc_Jun2016,c(10,50,90),interval="delta")

############# Seperate Top and Bottom ##########

# Prepare dataset ---------------------------------------------------------

Oct14_Enchy<-subset(Enchy, Time=="Oct-14")
Jun14_Enchy<-subset(Enchy, Time=="Jun-14")
Jun15_Enchy<-subset(Enchy, Time=="Jun-15")
Jun16_Enchy<-subset(Enchy, Time=="Jun-16")

# Pair wise test ----------------------------------------------------------
#### Jun14_Enchy ####
Jun14_Enchy_Top<-subset(Enchy_Jun2014,Position=="Top")
Jun14_Enchy_Bottom<-subset(Enchy_Jun2014,Position=="Bottom")
shapiro.test(log(Jun14_Enchy_Top$Number+1000)) #Shapiro-Wilk normality yes
shapiro.test(log(Jun14_Enchy_Bottom$Number+1000)) #Shapiro-Wilk normality yes
shapiro.test(log(Jun14_Enchy$Total+1000)) #Shapiro-Wilk normality yes

Jun14_Enchy_Top_aov_1<-lmer(log(Jun14_Enchy_Top$Number+1000)~Treatment+(1|Block/Subsamples),data=Jun14_Enchy_Top)
shapiro.test(residuals(Jun14_Enchy_Top_aov_1))
Jun14_Enchy_Top_aov_2<-lmer(log(Jun14_Enchy_Top$Number+1000)~(1|Block/Subsamples),data=Jun14_Enchy_Top)
shapiro.test(residuals(Jun14_Enchy_Top_aov_2))
anova(Jun14_Enchy_Top_aov_1,Jun14_Enchy_Top_aov_2) # not significant - moderately affected?

ls_Jun14_Enchy_Top<-lsmeans(Jun14_Enchy_Top_aov_1,pairwise~Treatment,data=Jun14_Enchy_Top,adjust="tukey")
ls_Jun14_Enchy_Top


Jun14_Enchy_Bottom_aov_1<-lmer(log(Jun14_Enchy_Bottom$Number+1000)~Treatment+(1|Block/Subsamples),data=Jun14_Enchy_Bottom)
shapiro.test(residuals(Jun14_Enchy_Bottom_aov_1))
Jun14_Enchy_Bottom_aov_2<-lmer(log(Jun14_Enchy_Bottom$Number+1000)~(1|Block/Subsamples),data=Jun14_Enchy_Bottom)
shapiro.test(residuals(Jun14_Enchy_Bottom_aov_2))

anova(Jun14_Enchy_Bottom_aov_1,Jun14_Enchy_Bottom_aov_2)  # not significant

ls_Jun14_Enchy_Bottom<-lsmeans(Jun14_Enchy_Bottom_aov_1,pairwise~Treatment,data=Jun14_Enchy_Bottom,adjust="tukey")
ls_Jun14_Enchy_Bottom


#### Oct14_Enchy ####
Oct14_Enchy_Top<-subset(Oct14_Enchy,Position=="Top")
Oct14_Enchy_Bottom<-subset(Oct14_Enchy,Position=="Bottom")

Oct14_Enchy_Top_aov_1<-lmer(log(Oct14_Enchy_Top$Number+2000)~Treatment+(1|Block/Subsamples),data=Oct14_Enchy_Top)
shapiro.test(residuals(Oct14_Enchy_Top_aov_1))
Oct14_Enchy_Top_aov_2<-lmer(log(Oct14_Enchy_Top$Number+2000)~(1|Block/Subsamples),data=Oct14_Enchy_Top)
shapiro.test(residuals(Oct14_Enchy_Top_aov_2))
anova(Oct14_Enchy_Top_aov_1,Oct14_Enchy_Top_aov_2) # significant

ls_Oct14_Enchy_Top<-lsmeans(Oct14_Enchy_Top_aov_1,pairwise~Treatment,data=Oct14_Enchy_Top,adjust="tukey")
ls_Oct14_Enchy_Top

Oct14_Enchy_Bottom_aov_1<-lmer(log(Oct14_Enchy_Bottom$Number+500)~Treatment+(1|Block/Subsamples),data=Oct14_Enchy_Bottom)
shapiro.test(residuals(Oct14_Enchy_Bottom_aov_1))
Oct14_Enchy_Bottom_aov_2<-lmer(log(Oct14_Enchy_Bottom$Number+500)~(1|Block/Subsamples),data=Oct14_Enchy_Bottom)
shapiro.test(residuals(Oct14_Enchy_Bottom_aov_2))
anova(Oct14_Enchy_Bottom_aov_1,Oct14_Enchy_Bottom_aov_2) # not significant

#
Oct14_Enchy_Total_aov_1<-lmer(log(Oct14_Enchy_Top$Total+1000)~Treatment+(1|Block/Subsamples),data=Oct14_Enchy_Top)
shapiro.test(residuals(Oct14_Enchy_Total_aov_1))
Oct14_Enchy_Total_aov_2<-lmer(log(Oct14_Enchy_Top$Total+1000)~(1|Block/Subsamples),data=Oct14_Enchy_Top)
shapiro.test(residuals(Oct14_Enchy_Total_aov_2))
anova(Oct14_Enchy_Total_aov_1,Oct14_Enchy_Total_aov_2) # not significant

lsmeans(Oct14_Enchy_Total_aov_1,pairwise~Treatment,data=Oct14_Enchy_Top,adjust="tukey") # no significant pairwise comparison
Dunnett.Oct14_Enchy_Total<-glht(Oct14_Enchy_Total_aov_1,linfct=mcp(Treatment="Dunnett"),alternative = "less")
summary(Dunnett.Oct14_Enchy_Total)

Tukey.Oct14_Enchy_Total<-glht(Oct14_Enchy_Total_aov_1,linfct=mcp(Treatment="Tukey"),alternative = "less")
summary(Tukey.Oct14_Enchy_Total)


#### Jun15_Enchy ####
Jun15_Enchy_Top<-subset(Jun15_Enchy,Position=="Top")
Jun15_Enchy_Bottom<-subset(Jun15_Enchy,Position=="Bottom")

Jun15_Enchy_Top_aov_1<-lmer(log(Jun15_Enchy_Top$Number+2000)~Treatment+(1|Block/Subsamples),data=Jun15_Enchy_Top)
shapiro.test(residuals(Jun15_Enchy_Top_aov_1))
Jun15_Enchy_Top_aov_2<-lmer(log(Jun15_Enchy_Top$Number+2000)~(1|Block/Subsamples),data=Jun15_Enchy_Top)
shapiro.test(residuals(Jun15_Enchy_Top_aov_2))
anova(Jun15_Enchy_Top_aov_1,Jun15_Enchy_Top_aov_2) # not significant

Jun15_Enchy_Bottom_aov_1<-lmer(log(Jun15_Enchy_Bottom$Number+1000)~Treatment+(1|Block/Subsamples),data=Jun15_Enchy_Bottom)
shapiro.test(residuals(Jun15_Enchy_Bottom_aov_1))

Jun15_Enchy_Bottom_aov_2<-lmer(log(Jun15_Enchy_Bottom$Number+1000)~(1|Block/Subsamples),data=Jun15_Enchy_Bottom)
shapiro.test(residuals(Jun15_Enchy_Top_aov_2))

anova(Jun15_Enchy_Bottom_aov_1,Jun15_Enchy_Bottom_aov_2)# not significant


#### Jun16_Enchy ####

Jun16_Enchy_Top<-subset(Enchy_Jun2016,Position=="Top")
Jun16_Enchy_Bottom<-subset(Enchy_Jun2016,Position=="Bottom")

Jun16_Enchy_Top_aov_1<-lmer(log(Jun16_Enchy_Top$Number+1000)~Treatment+(1|Block/Subsamples),data=Jun16_Enchy_Top)
shapiro.test(residuals(Jun16_Enchy_Top_aov_1))
Jun16_Enchy_Top_aov_2<-lmer(log(Jun16_Enchy_Top$Number+1000)~(1|Block/Subsamples),data=Jun16_Enchy_Top)
shapiro.test(residuals(Jun16_Enchy_Top_aov_2))
anova(Jun16_Enchy_Top_aov_1,Jun16_Enchy_Top_aov_2)  # not significant

Jun16_Enchy_Bottom_aov_1<-lmer(log(Jun16_Enchy_Bottom$Number+1000)~Treatment+(1|Block/Subsamples),data=Jun16_Enchy_Bottom)
shapiro.test(residuals(Jun16_Enchy_Bottom_aov_1))

Jun16_Enchy_Bottom_aov_2<-lmer(log(Jun16_Enchy_Bottom$Number+1000)~(1|Block/Subsamples),data=Jun16_Enchy_Bottom)
shapiro.test(residuals(Jun16_Enchy_Bottom_aov_2))
anova(Jun16_Enchy_Bottom_aov_1,Jun16_Enchy_Bottom_aov_2)
ls_Jun16_Enchy_Bottom<-lsmeans(Jun16_Enchy_Bottom_aov,pairwise~Treatment,data=Jun16_Enchy_Bottom,adjust="tukey")
ls_Jun16_Enchy_Bottom
summary(Jun16_Enchy_Bottom_aov) # not significant

Jun16_Enchy_Total_aov_1<-lmer(log(Jun16_Enchy_Top$Total+10)~Treatment+(1|Block/Subsamples),data=Jun16_Enchy_Top)
shapiro.test(residuals(Jun16_Enchy_Total_aov_1))
Jun16_Enchy_Total_aov_2<-lmer(log(Jun16_Enchy_Top$Total+10)~(1|Block/Subsamples),data=Jun16_Enchy_Top)
shapiro.test(residuals(Jun16_Enchy_Total_aov_2))
anova(Jun16_Enchy_Total_aov_1,Jun16_Enchy_Total_aov_2)
# follow with multicomparison
lsmeans(Jun16_Enchy_Total_aov_1,pairwise~Treatment,data=Jun16_Enchy_Top,adjust="tukey")
library(multcomp)
Dunnett.Jun16_Enchy_Total<-glht(Jun16_Enchy_Total_aov_1,linfct=mcp(Treatment="Dunnett"),alternative = "less")
summary(Dunnett.Jun16_Enchy_Total)

Tukey.Jun16_Enchy_Total<-glht(Jun16_Enchy_Total_aov_1,linfct=mcp(Treatment="Tukey"),alternative = "less")
summary(Tukey.Jun16_Enchy_Total)



# Graph_1 -----------------------------------------------------------------


Oct14_Enchy<-subset(Enchy, Time=="Oct-14")
Jun14_Enchy<-subset(Enchy, Time=="Jun-14")
Jun15_Enchy<-subset(Enchy, Time=="Jun-15")
Jun16_Enchy<-subset(Enchy, Time=="Jun-16")

library(Rmisc)

Oct14_Enchy_Top_dt <- summarySE(data = Oct14_Enchy_Top, 
                               measurevar = "Number", 
                               groupvars = c("Treatment", "Time"))
Jun14_Enchy_Top_dt <- summarySE(data = Jun14_Enchy_Top, 
                                measurevar = "Number",
                                groupvars = c("Treatment", "Time"))
Jun15_Enchy_Top_dt <- summarySE(data = Jun15_Enchy_Top, 
                                measurevar = "Number",
                                groupvars = c("Treatment", "Time"))
Jun16_Enchy_Top_dt <- summarySE(data = Jun16_Enchy_Top, 
                                measurevar = "Number",
                                groupvars = c("Treatment", "Time"))


Oct14_Enchy_Subtop_dt <- summarySE(data = Oct14_Enchy_Bottom, 
                                measurevar = "Number",
                                groupvars = c("Treatment", "Time"))
Jun14_Enchy_Subtop_dt <- summarySE(data = Jun14_Enchy_Bottom, 
                                measurevar = "Number",
                                groupvars = c("Treatment", "Time"))
Jun15_Enchy_Subtop_dt <- summarySE(data = Jun15_Enchy_Bottom, 
                                measurevar = "Number", 
                                groupvars = c("Treatment", "Time"))
Jun16_Enchy_Subtop_dt <- summarySE(data = Jun16_Enchy_Bottom, 
                                   measurevar = "Number", 
                                groupvars = c("Treatment", "Time"))

Oct14_Enchy_Total_dt <- summarySE(data = Oct14_Enchy_Top, 
                                   measurevar = "Total", 
                                   groupvars = c("Treatment", "Time"))
Jun14_Enchy_Total_dt <- summarySE(data = Jun14_Enchy_Top, 
                                   measurevar = "Total", 
                                   groupvars = c("Treatment", "Time"))
Jun15_Enchy_Total_dt <- summarySE(data = Jun15_Enchy_Top, 
                                   measurevar = "Total", 
                                   groupvars = c("Treatment", "Time"))
Jun16_Enchy_Total_dt <- summarySE(data = Jun16_Enchy_Top, 
                                   measurevar = "Total", 
                                   groupvars = c("Treatment", "Time"))



par (mar = c(2,4,4,2))
#template
library(gplots)
{
  plotCI(x = Oct14_Enchy_Top_dt$Top, 
       uiw = Oct14_Enchy_Top_dt$se, 
       xaxt ="n", 
       las = 1,
       #xlim = c(0.5,2.5), 
       ylim = c(0,100), 
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "blue")
}
##########################################################
# ouput graphic
postscript(file = "/Users/JiayiQin/Dropbox/2nd_Miseq/ash_effect_on_enchy_manus/Fig3_stackbarplot.eps",
           width = 5.5, height = 4,pointsize = 7,
           bg = "white")

par(mfrow=c(3,4))

par(cex = 0.8)
par(mar = c(0, 0, 0, 0), oma = c(4, 5, 4, 0.5))
par(tcl = -0.25)

plotCI(x = Jun14_Enchy_Top_dt$Number, 
       uiw = Jun14_Enchy_Top_dt$se, 
       xaxt ="n", 
       las = 1,
       #xlim = c(0.5,2.5), 
       ylim = c(0,50000), 
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,
       col = "black")
mtext("Upper", font=2,side=2,line=3.5,cex =1)
mtext("June,2014",font=2,side=3,line=1.5,cex =1)

plotCI(x = Oct14_Enchy_Top_dt$Number, 
       uiw = Oct14_Enchy_Top_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim = c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 
#text(Oct14_Enchy_Top_dt$Treatment,Oct14_Enchy_Top_dt$Top+30,labels=c("a","a","b","ab"))
mtext("October,2014",font=2,side=3,line=1.5,cex =1)

plotCI(x = Jun15_Enchy_Top_dt$Number, 
       uiw = Jun15_Enchy_Top_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 
mtext("June,2015",font=2,side=3,line=1.5,cex =1)

plotCI(x = Jun16_Enchy_Top_dt$Number, 
       uiw = Jun16_Enchy_Top_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 
mtext("June,2016",font=2,side=3,line=1.5,cex =1)

plotCI(x = Jun14_Enchy_Subtop_dt$Number, 
       uiw = Jun14_Enchy_Subtop_dt$se, 
       las = 1,
       xaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 
mtext("Lower", font=2,side=2,line=3.5,cex =1)

plotCI(x = Oct14_Enchy_Subtop_dt$Number, 
       uiw = Oct14_Enchy_Subtop_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 
plotCI(x = Jun15_Enchy_Subtop_dt$Number, 
       uiw = Jun15_Enchy_Subtop_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 

plotCI(x = Jun16_Enchy_Subtop_dt$Number, 
       uiw = Jun16_Enchy_Subtop_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 

plotCI(x = Jun14_Enchy_Total_dt$Total, 
       uiw = Jun14_Enchy_Total_dt$se, 
       las = 1,
       xaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 
axis(1, at = c(1,2,3,4), label = c(0,3,4.5,6),cex=1)
mtext("Total", font=2,side=2,line=3.5,cex =1)

plotCI(x = Oct14_Enchy_Total_dt$Total, 
       uiw = Oct14_Enchy_Total_dt$se, 
       yaxt="n",
       xaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black")
axis(1, at = c(1,2,3,4), label = c(0,3,4.5,6),cex=1)
plotCI(x = Jun15_Enchy_Total_dt$Total, 
       uiw = Jun15_Enchy_Total_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black")
axis(1, at = c(1,2,3,4), label = c(0,3,4.5,6),cex=1)
plotCI(x = Jun16_Enchy_Total_dt$Total, 
       uiw = Jun16_Enchy_Total_dt$se,
       xaxt="n",
       yaxt="n",
       ylim=c(0,50000),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       #labels=c("a","ab","b","b"),
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black")
#text(Jun16_Enchy_Total_dt$Treatment,Jun16_Enchy_Total_dt$Total+30,labels=c("a","ab","b","b"))
axis(1, at = c(1,2,3,4), label = c(0,3,4.5,6),cex=1)



































