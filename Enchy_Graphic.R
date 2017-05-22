###############################################################################################
#                                                                                             #  
#                                                                                             #
#                                                                                             #
#                             Graphic for Definitive test                                     #
#                                                                                             #
#                               Export in Feb 22, 2016  
#
#   tiff command at the beginning of the script                                                                                       #
#   use dev.off() to close the edit and save image into an independent file
#   DPI are important, setting by tiff(res=)
###############################################################################################


# line for enchytraeid count ----------------------------------------------


library(gplots)
postscript(file = "/Users/JiayiQin/Dropbox/PhD report/Graph/Fig5_enchy_count.eps",
           width = 5.85, height = 4.15,pointsize = 7,
           bg = "white")
par(mfrow=c(3,4))

par(cex = 0.8)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 1.5))
par(tcl = -0.25)

plotCI(x = Jun14_Enchy_Top_dt$Number, 
       uiw = Jun14_Enchy_Top_dt$se, 
       xaxt ="n", 
       las = 1,
       #xlim = c(0.5,2.5), 
       ylim = c(0,150), 
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,
       col = "black")
mtext("Top", font=2,side=2,line=2.5,cex =1)
mtext("June,2014",font=2,side=3,line=1.5,cex =1)

plotCI(x = Oct14_Enchy_Top_dt$Number, 
       uiw = Oct14_Enchy_Top_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim=c(0,150),
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
       ylim=c(0,150),
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
       ylim=c(0,150),
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
       ylim=c(0,150),
       gap = 0, 
       ylab="", 
       xlab="",
       main = "",
       lty = 1,
       pch = 21,
       type = "l",
       lwd = 1,cex=1,cex.axis=1,
       col = "black") 
mtext("Bottom", font=2,side=2,line=2.5,cex =1)

plotCI(x = Oct14_Enchy_Subtop_dt$Number, 
       uiw = Oct14_Enchy_Subtop_dt$se, 
       xaxt="n",
       yaxt="n",
       ylim=c(0,150),
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
       ylim=c(0,150),
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
       ylim=c(0,150),
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
       ylim=c(0,150),
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
mtext("Total", font=2,side=2,line=2.5,cex =1)

plotCI(x = Oct14_Enchy_Total_dt$Total, 
       uiw = Oct14_Enchy_Total_dt$se, 
       yaxt="n",
       xaxt="n",
       ylim=c(0,150),
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
       ylim=c(0,150),
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
       ylim=c(0,150),
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

dev.off()
 

# stachbarplot ------------------------------------------------------------
library(ggplot2)
require(scales)
library(plyr)
postscript(file = "/Users/JiayiQin/Dropbox/PhD report/Graph/Fig6_enchy_sbp.eps",
           width =7, height = 4.15,pointsize = 7,
           bg = "white")

fill<-c("#33A02C","#B2DF8A","#1F78B4","#A6CEE3","#E31A1C","#FDBF6F")

Ged_barplot_2$Treat<-as.character(paste(Ged_barplot_2$Treat))
Ged_barplot_2$variable<-gsub("Chamaedrilus_aff._sphagnetorum_A","Chamaedrilus sphagnetorum",Ged_barplot_2$variable)
Ged_barplot_2$variable<-gsub("_"," ",Ged_barplot_2$variable)
colnames(Ged_barplot_2)<-c("Treat","Species","value")
stackbarplot<- ggplot(Ged_barplot_2, aes(x=Treat, y=value,fill=Species) )+ 
  geom_bar(stat="identity") +
  xlab("\nTreatment") +
  ylab("Percent\n") +
  scale_y_continuous(labels = percent_format(),
                     #limits=c(0.7,1),
                     oob = rescale_none) +
  theme(legend.key.size=unit(0.3,"inch"))+
  scale_fill_manual(values=fill) #+
  #theme_bw()
stackbarplot
dev.off()


# PCoA --------------------------------------------------------------------
tiff(res=600)
cairo_ps(file = "/Users/JiayiQin/Dropbox/PhD report/Graph/Fig7_enchy_pcoa.eps",
           width = 3.54331, height = 2.3622,pointsize = 7,
           fallback_resolution = 600)

par(cex = 0.8)
par(mfrow = c(1,1))
par(mgp=c(3, 0.7, 0), plt=c(0.2, 1, 0.2, 1), omd=c(0,0.9,0,0.9), family="Helvetica")

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
            lwd=0.8, 
            col= c("blue","green","red","orange")
)
legend("bottomleft", c("0 t/ha", "3 t/ha", "4.5 t/ha", "6 t/ha"), 
       col = c("blue","green","red","orange"),
       title="incidence dataset",
       text.col = "black", 
       # lty = c(1, 1, 1,1),
       pch = c(15,22,17,24),
       bty="n")
dev.off()
