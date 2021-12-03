rm(list = ls())


library(karyoploteR)
library(readxl)


#generate custom genome
custom.genome <- toGRanges("chr_list.txt")


#set graphical parameters
pp <- getDefaultPlotParams(plot.type=2)
pp$ideogramheight<-20
pp$data1inmargin<-10
pp$data2inmargin<-30


#zoom within the margins of the region of interest
zoom.region <- toGRanges(data.frame("chr6", 53260000, 54489999))

#plot karyotype 2
kp <- plotKaryotype(plot.type=2, genome = custom.genome, plot.params=pp, zoom=zoom.region, cex=1)

#add base pairs to chromosome region
kpAddBaseNumbers(kp, tick.dist=100000, tick.len=5, add.units=FALSE, digits=2, minor.ticks=FALSE, cex=.8, tick.col=NULL)

#import data with LD values (eryhtrogaster) creating a GRange object
LD_trend_Ame <- toGRanges("LD_values_Ame_m10.txt")

#import data with LD values (savignii) creating a GRange object
LD_trend_Egy <- toGRanges("LD_values_Egy_m10.txt")

#plot eryhtrogaster data
kpPoints(kp, data=LD_trend_Ame, y=LD_trend_Ame$V4, ymax=max(LD_trend_Ame$V4), cex=.5, data.panel = 1, r0=0.1, r1=0.9)
kpLines(kp, data=LD_trend_Ame, y=LD_trend_Ame$V4, ymax=max(LD_trend_Ame$V4), col="red", data.panel = 1, r0=0.1, r1=0.9)
kpAxis(kp, ymin=0, ymax=1, numticks = 5, data.panel = 1, r0=0.2, r1=0.9, cex=1.4)

#plot savignii data
kpPoints(kp, data=LD_trend_Egy, y=LD_trend_Egy$V4, ymax=max(LD_trend_Egy$V4), cex=.5, data.panel = 1, r0=0.1, r1=0.9)
kpLines(kp, data=LD_trend_Egy, y=LD_trend_Egy$V4, ymax=max(LD_trend_Egy$V4), col="green", data.panel = 1, r0=0.1, r1=0.9)

##add genes and relative markers indication on chromosome ideogram (according to NCBI annotation)
#adm gene
kpRect(kp, chr="chr6", x0=53264103, x1=53266254, y0=0, y1=1, data.panel="ideogram", col="black")

#adm marker
kpPlotMarkers(kp, chr="chr6", x=53266254, labels="ADM", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, ignore.chromosome.ends = TRUE, clipping=TRUE, label.margin=2)

#ano3 gene
kpRect(kp, chr="chr6", x0=53397976, x1=53519942, y0=0, y1=1, data.panel="ideogram", col="black")

#ano3 marker
kpPlotMarkers(kp, chr="chr6", x=53458959, labels="ANO3", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#SLC5A12 gene
kpRect(kp, chr="chr6", x0=53526610, x1=53544811, y0=0, y1=1, data.panel="ideogram", col="black")

#SLC5A12 marker
kpPlotMarkers(kp, chr="chr6", x=53535710, labels="SLC5A12", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#fibin gene
kpRect(kp, chr="chr6", x0=53625483, x1=53630333, y0=0, y1=1, data.panel="ideogram", col="black")

#fibin marker
kpPlotMarkers(kp, chr="chr6", x=53625000, labels ="FIBIN", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#bbox1 gene 
kpRect(kp, chr="chr6", x0=53635630, x1=53670240, y0=0, y1=1, data.panel="ideogram", col="black")

#bbox1 marker
kpPlotMarkers(kp, chr="chr6", x=53661900, labels="BBOX1", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#ccdc34 gene
kpRect(kp, chr="chr6", x0=53718663, x1=53740082, y0=0, y1=1, data.panel="ideogram", col="black")

#ccdc34 marker
kpPlotMarkers(kp, chr="chr6", x=53729372, labels="CCDC34", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#lgr4 gene
kpRect(kp, chr="chr6", x0=53746838, x1=53819814, y0=0, y1=1, data.panel="ideogram", col="black")

#lgr4 marker
kpPlotMarkers(kp, chr="chr6", x=53783326, labels="LGR4", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#lin7c gene
kpRect(kp, chr="chr6", x0=53838968, x1=53852117, y0=0, y1=1, data.panel="ideogram", col="black")

#lin7c marker
kpPlotMarkers(kp, chr="chr6", x=53845542, labels="LIN7C", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#bdnf gene
kpRect(kp, chr="chr6", x0=53886627, x1=53927580, y0=0, y1=1, data.panel="ideogram", col="black")

#bdnf marker
kpPlotMarkers(kp, chr="chr6", x=53907103, labels="BDNF", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#kif18a gene
kpRect(kp, chr="chr6", x0=54034054, x1=54065099, y0=0, y1=1, data.panel="ideogram", col="black")

#kif18a marker
kpPlotMarkers(kp, chr="chr6", x=54049576, labels="KIF18A", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#METTL15 gene
kpRect(kp, chr="chr6", x0=54065004, x1=54157704, y0=0, y1=1, data.panel="ideogram", col="black")

#METTL15 marker
kpPlotMarkers(kp, chr="chr6", x=54111354, labels="METTL15", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#plot low confidence region on the ideogram
low_confidence<-toGRanges("chr6_genomecov_high_low_cov_intervals.bed")
kpPlotRegions(kp, data.panel="ideogram", data=low_confidence, col="red")

#add low confidence region marker
kpPlotMarkers(kp, chr="chr6", x=53329000, labels="low confidence region", text.orientation = "horizontal", cex=.8, data.panel=1, r0=0, r1=0.1, label.margin=2)

#add r2 label on plot
kpAddLabels(kp, labels="r^2", data.panel = 1, r0=0.97, r1=0.99, cex=1.4)

#plot eryhtrogaster SNP counts as heatmaps
ame_snp_counts <- toGRanges("Ame_53260000-54489999_snp_counts")
kpHeatmap(kp, data.panel=2, r0=0, r1=0.15, data=ame_snp_counts, y=ame_snp_counts$V4, col=c("white", "red"))

#plot savignii SNP counts as heatmaps
egy_snp_counts <- toGRanges("Egy_53260000-54489999_snp_counts")
kpHeatmap(kp, data.panel=2, r0=0.15, r1=0.3, data=egy_snp_counts, y=egy_snp_counts$V4, col=c("white", "green"))

#add SNP counts label
kpAddLabels (kp, data.panel=2, labels = "SNP counts", r0=0.13, r1=0.17, cex=1.4)

#legend
legend(0.75, 1, fill = c("red", "green"), legend = c("H. r. erythrogaster", "H. r. savignii"), cex=1.1, bty="n") 


