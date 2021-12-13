rm(list = ls())

library(karyoploteR)
library(VariantAnnotation)
library("viridis")        



#generate custom genome
custom.genome <- toGRanges("chrom_list.ls")


#set graphical parameters
pp <- getDefaultPlotParams(plot.type=2)
pp$ideogramheight<-60
pp$data2height<-300
pp$data1inmargin<-20
pp$data2inmargin<-140
pp$data1outmargin<-10
pp$data1height<-1200
kp <- plotKaryotype(plot.type=2, genome = custom.genome, plot.params=pp, cex=1.8)


#add bases to chromosomes
kpAddBaseNumbers(kp, tick.dist = 5000000, tick.len = 6, tick.col="black", cex=1.5,
                 minor.tick.dist = 1000000, minor.tick.len = 4, minor.tick.col = "black")


#chr_14 
#read Hifi data
snps_Hifi <- read.table("chr_14_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr14", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_14_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr14", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_14_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr14", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_14_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_14_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_14_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_14_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_15 
#read Hifi data
snps_Hifi <- read.table("chr_15_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr15", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_15_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr15", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_15_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr15", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_15_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_15_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_15_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_15_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_16 
#read Hifi data
snps_Hifi <- read.table("chr_16_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr16", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_16_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr16", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_16_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr16", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_16_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_16_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_16_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_16_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_17 
#read Hifi data
snps_Hifi <- read.table("chr_17_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr17", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_17_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr17", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_17_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr17", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_17_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_17_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_17_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_17_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")




#chr_18 
#read Hifi data
snps_Hifi <- read.table("chr_18_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr18", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_18_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr18", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_18_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr18", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_18_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_18_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_18_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_18_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_19 
#read Hifi data
snps_Hifi <- read.table("chr_19_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr19", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_19_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr19", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_19_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr19", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_19_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_19_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_19_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_19_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_20 
#read Hifi data
snps_Hifi <- read.table("chr_20_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr20", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_20_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr20", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_20_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr20", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_20_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_20_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_20_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_20_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_21 
#read Hifi data
snps_Hifi <- read.table("chr_21_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr21", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_21_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr21", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_21_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr21", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_21_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_21_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_21_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_21_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_22 
#read Hifi data
snps_Hifi <- read.table("chr_22_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr22", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_22_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr22", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_22_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr22", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_22_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_22_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_22_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_22_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_23 
#read Hifi data
snps_Hifi <- read.table("chr_23_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr23", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_23_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr23", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_23_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border="#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr23", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_23_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_23_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_23_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_23_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_24 
#read Hifi data
snps_Hifi <- read.table("chr_24_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr24", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_24_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr24", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_24_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr24", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_24_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_24_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_24_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_24_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")




#chr_25 
#read Hifi data
snps_Hifi <- read.table("chr_25_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr25", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_25_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr25", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_25_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border="#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr25", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_25_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_25_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_25_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_25_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_26
#read Hifi data
snps_Hifi <- read.table("chr_26_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr26", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_26_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr26", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_26_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr26", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_26_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_26_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_26_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_26_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_27
#read Hifi data
snps_Hifi <- read.table("chr_27_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr27", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_27_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr27", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_27_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr27", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_27_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_27_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_27_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_27_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")




#chr_28
#read Hifi data
snps_Hifi <- read.table("chr_28_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr28", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_28_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr28", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_28_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr28", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_28_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_28_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_28_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_28_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_29
#read Hifi data
snps_Hifi <- read.table("chr_29_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr29", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_29_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr29", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_29_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr29", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_29_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_29_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_29_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_29_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")



#chr_30
#read Hifi data
snps_Hifi <- read.table("chr_30_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr30", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_30_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr30", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_30_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr30", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_30_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_30_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_30_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_30_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_31
#read Hifi data
snps_Hifi <- read.table("chr_31_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr31", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_31_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr31", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_31_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr31", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_31_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_31_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_31_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_31_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_32
#read Hifi data
snps_Hifi <- read.table("chr_32_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr32", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_32_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr32", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_32_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])


#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr32", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_32_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_32_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_32_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_32_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_33
#read Hifi data
snps_Hifi <- read.table("chr_33_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr33", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_33_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr33", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_33_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr33", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_33_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_33_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_33_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_33_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_34
#read Hifi data
snps_Hifi <- read.table("chr_34_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr34", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_34_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr34", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_34_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr34", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_34_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_34_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_34_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_34_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_35
#read Hifi data
snps_Hifi <- read.table("chr_35_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr35", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_35_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr35", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_35_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr35", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_35_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_35_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_35_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_35_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_36
#read Hifi data
snps_Hifi <- read.table("chr_36_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr36", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_36_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr36", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_36_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr36", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_36_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_36_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_36_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_36_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_37
#read Hifi data
snps_Hifi <- read.table("chr_37_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr37", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_37_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr37", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_37_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr37", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_37_masked_1000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_37_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_37_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_37_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_38
#read Hifi data
snps_Hifi <- read.table("chr_38_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr38", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_38_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr38", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_38_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col="#DC3220")
kpAxis(kp, data.panel=1, chr="chr38", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_38.masked_1000bp_no_gaps.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_38_nuc_nogaps_def.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_38_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_38_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")
