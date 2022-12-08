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
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 6, tick.col="black", cex=1.5,
                 minor.tick.dist = 1000000, minor.tick.len = 4, minor.tick.col = "black")


#chr_1 
#read Hifi data
snps_Hifi <- read.table("chr_1_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr1", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_1_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr1", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_1_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr1", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_1_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_1_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_1_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))


#plot gaps on ideogram as regions
gaps<-toGRanges("chr_1_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_2 
#read Hifi data
snps_Hifi <- read.table("chr_2_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr2", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_2_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr2", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_2_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr2", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_2_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_2_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_2_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_2_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_3 
#read Hifi data
snps_Hifi <- read.table("chr_3_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr3", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_3_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr3", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_3_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr3", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_3_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_3_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_3_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_3_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_Z
#read Hifi data
snps_Hifi <- read.table("chr_Z_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chrZ", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_Z_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chrZ", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_Z_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chrZ", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_Z_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_Z_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_Z_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_Z_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_4
#read Hifi data
snps_Hifi <- read.table("chr_4_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr4", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_4_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr4", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_4_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity 
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr4", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_4_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_4_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_4_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_4_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_5
#read Hifi data
snps_Hifi <- read.table("chr_5_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr5", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_5_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr5", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_5_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr5", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_5_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_5_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_5_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_5_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_6
#read Hifi data
snps_Hifi <- read.table("chr_6_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr6", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_6_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr6", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_6_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr6", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_6_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_6_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_6_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_6_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_7
#read Hifi data
snps_Hifi <- read.table("chr_7_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr7", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_7_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr7", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_7_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr7", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)

#plot repeated regions as regions
repeats <- toGRanges("chr_7_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_7_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_7_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_7_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_8
#read Hifi data
snps_Hifi <- read.table("chr_8_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr8", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_8_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr8", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_8_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr8", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_8_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_8_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_8_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_8_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_W
#read Hifi data
snps_Hifi <- read.table("chr_W_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chrW", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_W_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chrW", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_W_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chrW", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_W_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_W_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_W_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_W_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_9
#read Hifi data
snps_Hifi <- read.table("chr_9_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr9", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_9_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr9", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_9_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr9", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_9_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_9_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_9_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_9_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_10
#read Hifi data
snps_Hifi <- read.table("chr_10_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr10", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_10_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr10", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_10_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr10", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_10_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_10_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_10_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_10_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_11
#read Hifi data
snps_Hifi <- read.table("chr_11_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr11", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_11_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr11", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)


#read ddRAD data
snps_ddRAD <- read.table("chr_11_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr11", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_11_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_11_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_11_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_11_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")


#chr_12
#read Hifi data
snps_Hifi <- read.table("chr_12_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr12", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_12_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr12", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)

#read ddRAD data
snps_ddRAD <- read.table("chr_12_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr12", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_12_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_12_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_12_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_12_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")

#chr_13
#read Hifi data
snps_Hifi <- read.table("chr_13_snp_pos_def_Hifi20X", sep="", header=FALSE)

#create GRange object
snps_Hifi <- toGRanges(snps_Hifi[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_Hifi, r0=0, r1=0.25, window.size=40000, border="#85C0F9", col="#85C0F9")
kpAxis(kp, data.panel=1, chr="chr13", ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.25, cex=1.8)

#read WGS data
snps_WGS <- read.table("chr_13_snp_pos_def_WGS", sep="", header=FALSE)

#create GRange object
snps_WGS <- toGRanges(snps_WGS[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_WGS, r0=0.35, r1=0.6, window.size=40000, border="#005AB5", col="#005AB5")
kpAxis(kp, data.panel=1, chr="chr13", ymax=kp$latest.plot$computed.values$max.density, r0=0.35, r1=0.6, cex=1.8)

#read ddRAD data
snps_ddRAD <- read.table("chr_13_snp_pos_def_ddRAD", sep="", header=FALSE)

#create GRange object
snps_ddRAD <- toGRanges(snps_ddRAD[,c(1,2,2)])

#kpPlotDensity
kp <- kpPlotDensity(kp, data.panel=1, data=snps_ddRAD, r0=0.7, r1=0.95, window.size=40000, border= "#DC3220", col= "#DC3220")
kpAxis(kp, data.panel=1, chr="chr13", ymax=kp$latest.plot$computed.values$max.density, r0=0.7, r1=0.95, cex=1.8)


#plot repeated regions as regions
repeats <- toGRanges("chr_13_masked_3000bp.bed")
kpPlotRegions(kp, data.panel=2, r0=autotrack(1,3), data=repeats, col="darkviolet")

#plot GC content as heatmaps
gc.cont <- toGRanges("chr_13_GC.txt")
kpHeatmap(kp, data.panel=2, r0=autotrack(2,3), data=gc.cont, y=gc.cont$V4, col=c("black","orange","#ed21ff"))

#plot Pacbiocoverage as heatmap using the viridis palette
Pacbio_cov <- toGRanges("chr_13_genomecov_500bp")
kpHeatmap(kp, data.panel=2, data=Pacbio_cov, y=Pacbio_cov$V4, r0=autotrack(3,3), col=viridis(5))

#plot gaps on ideogram as regions
gaps<-toGRanges("chr_13_gaps.bed")
kpPlotRegions(kp, data.panel="ideogram", data=gaps, col="black")

