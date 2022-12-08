setwd("path/to/folder")

#load chromosome coordinates
chrs = read.table("CHR_coords.bed", colClasses = c("character", "numeric", "numeric"), sep = "\t")

library(circlize)

#set text and tracks color
col_text <- "grey40"
col_track1 <- "grey40"
col_track2 <- "#ff0900"
col_track3 <- "#ff8400"
col_track4 <- "#ffd000"
col_track5 <- "#86eb03"
col_track6 <- "#01f8b5"
col_track7 <- "#01d8ff"
col_track8 <- "#0066ff"
col_track9 <- "#9421ff"
col_track10 <- "#ed21ff"

pdf("FIGURE1_panelC.pdf", width = 6, height = 6) 

circos.clear()
circos.par(points.overflow.warning=FALSE, "start.degree" = 90, canvas.ylim = c(-1,1), canvas.xlim = c(-1.1,1.1), "cell.padding" = c(0,0,0,0), "gap.degree" = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 11))
circos.initialize(factors=chrs$V1, xlim=matrix(c(chrs$V2,chrs$V3),ncol=2))

#TRACK 1 - CHROMOSOMES
circos.track(ylim=c(0,1),panel.fun=function(x,y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  length= (max(xlim))
  length_label = paste((format(round(length/(10^6), 1), nsmall = 1)), "Mb", sep="")
  ylim=CELL_META$ylim
  
  if(length > 63258489) {
    ticks = c(0, 20*10^6, 40*10^6, 60*10^6, 80*10^6, 100*10^6, 120*10^6, 140*10^6, 160*10^6, length)
    labels =  c("0Mb", "20Mb", "40Mb", "60Mb", "80Mb", "100Mb", "120Mb", "140Mb", "160Mb", length_label)
    circos.genomicAxis(h = "top", lwd = 0.2, major.at = ticks, labels = labels, 
                       labels.cex = 0.5, col=col_text, labels.col=col_text,
                       labels.facing="clockwise",
                       major.tick.length = convert_y(0.7, unit=c("mm")))
  } else if (length > 25880253){
    ticks = c(0, 20*10^6, 40*10^6, length)
    labels =  c("0Mb", "20Mb", "40Mb", length_label)
    circos.genomicAxis(h = "top", lwd = 0.2, major.at = ticks, labels = labels, 
                       labels.cex = 0.5, col=col_text, labels.col=col_text,
                       labels.facing="clockwise",
                       major.tick.length = convert_y(0.7, unit=c("mm")))
  }else if (length > 16541138){
    ticks = c(0, length)
    labels = c("0Mb", length_label)
    circos.genomicAxis(h = "top", lwd = 0.2, major.at = ticks, labels = labels, 
                       labels.cex = 0.5, col=col_text, labels.col=col_text,
                       labels.facing="clockwise",
                       major.tick.length = convert_y(0.7, unit=c("mm")))
  } else if (length > 2102120){
    ticks = c(0, length)
    labels = c("", length_label)
    circos.genomicAxis(h = "top", lwd = 0.2, major.at = ticks, labels = labels, 
                       labels.cex = 0.5, col=col_text, labels.col=col_text,
                       labels.facing="clockwise",
                       major.tick.length = convert_y(0.7, unit=c("mm")))
  } else {
    ticks = c(length)
    labels = c(length_label)
    circos.genomicAxis(h = "top", lwd = 0.2, major.at = ticks, labels = labels, 
                       labels.cex = 0.5, col=col_text,labels.col=col_text,
                       labels.facing="clockwise", minor.ticks=0, 
                       major.tick.length = convert_y(0.7, unit=c("mm")))}
  
  circos.text(mean(xlim), mean(ylim), chr, font = c(2), cex=0.48, col=col_track1, facing = "down")
  
},bg.col="grey90",bg.border=NA,track.height=0.05)

#LEGEND
library(ComplexHeatmap)
lgd = Legend(labels = c("PacBio cov.", "Repeats", "GC content", "CpG islands", "Genes", "PhyloP acc.", "PhyloP cons.", "PhastCons CEs", "HAL cov."),
             labels_gp = gpar(fontsize=6), legend_gp = gpar(fill = c(col_track2,col_track3,col_track4,col_track5,col_track6, col_track7, col_track8, col_track9, col_track10, lwd = 2)),
             grid_height = unit(3, "mm"), grid_width = unit(3, "mm"))

draw(lgd, x = unit(0.66, "cm"), y = unit(2.2, "cm"), just = c("left"))

#TRACK 2 - PACBIO READS COVERAGE

#load PacBio coverage file generated with Mosdepth on 200 kbp windows
PacBio <- read.table("PacBio_cov_200k.bed")
colnames(PacBio) <- c("chr", "start", "end", "cov")

#calculate 99th percentile
quantile(PacBio$cov, c(.99)) #88.6438 

#cap the values > 100 
library(tidyverse)
PacBio_sub <- PacBio %>% mutate(cov = ifelse(cov > 100, 100, cov)) 

#add track
circos.genomicTrackPlotRegion(PacBio_sub, track.height = 0.06, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, border = col_track2, col = col_track2, ...)
}, ylim = range(PacBio_sub$cov), bg.border = NA)

#add y axis
circos.yaxis(at=c(0, (max(PacBio_sub$cov)/2), max(PacBio_sub$cov)), labels = c(0,50,"100+"), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))

#TRACK 3 - REPEATS DENSITY

#in percentage

#load file with intersection between repeats bases and 200 kbp windows (bedtools intersect) 
repeat_density_p <- read.table("REPEATS_200k.bed")
colnames(repeat_density_p) <- c("chr", "start", "end", "ovl")

#aggregate windows when chr, start, ends coincide
repeat_desity_p_agg <- aggregate(.~chr+start+end, repeat_density_p, sum)

#calculate fraction of the windows covered by repeats
repeat_desity_p_agg_percent <- data.frame(repeat_desity_p_agg$chr, repeat_desity_p_agg$start, repeat_desity_p_agg$end, (repeat_desity_p_agg$ovl/(repeat_desity_p_agg$end-repeat_desity_p_agg$start))*100)
colnames(repeat_desity_p_agg_percent) <- c("chr", "start", "end", "repeats")

#add track
circos.genomicTrackPlotRegion(repeat_desity_p_agg_percent, track.height = 0.06, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, border= col_track3, col = col_track3,  ...)
}, ylim = range(repeat_desity_p_agg_percent$repeats), bg.border = NA)

#add y axis
circos.yaxis(at=c(0, (max(repeat_desity_p_agg_percent$repeats)/2), max(repeat_desity_p_agg_percent$repeats)), labels = c(0,"44","88"), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))

#TRACK 4 - GC CONTENT

#in percentage

#load GC content calculated with bedtools nuc on 200 kbp windows
GC <- read.table("GC_CONTENT_200k.bed")
colnames(GC) <- c("chr", "start", "end", "GC") 

#add track
circos.genomicTrackPlotRegion(GC, track.height = 0.06, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, border= col_track4, col = col_track4, ...)
}, ylim = range(GC$GC), bg.border = NA)

#add y axis
circos.yaxis(at=c(0, (max(GC$GC)/2), max(GC$GC)), labels = c(0,"33.2","66.3"), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))

#TRACK 5 - CpG ISLANDS DENSITY

#number

#load file with intersection between CpG islands bases and 200 kbp windows (bedtools intersect) 
CpG_density_n <- read.table("CpG_ISLANDS_200k.bed")
colnames(CpG_density_n) <- c("chr", "start", "end", "CpG") 

#calculate 99th percentile
quantile(CpG_density_n$CpG, c(.99)) #42

#cap values > 42
library(tidyverse)
CpG_density_n_sub <- CpG_density_n %>% mutate(CpG = ifelse(CpG > 42, 42, CpG)) #uso 42 come limite

#add track
circos.genomicTrackPlotRegion(CpG_density_n_sub, track.height = 0.06, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, col = col_track5, border = col_track5, ...)
}, ylim = range(CpG_density_n_sub$CpG), bg.border = NA)

#add y axis
circos.yaxis(at=c(0, (max(CpG_density_n_sub$CpG)/2), max(CpG_density_n_sub$CpG)), labels = c(0,"21","42+"), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))

#TRACK 6 - GENE DENSITY

#number

#load file with intersection between genes bases and 200 kbp windows (bedtools intersect) 
gene_density_n <- read.table("genes_CHRs_merged_200k_n_3839.bed")
colnames(gene_density_n) <- c("chr", "start", "end", "genes") 

#calculate 99th percentile
quantile(gene_density_n$genes, c(.99)) #16

#cap values > 16
library(tidyverse)
gene_density_n_sub <- gene_density_n %>% mutate(genes = ifelse(genes > 16, 16, genes))

#add track
circos.genomicTrackPlotRegion(gene_density_n_sub, track.height = 0.06, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, border= col_track6, col = col_track6,  ...)
}, ylim = range(gene_density_n_sub$genes), bg.border = NA)

#add y axis
circos.yaxis(at=c(0, (max(gene_density_n_sub$genes)/2), max(gene_density_n_sub$genes)), labels = c(0,"8","16+"), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))


#TRACK 7 - ACCELERATED SITES 

#number 

#load file with intersection between phyloP accelerated sites and 200 kbp windows (bedtools intersect) 
acc_density_n <- read.table("PhyloP_10bp_FDR_ACC_200k.bed")
colnames(acc_density_n) <- c("chr", "start", "end", "acc")

#calculate 99th percentile
quantile(acc_density_n$acc, c(.99)) #813.66 

#cap values > 814
library(tidyverse)
acc_density_n_sub <- acc_density_n %>% mutate(acc = ifelse(acc > 814, 814, acc))

#add track
circos.genomicTrackPlotRegion(acc_density_n_sub, track.height = 0.06, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, border= col_track7, col = col_track7,  ...)
}, ylim = range(acc_density_n_sub$acc), bg.border = NA)

#add y axis
circos.yaxis(at=c(0, (max(acc_density_n_sub$acc)/2), max(acc_density_n_sub$acc)), labels = c("0","407","814+"), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))

#TRACK 8 - CONSERVED PHYLOP SITES 

#number

#load file with intersection between phyloP conserved sites and 200 kbp windows (bedtools intersect) 
cons_density_n <- read.table("PhyloP_10bp_FDR_CONS_200k.bed")
colnames(cons_density_n) <- c("chr", "start", "end", "cons")

#calculate 99th percentile
quantile(cons_density_n$cons, c(.99)) #2106.82  

#cap values > 2106
library(tidyverse)
cons_density_n_sub <- cons_density_n %>% mutate(cons = ifelse(cons > 2106, 2106, cons))

#add track
circos.genomicTrackPlotRegion(cons_density_n_sub, track.height = 0.06, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, border= col_track8, col = col_track8,  ...)
}, ylim = range(cons_density_n_sub$cons), bg.border = NA)

#add y axis
circos.yaxis(at=c(0, (max(cons_density_n_sub$cons)/2), max(cons_density_n_sub$cons)), labels = c("0","1053","2106+"), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))

#TRACK 9 - PhastCons CEs

#number 

#load file with intersection between phastCons CEs sites and 200 kbp windows (bedtools intersect) 
most_cons_density_n <- read.table("phastCons_CEs_200k.bed")
colnames(most_cons_density_n) <- c("chr", "start", "end", "mostcons")

#add track
circos.genomicTrackPlotRegion(most_cons_density_n, track.height = 0.06, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, border= col_track9, col = col_track9,  ...)
}, ylim = range(most_cons_density_n$mostcons), bg.border = NA)

#add y axis
circos.yaxis(at=c(0, (max(most_cons_density_n$mostcons)/2), max(most_cons_density_n$mostcons)), labels = c("0","619","1238"), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))


#TRACK 10 - HAL COVERAGE

#number

#load file with number of species aligned to 200 kbp windows (halCoverage)
HAL <- read.table("HAL_cov_200k.bed")
colnames(HAL) <- c("chr", "start", "end", "id", "cov")
HAL$end <- HAL$end+199999 #modify end

#add track
circos.genomicTrackPlotRegion(HAL, track.height = 0.05, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", lwd = 0.1, area=TRUE, border= col_track10, col = col_track10, ...)
}, ylim = range(HAL$cov), bg.border = NA)

#add y axis
circos.yaxis(at=c(min(HAL$cov), (max(HAL$cov)/2), max(HAL$cov)), labels = c(0,3.5,7), sector.index = "1", track.index = get.current.track.index(), 
             labels.cex=0.31, lwd=0, labels.col=col_text, col=col_text, tick = TRUE, 
             tick.length = convert_x(0.4, "mm", sector.index = get.current.sector.index(), track.index=get.current.track.index()))

dev.off()