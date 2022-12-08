setwd("E:/path/to/folder/")
chrs = read.table("CHR_coords.bed", colClasses = c("character", "numeric", "numeric"), sep = "\t")
library(circlize)

#pangenome colors as panels a and b
col_text <- "grey40"
col_track1 <- "grey40"
col_bHirRus1p <- "#000000"
col_bHirRus1a <- "#edb497"
col_Hr2p <- "#ff7f0e"
col_Hr2a <- "#1f77b4"
col_Hr3p <- "#d62728"
col_Hr3a <- "#2ca02c"
col_Hr4p <- "#8c564b"
col_Hr4a <- "#9467bd"
col_HrA1p <- "#7f7f7f"
col_HrA1a <- "#e377c2"
col_HrA2p <- "#17becf"
col_HrA2a <- "#f9cf00"

pdf("FIGURE5_panelA.pdf", width = 6, height = 6)

circos.clear()
circos.par(points.overflow.warning=FALSE, "start.degree" = 90, canvas.ylim = c(-1.5,1), canvas.xlim = c(-1.1,1.1), "cell.padding" = c(0,0,0,0), "gap.degree" = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,2,2,2))
circos.initialize(factors=chrs$V1, xlim=matrix(c(chrs$V2,chrs$V3),ncol=2))

#generate legend
library(ComplexHeatmap)
lgd = Legend(labels = c("bHirRus1p","bHirRus1a","Hr2p","Hr2a","Hr3p","Hr3a","Hr4p","Hr4a","HrA1p","HrA1a","HrA2p","HrA2a"),
             labels_gp = gpar(fontsize=6), ncol=6,legend_gp = gpar(fill = c(col_bHirRus1p,col_bHirRus1a,col_Hr2p,col_Hr2a,col_Hr3p,col_Hr3a,col_Hr4p,col_Hr4a,col_HrA1p,col_HrA1a,col_HrA2p,col_HrA2a, lwd = 2)),
             grid_height = unit(2, "mm"), grid_width = unit(3, "mm"))
draw(lgd, x = unit(3, "cm"), y = unit(2, "cm"), just = c("left"))

#circos.par(points.overflow.warning=FALSE, "start.degree" = 90, canvas.ylim = c(-1.5,1), canvas.xlim = c(-1.1,1.1), "cell.padding" = c(0,0,0,0), "gap.degree" = c(2))
#circos.initialize(factors=chrs$V1, xlim=matrix(c(chrs$V2,chrs$V3),ncol=2))

#TRACK 1 - CHROMOSOMES
circos.track(ylim=c(0,1),panel.fun=function(x,y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  length= (max(xlim))
  length_label = paste((format(round(length/(10^6), 1), nsmall = 1)), "Mb", sep="")
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, font = c(2), cex=0.5, col=col_track1, facing = "down")
},bg.col="white",bg.border=NA,track.height=0.05)

#TRACK 2 - bHirRu1 genes

bHirRus1 <- read.table("genes_PANGENOME_4_col_CIRLIZE.bed")
colnames(bHirRus1) <- c("chr", "start", "end", "gene")

circos.genomicTrackPlotRegion(bHirRus1, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_bHirRus1p, border = NA, ...)
                              })

#TRACK 3 - Hr2 pri genes

Hr2p <- read.table("Hr2p_orthologFind.bed")
colnames(Hr2p) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(Hr2p, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_Hr2p, border = NA, ...)
                              })
#TRACK 4 - Hr2 alt genes

Hr2a <- read.table("Hr2a_orthologFind.bed")
colnames(Hr2a) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(Hr2a, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_Hr2a, border = NA, ...)
                              })

#TRACK 5 - Hr3 pri genes

Hr3p <- read.table("Hr3p_orthologFind.bed")
colnames(Hr3p) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(Hr3p, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_Hr3p, border = NA, ...)
                              })
#TRACK 6 - Hr3 alt genes

Hr3a <- read.table("Hr3a_orthologFind.bed")
colnames(Hr3a) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(Hr3a, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_Hr3a, border = NA, ...)
                              })

#TRACK 7 - Hr4 pri genes

Hr4p <- read.table("Hr4p_orthologFind.bed")
colnames(Hr4p) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(Hr4p, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_Hr4p, border = NA, ...)
                              })
#TRACK 8 - Hr4 alt genes

Hr4a <- read.table("Hr4a_orthologFind.bed")
colnames(Hr4a) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(Hr4a, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_Hr4a, border = NA, ...)
                              })

#TRACK 9 - HrA1 pri genes

HrA1p <- read.table("HrA1p_orthologFind.bed")
colnames(HrA1p) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(HrA1p, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_HrA1p, border = NA, ...)
                              })
#TRACK 10 - Hr4 alt genes

HrA1a <- read.table("HrA1a_orthologFind.bed")
colnames(HrA1a) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(HrA1a, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_HrA1a, border = NA, ...)
                              })

#TRACK 11 - HrA2 pri genes

HrA2p <- read.table("HrA2p_orthologFind.bed")
colnames(HrA2p) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(HrA2p, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_HrA2p, border = NA, ...)
                              })
#TRACK 12 - HrA2 alt genes

HrA2a <- read.table("HrA2a_orthologFind.bed")
colnames(HrA2a) <- c("chr", "start", "end", "genes")

circos.genomicTrackPlotRegion(HrA2a, track.height = 0.05, bg.border = NA, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = col_HrA2a, border = NA, ...)
                              })

dev.off()
