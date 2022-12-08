setwd("/path/to/folder/")
library(circlize)

col_text <- "grey40"
col_track1 <- "grey40"

name <- "mito"
start <- as.numeric("0")
end <- as.numeric("16277") #mitogenome length
size <- data.frame(name, start, end)

pdf("Figure S1.pdf", width = 15, height = 15) 

#chromosome track
circos.clear()
circos.par(points.overflow.warning=FALSE, gap.degree=0, "start.degree" = 90, "canvas.ylim" = c(0, 0), "canvas.xlim" = c(-1.3, 1.3), cell.padding = c(0, 0, 0, 0))
circos.initialize(factors=size$name, xlim=matrix(c(size$start,size$end),ncol=2))
#canvas.ylim = c(-1.5,1)
#canvas.xlim = c(-1.1,1.1)

circos.track(ylim=c(0,1),panel.fun=function(x,y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  length= (max(xlim))
  length_label = paste((format(round(length/(10^3), 1), nsmall = 1)), "kb", sep="")
  ylim=CELL_META$ylim
  
  ticks = c(0, 1*10^3, 2*10^3, 3*10^3, 4*10^3, 5*10^3, 6*10^3, 7*10^3, 8*10^3, 9*10^3, 10*10^3, 11*10^3, 12*10^3, 13*10^3, 14*10^3, 15*10^3,length)
  labels =  c("0kb", "1kb", "2kb", "3kb", "4kb", "5kb", "6kb", "7kb", "8kb", "9kb", "10kb", "11kb", "12kb", "13kb", "14kb", "15kb", length_label)
  circos.genomicAxis(h = "top", lwd = 1.7, major.at = ticks, labels = labels, 
                     labels.cex = 2, col=col_text, labels.col=col_text,
                     labels.facing="clockwise",
                     major.tick.length = convert_y(2, unit=c("mm")))
},bg.col="grey90",bg.border=NA,track.height=0.05)

circos.clear()

#first circle
bHirRus1 <- read.table("Mitos2.txt", sep="\t", header = FALSE)
names(bHirRus1) <- c("dataset", "type","name","start","end","a","b")
sectors = cbind(bHirRus1$start, bHirRus1$end)

bHirRus1$type=as.factor(bHirRus1$type)
color_easy = c("#FBBC05", "#4285F4", "#EA4335", "#34A853")[bHirRus1$type]

par(new = TRUE) # <- magic
circos.par(gap.degree=0, "canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(0, 0), "start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize(bHirRus1$name, xlim = sectors)

circos.track(ylim = c(1, 1.4), panel.fun = function(x, y, ...) {
  circos.arrow(CELL_META$xlim[1], CELL_META$xlim[2], 
               arrow.head.width = CELL_META$yrange*0.6, arrow.head.length = ux(0.31, "cm"),
               col = color_easy[CELL_META$sector.numeric.index])
}, bg.border = NA, "track.height" = 0.3)

dev.off()
