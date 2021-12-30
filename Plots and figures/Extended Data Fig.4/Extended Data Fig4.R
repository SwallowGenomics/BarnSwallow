setwd("E:/PhD/Genomica/PAPER BARN SWALLOW GENOMICS/CONFRONTO/CACTUS/PROVA3 - AURELIANO/INTERSECT")

#load phyloP files without GAPS and with ID
acc_CM <- read.table("PhyloP_10bp_ACC_NO_GAPS_ID.bed.bed")
cons_CM <- read.table("PhyloP_10bp_CONS_NO_GAPS_ID.bed.bed")
#renmae columns
colnames(acc_CM) <- c("chr", "start", "end", "logP", "id")
colnames(cons_CM) <- c("chr", "start", "end", "logP", "id") 

#change columns order
acc_CM1 <- data.frame(acc_CM$id, acc_CM$chr, acc_CM$end, acc_CM$logP)
cons_CM1 <- data.frame(cons_CM$id, cons_CM$chr, cons_CM$end, cons_CM$logP)
#rename columns
colnames(acc_CM1) <- c("id", "chr", "end", "logP")
colnames(cons_CM1) <- c("id", "chr", "end", "logP") 

#absolute values on accelerated sites (they are negative)
acc_CM1_abs <- data.frame(acc_CM1$id, acc_CM1$chr, acc_CM1$end, abs(acc_CM1$logP) )
colnames(acc_CM1_abs) <- c("id", "chr", "end", "logP") 
range(acc_CM1_abs$logP ) #2.698 20.000

#Bonferroni tresholds (total sites tested, both accelerated and conserved: 98820052)
-log((0.05/98820052),10) #9.295875 logP
10^-(9.295875) #5.059703e-10 Pvalue

library("CMplot")

col_all = c("#d10000", "#d21c00", "#d33500", "#d55200", "#d66b00", "#d78800", "#d8a200", "#d9c000", "#dada00", "#bedb00", "#a6dd00", "#89de00", "#70df00", "#52e000", "#38e100", "#1ae200", "#00e300", "#00e51f", "#00e63a", "#00e759", "#00e874", "#00e994", "#00eab0", "#00ebd0", "#00eded", "#00ceee", "#00b3ef", "#0094f0", "#0079f1", "#0059f2", "#003df4", "#001df5", "#0000f6", "#2100f7", "#3e00f8", "#5f00f9", "#7d00fa", "#a000fc", "#be00fd", "#e000fe", "#ff00ff")

#panel a
CMplot(acc_CM1_abs, type = "p", plot.type = "m", LOG10 = FALSE, col=col_all, lwd.axis=3, cex.axis=4,cex.lab=5,ylab.pos=4.5,xticks.pos=3,chr.labels.angle=45, amplify=TRUE, threshold=9.295875, threshold.col="black", threshold.lwd=4, threshold.lty=1, band=0.5, memo= "panelA", dpi=600,  file = "jpg", file.output=TRUE, verbose = TRUE, width = 70, height = 20, mar=c(10,10,6,6))

#panel b
CMplot(cons_CM1, type = "p", plot.type = "m", ylim=c(2.69,3.6), LOG10 = FALSE, col=col_all, lwd.axis=3, cex.axis=4,cex.lab=5,ylab.pos=6,xticks.pos=3,chr.labels.angle=45, amplify=TRUE, band=0.5, memo= "panelB", dpi=600,  file = "jpg", file.output=TRUE, verbose = TRUE, width = 70, height = 20, mar=c(10,11,6,6))

#panel c
pdf("panelC.pdf", width = 9, height = 5)
hist(acc_CM1$logP, main=NA, breaks=100, xlim=c(-20,-2), xlab="Accelerated logP", col="#a6cee3", border=NA)
dev.off()

#panel d
pdf("panelD.pdf", width = 9, height = 5)
hist(cons_CM1$logP, main=NA,breaks=100,xlim =c(2.6,3.6), xlab="Conserved logP", col="#1f78b4", border=NA)
dev.off()

#panel e 
perc <- read.table("percentage_table.txt", header=TRUE)

library(reshape2)

perc_melt <- melt(perc)

library(ggplot2)

pdf("panelE.pdf", width = 7, height = 5)

ggplot(data=perc_melt, aes(x=Feature, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=value), vjust=0.6,hjust=-0.2, color="black",position = position_dodge(0.9), size=2.5)+
  #geom_text(aes(label=value), vjust=-0.5, color="black",position = position_dodge(0.9), size=3)+
  scale_x_discrete(limits = unique(perc_melt$Feature)) +
  scale_fill_brewer(palette="Paired", direction = -1) +
  labs(x ="Functional feature fraction              Site functional category", y = "Bases (%)", size =1) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  theme_classic() +
  coord_flip()

dev.off()
