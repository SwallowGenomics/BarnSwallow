#load chr. sizes and type file
sizes = read.table("CHR_sizes_type.txt", colClasses = c("character", "numeric", "character"), sep = "\t")
sizes <- sizes[order(sizes$V1),]
colnames(sizes) <- c("chr", "size", "type")

#panel a

GC <- read.table("GC_CONTENT.bed")
colnames(GC) <- c("chr", "size", "GC", "type") 

library(dplyr)
GC_scatter <- GC %>% mutate(type = case_when(chr == 1 ~ "Macro.", chr == 2 ~ "Macro.", chr == 3 ~ "Macro.", chr == 4 ~ "Macro.", chr == 5 ~ "Macro.", chr == 6 ~ "Macro.", chr == 7 ~ "Macro.", chr == 8 ~ "Macro.", chr == 9 ~ "Inter.", chr == 10 ~ "Inter.", chr == 11 ~ "Inter.", chr == 12 ~ "Inter.", chr == 13 ~ "Inter.", chr == 14 ~ "Micro.", chr == 15 ~ "Micro.", chr == 16 ~ "Micro.", chr == 17 ~ "Micro.", chr == 18 ~ "Micro.", chr == 19 ~ "Micro.", chr == 20 ~ "Micro.", chr == 21 ~ "Micro.", chr == 22 ~ "Micro.", chr == 23 ~ "Micro.", chr == 24 ~ "Micro.", chr == 25 ~ "Micro.", chr == 26 ~ "Micro.", chr == 27 ~ "Micro.", chr == 28 ~ "Micro.", chr == 29 ~ "Micro.", chr == 30 ~ "Micro.", chr == 31 ~ "Excluded", chr == 32 ~ "Micro.", chr == 33 ~ "Excluded", chr == 34 ~ "Excluded", chr == 35 ~ "Micro.", chr == 36 ~ "Micro.", chr == 37 ~ "Micro.", chr == 38 ~ "Micro.", chr == 39 ~ "Micro.", chr == "W" ~ "Inter.", chr == "Z" ~ "Macro."))

library(ggrepel)
library(ggplot2)
pdf("panelA.pdf", width = 8, height = 5)
ggplot(GC_scatter, aes(log(size), GC, color=type, label=chr)) +
  geom_point() +
  scale_color_manual(values=c("grey", "#ff006d", "#ffb500", "#007eff")) +
  geom_text_repel(min.segment.length=0.5,box.padding=0.3,size=4, nudge_x = 0, nudge_y = 0.1 ) +
  theme_light() +
  labs (y= "GC",  x="Chr size (log)") 
dev.off()

#panel b

CpG <- read.table("CpG_ISLANDS.bed")
colnames(CpG) <- c("chr", "start", "end") 
#number bases per chr
CpG$length <- CpG$end-CpG$start
CpG_agg <- aggregate(CpG[, 4], list(CpG$chr), sum)
colnames(CpG_agg) <- c("chr", "sumCpG")
#add chr. size
CpG1 <- data.frame(CpG_agg$chr, sizes$size, CpG_agg$sumCpG)
colnames(CpG1) <- c("chr", "size", "sumCpG")
CpG1$percent <- (CpG1$sumCpG/CpG1$size)*100

#add type
CpG1 <- CpG1 %>% mutate(type = case_when(chr == 1 ~ "Macro.", chr == 2 ~ "Macro.", chr == 3 ~ "Macro.", chr == 4 ~ "Macro.", chr == 5 ~ "Macro.", chr == 6 ~ "Macro.", chr == 7 ~ "Inter.", chr == 8 ~ "Inter.", chr == 9 ~ "Inter.", chr == 10 ~ "Inter.", chr == 11 ~ "Inter.", chr == 12 ~ "Inter.", chr == 13 ~ "Inter.", chr == 14 ~ "Micro.", chr == 15 ~ "Micro.", chr == 16 ~ "Micro.", chr == 17 ~ "Micro.", chr == 18 ~ "Micro.", chr == 19 ~ "Micro.", chr == 20 ~ "Micro.", chr == 21 ~ "Micro.", chr == 22 ~ "Micro.", chr == 23 ~ "Micro.", chr == 24 ~ "Micro.", chr == 25 ~ "Micro.", chr == 26 ~ "Micro.", chr == 27 ~ "Micro.", chr == 28 ~ "Micro.", chr == 29 ~ "Micro.", chr == 30 ~ "Micro.", chr == 31 ~ "Micro.", chr == 32 ~ "Micro.", chr == 33 ~ "Micro.", chr == 34 ~ "Micro.", chr == 35 ~ "Micro.", chr == 36 ~ "Micro.", chr == 37 ~ "Micro.",  chr == 38 ~ "Micro.",  chr == 39 ~ "Micro.", chr == "W" ~ "Inter.", chr == "Z" ~ "Macro."))

CpG_scatter <- CpG1 %>% mutate(type = case_when(chr == 1 ~ "Macro.", chr == 2 ~ "Macro.", chr == 3 ~ "Macro.", chr == 4 ~ "Macro.", chr == 5 ~ "Macro.", chr == 6 ~ "Macro.", chr == 7 ~ "Inter.", chr == 8 ~ "Inter.", chr == 9 ~ "Inter.", chr == 10 ~ "Inter.", chr == 11 ~ "Inter.", chr == 12 ~ "Inter.", chr == 13 ~ "Inter.", chr == 14 ~ "Micro.", chr == 15 ~ "Micro.", chr == 16 ~ "Micro.", chr == 17 ~ "Micro.", chr == 18 ~ "Micro.", chr == 19 ~ "Micro.", chr == 20 ~ "Micro.", chr == 21 ~ "Micro.", chr == 22 ~ "Micro.", chr == 23 ~ "Micro.", chr == 24 ~ "Micro.", chr == 25 ~ "Micro.", chr == 26 ~ "Micro.", chr == 27 ~ "Micro.", chr == 28 ~ "Micro.", chr == 29 ~ "Micro.", chr == 30 ~ "Micro.", chr == 31 ~ "Excluded", chr == 32 ~ "Micro.", chr == 33 ~ "Excluded", chr == 34 ~ "Excluded", chr == 35 ~ "Micro.", chr == 36 ~ "Micro.", chr == 37 ~ "Micro.", chr == 38 ~ "Micro.", chr == 39 ~ "Micro.", chr == "W" ~ "Inter.", chr == "Z" ~ "Macro."))

pdf("panelB.pdf", width = 8, height = 5)
ggplot(CpG_scatter, aes(log(size), percent, color=type, label=chr)) +
  geom_point() +
  scale_color_manual(values=c("grey", "#ff006d", "#ffb500", "#007eff")) +
  geom_text_repel(min.segment.length=0.5,box.padding=0.3,size=4, nudge_x = 0, nudge_y = 0.1 ) +
  theme_light() +
  labs (y= "CpG islands",  x="Chr size (log)") 
dev.off()

#panel c

genes <- read.table("GENES.bed")
colnames(genes) <- c("chr", "start", "end") 
#number bases per chr
genes$length <- genes$end-genes$start
genes_agg <- aggregate(genes[, 4], list(genes$chr), sum)
colnames(genes_agg) <- c("chr", "sumgenes")
#add chr. size
genes1 <- data.frame(genes_agg$chr, sizes$size, genes_agg$sumgenes)
colnames(genes1) <- c("chr", "size", "sumgenes")
genes1$percent <- (genes1$sumgenes/genes1$size)*100

#add chr. type
genes1 <- genes1 %>% mutate(type = case_when(chr == 1 ~ "Macro.", chr == 2 ~ "Macro.", chr == 3 ~ "Macro.", chr == 4 ~ "Macro.", chr == 5 ~ "Macro.", chr == 6 ~ "Macro.", chr == 7 ~ "Inter.", chr == 8 ~ "Inter.", chr == 9 ~ "Inter.", chr == 10 ~ "Inter.", chr == 11 ~ "Inter.", chr == 12 ~ "Inter.", chr == 13 ~ "Inter.", chr == 14 ~ "Micro.", chr == 15 ~ "Micro.", chr == 16 ~ "Micro.", chr == 17 ~ "Micro.", chr == 18 ~ "Micro.", chr == 19 ~ "Micro.", chr == 20 ~ "Micro.", chr == 21 ~ "Micro.", chr == 22 ~ "Micro.", chr == 23 ~ "Micro.", chr == 24 ~ "Micro.", chr == 25 ~ "Micro.", chr == 26 ~ "Micro.", chr == 27 ~ "Micro.", chr == 28 ~ "Micro.", chr == 29 ~ "Micro.", chr == 30 ~ "Micro.", chr == 31 ~ "Micro.", chr == 32 ~ "Micro.", chr == 33 ~ "Micro.", chr == 34 ~ "Micro.", chr == 35 ~ "Micro.", chr == 36 ~ "Micro.", chr == 37 ~ "Micro.", chr == 38 ~ "Micro.", chr == 39 ~ "Micro.",chr == "W" ~ "Inter.", chr == "Z" ~ "Macro."))

genes_scatter <- genes1 %>% mutate(type = case_when(chr == 1 ~ "Macro.", chr == 2 ~ "Macro.", chr == 3 ~ "Macro.", chr == 4 ~ "Macro.", chr == 5 ~ "Macro.", chr == 6 ~ "Macro.", chr == 7 ~ "Inter.", chr == 8 ~ "Inter.", chr == 9 ~ "Inter.", chr == 10 ~ "Inter.", chr == 11 ~ "Inter.", chr == 12 ~ "Inter.", chr == 13 ~ "Inter.", chr == 14 ~ "Micro.", chr == 15 ~ "Micro.", chr == 16 ~ "Micro.", chr == 17 ~ "Micro.", chr == 18 ~ "Micro.", chr == 19 ~ "Micro.", chr == 20 ~ "Micro.", chr == 21 ~ "Micro.", chr == 22 ~ "Micro.", chr == 23 ~ "Micro.", chr == 24 ~ "Micro.", chr == 25 ~ "Micro.", chr == 26 ~ "Micro.", chr == 27 ~ "Micro.", chr == 28 ~ "Micro.", chr == 29 ~ "Micro.", chr == 30 ~ "Micro.", chr == 31 ~ "Excluded", chr == 32 ~ "Micro.", chr == 33 ~ "Excluded", chr == 34 ~ "Excluded", chr == 35 ~ "Micro.", chr == 36 ~ "Micro.", chr == 37 ~ "Micro.", chr == 38 ~ "Micro.", chr == 39 ~ "Micro.", chr == "W" ~ "Inter.", chr == "Z" ~ "Macro."))

pdf("panelC.pdf", width = 8, height = 5)
ggplot(genes_scatter, aes(log(size), percent, color=type, label=chr)) +
  geom_point() +
  scale_color_manual(values=c("grey", "#ff006d", "#ffb500", "#007eff")) +
  geom_text_repel(min.segment.length=0.5,box.padding=0.3,size=4, nudge_x = 0, nudge_y = 0.1 ) +
  theme_light() +
  labs (y= "Genes",  x="Chr size (log)") 
dev.off()
#panel d

repeats <- read.table("REPEATS.bed")

colnames(repeats) <- c("chr", "start", "end")

#number bases per chr
repeats$length <- repeats$end-repeats$start

repeats_agg <- aggregate(repeats[, 4], list(repeats$chr), sum)

colnames(repeats_agg) <- c("chr", "sumrepeats")

#aggiungo chr. sizes
repeats1 <- data.frame(repeats_agg$chr, sizes$size, repeats_agg$sumrepeats)

colnames(repeats1) <- c("chr", "size", "sumrepeats")

repeats1$percent <- (repeats1$sumrepeats/repeats1$size)*100

#add chr. type
repeats1 <- repeats1 %>% mutate(type = case_when(chr == 1 ~ "Macro.", chr == 2 ~ "Macro.", chr == 3 ~ "Macro.", chr == 4 ~ "Macro.", chr == 5 ~ "Macro.", chr == 6 ~ "Macro.", chr == 7 ~ "Inter.", chr == 8 ~ "Inter.", chr == 9 ~ "Inter.", chr == 10 ~ "Inter.", chr == 11 ~ "Inter.", chr == 12 ~ "Inter.", chr == 13 ~ "Inter.", chr == 14 ~ "Micro.", chr == 15 ~ "Micro.", chr == 16 ~ "Micro.", chr == 17 ~ "Micro.", chr == 18 ~ "Micro.", chr == 19 ~ "Micro.", chr == 20 ~ "Micro.", chr == 21 ~ "Micro.", chr == 22 ~ "Micro.", chr == 23 ~ "Micro.", chr == 24 ~ "Micro.", chr == 25 ~ "Micro.", chr == 26 ~ "Micro.", chr == 27 ~ "Micro.", chr == 28 ~ "Micro.", chr == 29 ~ "Micro.", chr == 30 ~ "Micro.", chr == 31 ~ "Micro.", chr == 32 ~ "Micro.", chr == 33 ~ "Micro.", chr == 34 ~ "Micro.", chr == 35 ~ "Micro.", chr == 36 ~ "Micro.", chr == 37 ~ "Micro.", chr == 38 ~ "Micro.", chr == 39 ~ "Micro.", chr == "W" ~ "Inter.", chr == "Z" ~ "Macro."))

repeats_scatter <- repeats1 %>% mutate(type = case_when(chr == 1 ~ "Macro.", chr == 2 ~ "Macro.", chr == 3 ~ "Macro.", chr == 4 ~ "Macro.", chr == 5 ~ "Macro.", chr == 6 ~ "Macro.", chr == 7 ~ "Inter.", chr == 8 ~ "Inter.", chr == 9 ~ "Inter.", chr == 10 ~ "Inter.", chr == 11 ~ "Inter.", chr == 12 ~ "Inter.", chr == 13 ~ "Inter.", chr == 14 ~ "Micro.", chr == 15 ~ "Micro.", chr == 16 ~ "Micro.", chr == 17 ~ "Micro.", chr == 18 ~ "Micro.", chr == 19 ~ "Micro.", chr == 20 ~ "Micro.", chr == 21 ~ "Micro.", chr == 22 ~ "Micro.", chr == 23 ~ "Micro.", chr == 24 ~ "Micro.", chr == 25 ~ "Micro.", chr == 26 ~ "Micro.", chr == 27 ~ "Micro.", chr == 28 ~ "Micro.", chr == 29 ~ "Micro.", chr == 30 ~ "Micro.", chr == 31 ~ "Excluded", chr == 32 ~ "Micro.", chr == 33 ~ "Excluded", chr == 34 ~ "Excluded", chr == 35 ~ "Micro.", chr == 36 ~ "Micro.", chr == 37 ~ "Micro.", chr == 38 ~ "Micro.", chr == 39 ~ "Micro.", chr == "W" ~ "Inter.", chr == "Z" ~ "Macro."))

library(ggplot2)
pdf("panelD.pdf", width = 8, height = 5)
ggplot(repeats_scatter, aes(log(size), percent, color=type, label=chr)) +
  geom_point() +
  scale_color_manual(values=c("grey", "#ff006d", "#ffb500", "#007eff")) +
  geom_text_repel(min.segment.length=0.5,box.padding=0.3,size=4, nudge_x = 0, nudge_y = 0.1 ) +
  theme_light() +
  labs (y= "%Repeats",  x="Chr size (log)") 
dev.off()