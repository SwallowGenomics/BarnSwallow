setwd("E:/path/to/folder")

library(ggplot2)
library(reshape2)
library(tidyr)

#upper part
barplot <- read.table("barplot.txt", header=TRUE)

pdf("barplot.pdf", width = 31, height = 6)

ggplot(data=barplot, aes(x=variable, y=genes)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=genes), vjust=-1, color="black", size=10)+
  scale_x_discrete(limits=barplot$variable, expand=c(0,0)) +
  scale_y_continuous(trans='log10',name = "Genes")+
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

dev.off()

#lower
table1 <- read.table("matrix.txt", header=TRUE)

table_melt1 <- melt(table1)

pdf("matrix_black.pdf", width = 31, height = 5)

ggplot() +
  geom_tile(data=table_melt1, aes(x=variable, y=indiv,fill=as.factor(value)), colour="white", show.legend = FALSE) +
  coord_equal(expand = 0) +
  scale_fill_manual(values=c("white", "black", "#dad7d7")) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") +
  scale_y_discrete(limits = unique(rev(table$indiv))) +
  theme( axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(color = "black", size = 40, angle = 0, hjust = 1, vjust = 0.4, face = "plain"))

dev.off()