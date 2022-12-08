setwd("E:/path/to/files/")

table <- read.table("sizes_and_estimations.txt", header=TRUE)

library(ggplot2)
library(ggrepel)

pdf("FIGURE2_panelB.pdf", width = 8, height = 5)

ggplot(stime, aes(x=size, y=mean, label=chr, color=type)) +
  geom_smooth(method=lm, color="dark grey", size=0.6) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=1500000, size=0.4) +
  geom_point(aes(stroke=0.4), size=2.5, shape=21, fill="black") +
  scale_color_manual(values=c("#ff006d", "#ffb500", "#007eff")) +
  geom_text_repel(color="black",min.segment.length=5,box.padding=0.2,size=3.5, nudge_x = 0.1, nudge_y = 0.2) +
  theme_light() +
  labs(y= "Chr. karyotype size (bp)",  x="Chr. assembly size (bp)") 

dev.off()

