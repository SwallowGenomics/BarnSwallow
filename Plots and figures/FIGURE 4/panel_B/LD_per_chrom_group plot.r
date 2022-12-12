rm(list = ls())

library(ggplot2)
library(readxl)


##LD_per_chrom_group plot

#import data
LD_per_group <- read_excel("LD_values_boxplot.xlsx")

#convert to data frame
LD_per_group <- as.data.frame(LD_per_group)

#convert distance bins to factors
LD_per_group$distance <- as.factor(LD_per_group$distance)

#set the correct order for plotting distance bins
level_order <- c('0-1', '1-2', '2-3', '3-4', '4-5', '5-7', '7-9', '9-11', '11-13', '13-15', '15-20', '20-25', '25-30', '30-40', '40-50', '50-60')

#set the correct order for chromosome groups: turn the 'type' column into a character vector
LD_per_group$chr_type <- as.character(LD_per_group$chr_type)

#then turn it back into a factor with the levels in the correct order
LD_per_group$chr_type <- factor(LD_per_group$chr_type, levels=unique(LD_per_group$chr_type))

#assign colors to chromosome groups
group.colors <- c(macro = "#ffb500", intermediate = "#ff006d", micro ="#007eff")

#boxplot
b <- ggplot(LD_per_group, aes(x=factor(distance,level=level_order), y=LD, fill=chr_type)) + geom_boxplot() + scale_fill_manual(values=group.colors)

#mark more y values
b <- b + scale_y_continuous(breaks=seq(0,0.2,0.01))

#customize axis titles; graphical adjustments
b <- b + labs(x= "distance (kbp)", y = "LD (r^2)")
b <- b + theme_light()

#set black and white theme
b <- b + theme_bw()

#remove grid lines
b <- b + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), panel.border=element_blank(), axis.line = element_line(colour = "black"))


#change legend position and orientation
b <- b + theme(legend.position="top", legend.direction="horizontal")

#adjust legend size
b <- b + theme(
  legend.text = element_text(size=18)
)

b <- b + theme(
  legend.title = element_text(size=18)
)


#adjust axis title size
b <- b + theme(
  axis.title = element_text(size=17)
)

#increase axis values size
b <- b + theme(axis.text = element_text(size = 15))


b