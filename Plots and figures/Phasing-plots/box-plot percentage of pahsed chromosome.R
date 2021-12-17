library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(viridis)


#import dataset
phasing <- read_excel("phasing con chr in numeri e dimensione.xlsx")

#convert to data.frame
df <- as.data.frame(phasing3)


#Order chromosomes 
level_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11","12","13","14","15","16", "17", "18", "19", "20", "21", "22","23","24","25","26","27","28","29","30", "31","32","33","34","35","36","37", "38","Z","W")





ggplot(df, aes(x=factor(chromosome, levels = level_order), y=percentage, fill=Type)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#ff006d", "#ffb500", "#007eff", "purple")) + theme_light() + 
  labs(y='Percentage of phased chromosome', x='Chromosome') +
  scale_y_continuous(labels = scales::percent)+				#put y axis in percentage
  theme(text = element_text(size=40))
  axis.text = element_text(size = 35)





