rm(list = ls())

library(ggplot2)
library(readxl)
library(ggtext)


## LD decay plot

#import data
LD_whole<-read_excel("LD_decay_all_data.xlsx")

#convert to data frame
LD_whole<-as.data.frame(LD_whole)

#geom_line plot
a <- ggplot()+
  geom_line(data=LD_whole,aes(y=erythrogaster,x=window,colour="erythrogaster"),size=2 )+
  geom_line(data=LD_whole,aes(y=savignii,x=window,colour="savignii"),size=2 ) +
  geom_line(data=LD_whole,aes(y=gutturalis,x=window,colour="gutturalis"),size=2 ) +
  geom_line(data=LD_whole,aes(y=rustica,x=window,colour="rustica"),size=2 ) +
  geom_line(data=LD_whole,aes(y=rustica_gutturalis,x=window,colour="rustica_gutturalis"),size=2 ) +
  geom_line(data=LD_whole,aes(y=rustica_tytleri,x=window,colour="rustica_tytleri"),size=2 ) +
  geom_line(data=LD_whole,aes(y=transitiva,x=window,colour="transitiva"),size=2 ) +
  geom_line(data=LD_whole,aes(y=tytleri,x=window,colour="tytleri"),size=2 ) +
  geom_line(data=LD_whole,aes(y=gutturalis_tytleri,x=window,colour="gutturalis_tytleri"),size=2 ) +
  scale_color_manual(values = c("erythrogaster" = "#e6194B", "savignii" = "#000000", "gutturalis" = "#800000", "rustica" = "#808000", "rustica_gutturalis" = "#FF00AA", "rustica_tytleri" = "#911eb4", "transitiva" = "#4363d8", "tytleri"= "#42d4f4", "gutturalis_tytleri"= "#469990"))

#customize axis titles; graphical adjustments
a <- a + labs(y = "LD (r^2)")
a <- a + labs(x = "distance (kbp)")
a <- a + theme_light()

#mark more y values
a <- a + scale_y_continuous(breaks=seq(0,0.2,0.01))

#set black and white theme
a <- a + theme_bw()

#remove grid lines
a <- a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), panel.border=element_blank(), axis.line = element_line(colour = "black"))

#edit legend labels content; specify italic font for subspecies name using ggtext extension
a <- a + scale_colour_manual(name = "Sample size", breaks = c("savignii", "erythrogaster", "transitiva", "tytleri", "rustica_tytleri", "rustica_gutturalis", "gutturalis", "rustica", "gutturalis_tytleri"), values=c("#000000", "#e6194B", "#4363d8", "#42d4f4", "#911eb4", "#FF00AA", "#800000", "#808000", "#469990"), labels = c("*H. r. savignii* (N=8)", "*H. r. erythrogaster* (N=8)", "*H. r. transitiva* (N=8)", "*H. r. tytleri* (N=10)", "*H. r. rusticaxtytleri* (N=16)", "*H. r. rusticaxgutturalis* (N=21)", "*H. r. gutturalis* (N=34)", "*H. r. rustica* (N=25)", "*H. r. gutturalisxtyleri* (N=29)"), guide=guide_legend(ncol=3))

#adjust legend text size
a <- a + theme(legend.text = element_markdown(size=18))

#change legend position and orientation
a<- a + theme(legend.position="top", legend.direction="horizontal")

#hide title legend
a <- a + theme(legend.title = element_blank())

#adjust axis title size
a <- a + theme(
  axis.title = element_text(size=17)
)


#increase axis values size
a <- a + theme(axis.text = element_text(size = 15))


a

