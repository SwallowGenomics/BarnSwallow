library(ggplot2)
library(readxl)
library(gridExtra)
library(grid)
library(scales)


rm(list = ls())


#import datasets
SV_down<-read_excel("SV_down.xlsx")
SV_full<-read_excel("SV_full.xlsx")
SNVsdown<-read_excel("SNVs_down.xlsx")
SNVsfull<-read_excel("SNVs_full.xlsx")
Indelsdown<-read_excel("INdels_down.xlsx")
Indelsfull<-read_excel("INdels_full.xlsx")


#convert to dataframe
SV_down<-as.data.frame(SV_down)
SV_full<-as.data.frame(SV_full)
SNVsdown<-as.data.frame(SNVsdown)
SNVsfull<-as.data.frame(SNVsfull)
Indelsdown<-as.data.frame(Indelsdown)
Indelsfull<-as.data.frame(Indelsfull)


##plot SVs
a <- ggplot()+
 geom_line(data=SV_full, aes(x=sample, y=variant_number, group=variant_type, color=variant_type))+
 geom_line(data=SV_down, linetype="dashed", aes(x=sample, y=variant_number, group=variant_type, color=variant_type))+
  scale_color_manual(values = c("red", "blue"))
  

#adjust background
a <- a + theme_bw()

#adjust y limits
a <- a + ylim(25000,125000)


#remove axis titles
a <- a + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#hide legend
a <- a + theme(legend.position = "none")

#increase axis values size
a <- a + theme(axis.text = element_text(size = 12))

#put title
a <- a + ggtitle("SVs")

#center title
a <- a + theme(plot.title = element_text(hjust = 0.5)) 

#change scale (to k-values)
a <- a + scale_y_continuous(labels = label_number(suffix = " k", scale = 1e-3), breaks = seq(20000,130000,10000)) 


##plot indels
b <- ggplot()+
  geom_line(data=Indelsfull, aes(x=sample, y=variant_number, group=variant_type, color=variant_type))+
  geom_line(data=Indelsdown, linetype="dashed", aes(x=sample, y=variant_number, group=variant_type, color=variant_type))+
  scale_color_manual(values = c("red", "blue"))


#adjust background
b <- b + theme_bw()

#adjust y limits
b <- b + ylim(500000,1750000)


#remove axis titles
b <- b + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#hide legend
b <- b + theme(legend.position = "none")

#increase axis values size
b <- b + theme(axis.text = element_text(size = 12))

#put title
b <- b + ggtitle("indels")

#center title
b <- b + theme(plot.title = element_text(hjust = 0.5)) 

#change scale (to k-values)
b <- b + scale_y_continuous(labels = label_number(suffix = " k", scale = 1e-3), breaks = seq(500000,2000000,200000)) 


##plot SNVs
c <- ggplot()+
  geom_line(data=SNVsfull, aes(x=sample, y=variant_number, group=variant_type, color=variant_type))+
  geom_line(data=SNVsdown, linetype="dashed", aes(x=sample, y=variant_number, group=variant_type, color=variant_type))+
  scale_color_manual(values = c("red", "blue"))


#adjust background
c <- c + theme_bw()

#adjust y limits
c <- c + ylim(2500000,12500000)

#remove axis titles
c <- c + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#hide legend
c <- c + theme(legend.position = "none")

#increase axis values size
c <- c + theme(axis.text = element_text(size = 12))

#put title
c <- c + ggtitle("SNVs")

#center title
c <- c + theme(plot.title = element_text(hjust = 0.5)) 

#change scale (to k-values)
c <- c + scale_y_continuous(labels = label_number(suffix = " k", scale = 1e-3), breaks = seq(2500000,15000000, 1000000)) 


#Put the three plots together
grid.arrange (c, b, a, ncol=2)


