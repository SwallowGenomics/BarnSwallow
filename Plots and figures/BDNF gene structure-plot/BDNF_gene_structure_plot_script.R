rm(list = ls())

library(ggplot2)
library(readxl)


#import datasets for the different isoforms
big_exons_X1<-read_excel("X1_big_exons.xlsx")
small_exons_X1<-read_excel("X1_small_exons.xlsx")
whole_X1<-read_excel("X1_whole.xlsx")

big_exons_X2<-read_excel("X2_big_exons.xlsx")
small_exons_X2<-read_excel("X2_small_exons.xlsx")
whole_X2<-read_excel("X2_whole.xlsx")

big_exons_X3<-read_excel("X3_big_exons.xlsx")
small_exons_X3<-read_excel("X3_small_exons.xlsx")
whole_X3<-read_excel("X3_whole.xlsx")

big_exons_X4<-read_excel("X4_big_exons.xlsx")
small_exons_X4<-read_excel("X4_small_exons.xlsx")
whole_X4<-read_excel("X4_whole.xlsx")

#convert all data to dataframe
big_exons_X1<-as.data.frame(big_exons_X1)
small_exons_X1<-as.data.frame(small_exons_X1)
whole_X1<-as.data.frame(whole_X1)

big_exons_X2<-as.data.frame(big_exons_X2)
small_exons_X2<-as.data.frame(small_exons_X2)
whole_X2<-as.data.frame(whole_X2)

big_exons_X3<-as.data.frame(big_exons_X3)
small_exons_X3<-as.data.frame(small_exons_X3)
whole_X3<-as.data.frame(whole_X3)

big_exons_X4<-as.data.frame(big_exons_X4)
small_exons_X4<-as.data.frame(small_exons_X4)
whole_X4<-as.data.frame(whole_X4)


#plot all isoforms together
p1 <- ggplot(big_exons_X1) + 
  geom_segment(data=whole_X1, aes(x=start,xend=end,y=1,yend=1)) +
  geom_rect(aes(xmin=start, xmax=end,ymin=0.9,ymax=1.1),
            colour="#282a73", fill="#282a73") +
  geom_rect(data=small_exons_X1, aes(xmin=start, xmax=end,ymin=0.95,ymax=1.05),
            colour="#282a73", fill="#282a73") +
  
  geom_segment(data=whole_X3, aes(x=start,xend=end,y=0.5,yend=0.5)) +
  geom_rect(data=big_exons_X3, aes(xmin=start, xmax=end,ymin=0.4,ymax=0.6),
            colour="#282a73", fill="#282a73") +
  geom_rect(data=small_exons_X3, aes(xmin=start, xmax=end,ymin=0.45,ymax=0.55),
            colour="#282a73", fill="#282a73") +

  geom_segment(data=whole_X2, aes(x=start,xend=end,y=0,yend=0)) +
  geom_rect(data=big_exons_X2, aes(xmin=start, xmax=end,ymin=-0.1,ymax=0.1),
            colour="#282a73", fill="#282a73") +
  geom_rect(data=small_exons_X2, aes(xmin=start, xmax=end,ymin=-0.05,ymax=0.05),
            colour="#282a73", fill="#282a73") +
  

  geom_segment(data=whole_X4, aes(x=start,xend=end,y=-0.5,yend=-0.5)) +
  geom_rect(data=big_exons_X4, aes(xmin=start, xmax=end,ymin=-0.6,ymax=-0.4),
            colour="#282a73", fill="#282a73") +
  geom_rect(data=small_exons_X4, aes(xmin=start, xmax=end,ymin=-0.55,ymax=-0.45),
            colour="#282a73", fill="#282a73") +
  
  ylim(c(-1.5,1.5)) 
  xlim(c(-6,32000))

#add isoform labels:
p1<-p1 + geom_text(
  label="X1", 
  x = 40905, y = 1.15, 
  check_overlap = T
)

p1<-p1 + geom_text(
  label="X3", 
  x = 38582, y = 0.65, 
  check_overlap = T
)

p1<-p1 + geom_text(
  label="X2", 
  x = 20723, y = 0.15, 
  check_overlap = T
)

p1<-p1 + geom_text(
  label="X4", 
  x = 17726, y = -0.35, 
  check_overlap = T
)


#add arrows (indicating the direction of transcription)
p1<-p1 + geom_segment(aes(x = 41304, y = 0.8, xend = 39950, yend = 0.8),
                      arrow = arrow(length = unit(0.2, "cm")))

p1<-p1 + geom_segment(aes(x = 39014, y = 0.3, xend = 37660, yend = 0.3),
                      arrow = arrow(length = unit(0.2, "cm")))

p1<-p1 + geom_segment(aes(x = 21178, y = -0.2, xend = 19824, yend = -0.2),
                      arrow = arrow(length = unit(0.2, "cm")))

p1<-p1 + geom_segment(aes(x = 18101, y = -0.7, xend = 16747, yend = -0.7),
                      arrow = arrow(length = unit(0.2, "cm")))


#add ATG lines
p1<-p1 + geom_segment(aes(x = 39394, y = 1.1, xend = 39494, yend = 1.2), col="#FFAA22")

p1<-p1 + geom_segment(aes(x = 1545, y = 0.6, xend = 1645, yend = 0.7), col="#FFAA22")

p1<-p1 + geom_segment(aes(x = 1545, y = 0.1, xend = 1645, yend = 0.2), col="#FFAA22")

p1<-p1 + geom_segment(aes(x = 1545, y = -0.4, xend = 1645, yend = -0.3), col="#FFAA22")

#add ATG labels
p1<-p1 + geom_text(
  label="ATG", 
  x = 39494, y = 1.25, color = "aquamarine4", 
  check_overlap = T
)

p1<-p1 + geom_text(
  label="ATG", 
  x = 1645, y = 0.75, color = "aquamarine4", 
  check_overlap = T
)

p1<-p1 + geom_text(
  label="ATG", 
  x = 1645, y = 0.25, color = "aquamarine4", 
  check_overlap = T
)

p1<-p1 + geom_text(
  label="ATG", 
  x = 1645, y = -0.25, color = "aquamarine4", 
  check_overlap = T
)



#empty all canvas from the plot and obtain an empty plot
p1 <- p1 + theme_void()


#mark the region containing CpG sites that will be magnified in the other panels
p1 <- p1+ geom_rect(aes(xmin=21349, xmax=21595,ymin=-0.85,ymax=1.1), colour="#282a73", fill=NA)


#import CpG islands coordinates
CpG<-read_excel("CpG_islands.xlsx")
CpG<-as.data.frame(CpG)


#plot CpG islands
p1<-p1 + geom_segment(data=CpG, aes(x=start,xend=end,y=1.5,yend=1.5), col="darkviolet", size=1.5)


#add CpG label
p1<-p1 + geom_text( label="CpG", x = -2, y = 1.5, color = "darkviolet",check_overlap = T)


p1

