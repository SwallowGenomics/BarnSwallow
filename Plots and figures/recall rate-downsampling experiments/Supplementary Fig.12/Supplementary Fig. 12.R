library(readxl)
library(ggplot2)
library(tidyverse)


#import dataset
downsampling <- read_excel("Downsampling.xlsx")


#convert to dataframe
df <- data.frame(downsampling)


#plot downsampling chart
ggplot(df, aes(x = Coverage,)) +			
  geom_line(aes(y=Recall, colour = "Recall"))+				#plot Recall line
  geom_line(aes(y=Precision, colour ="Precision"))+			#plot Precision line
  geom_point(data = df, aes(Coverage, Recall)) +			
  geom_point(data = df, aes(Coverage, Precision)) +
  scale_x_reverse()+										
  geom_errorbar(data=df, aes(Coverage, Recall, ymin = Recall-SD_recall, ymax = Recall+SD_recall), width=0.8)+					#add SD_information for Recall
  geom_errorbar(data=df, aes(Coverage, Precision, ymin=Precision-SD_precision, ymax=Precision+SD_precision), width=0.8) +		#add SD information for Precision	
  theme_bw()+
  theme(axis.text = element_text(size=10))+							#font size
  theme(axis.title = element_text(size=10))+						#font size
  scale_color_manual("Rate",
                     breaks = c("Recall", "Precision"),
                     values = c("dodgerblue", "green"))+			#colour setting and legend setting 
  xlab("Coverage")+
  ylab("Rate")+
  ylim(0,1)

