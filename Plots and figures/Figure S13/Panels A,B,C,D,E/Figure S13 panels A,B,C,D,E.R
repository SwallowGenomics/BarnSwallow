library(readxl)
library(ggplot2)
library(ggpubr)


#import dataset and convert to data.frame
df1 <- data.frame(read_excel("genomecovA1.xlsx"))
df2 <- data.frame(read_excel("genomecovA2.xlsx"))
df3 <- data.frame(read_excel("Genomecov2.xlsx"))
df4 <- data.frame(read_excel("Genomecov4.xlsx"))
df5 <- data.frame(read_excel("Genomecov3.xlsx"))


#Plot Genomecov Sample A1
a <- ggplot(df1, aes(x=Coverage, y=Frazione.di.basi.coverage)) +
    geom_bar(stat = "identity", width= 0.5, color="cornflowerblue", fill="cornflowerblue") + 
    theme_bw() +
    theme(axis.text = element_text(size=25))+
  theme(axis.title = element_text(size=20))+						#Font size 
  theme(plot.title = element_text(size=30))+						#Font size
    ggtitle("Sample A1", ) + labs(y="Fraction of bases")+			#Title and legend
    xlim(0, 65) +
    annotate("text", x = 10, y = 0.037, label = "Expected\ncoverage\n33x", size = 7) 	#Adding expected coverage 
  
  
#Plot Genomecov Sample A2
b <- ggplot(df2, aes(x=Coverage, y=Frazioni.di.basi.coverage)) +
    geom_bar(stat = "identity", width= 0.5, color="cornflowerblue", fill="cornflowerblue") + 
    theme_bw() + 
  theme(axis.text = element_text(size=25))+							#Font size 
  theme(axis.title = element_text(size=20))+                        #Font size
  theme(plot.title = element_text(size=30))+                        #Title and legend
    ggtitle("Sample A2") + labs(y="Fraction of bases")+
    xlim(0, 35)+
  annotate("text", x = 4.5, y = 0.060, label = "Expected\ncoverage\n15x", size = 7) 		#Adding expected coverage 

 
#Plot Genomecov Sample 2
c <- ggplot(df3, aes(x=Coverage, y=Frazione.di.basi.coverage)) +
    geom_bar(stat = "identity", width= 0.5, color="cornflowerblue", fill="cornflowerblue") + 
    theme_bw() + 
  theme(axis.text = element_text(size=25))+							#Font size 
  theme(axis.title = element_text(size=20))+                        #Font size
  theme(plot.title = element_text(size=30))+                        #Title and legend
    ggtitle("Sample 2") + labs(y="Fraction of bases")+
    xlim(0, 50)+
  annotate("text", x = 7, y = 0.045, label = "Expected\ncoverage\n25x", size = 7) 		#Adding expected coverage 


#Plot Genomecov Sample 3
d <- ggplot(df4, aes(x=Coverage, y=Frazione.di.basi.coverage)) +
    geom_bar(stat = "identity", width= 0.5, color="cornflowerblue", fill="cornflowerblue") + 
    theme_bw() + 
  theme(axis.text = element_text(size=25))+							#Font size 
  theme(axis.title = element_text(size=20))+                        #Font size
  theme(plot.title = element_text(size=30))+                        #Title and legend
    ggtitle("Sample 4") + labs(y="Fraction of bases")+
    xlim(0, 45)+
  annotate("text", x = 7, y = 0.050, label = "Expected\ncoverage\n19x", size = 7) 	#Adding expected coverage 


#Plot Genomecov Sample 4
e <- ggplot(df5, aes(x=Coverage, y=Frazione.di.basi)) +
    geom_bar(stat = "identity", width= 0.5, color="cornflowerblue", fill="cornflowerblue") + 
    theme_bw() + 
  theme(axis.text = element_text(size=25))+							#Font size 
  theme(axis.title = element_text(size=20))+                        #Font size
  theme(plot.title = element_text(size=30))+                        #Title and legend
    ggtitle("Sample 3") + labs(y="Fraction of bases")+
    xlim(0, 40)+
  annotate("text", x = 5, y = 0.050, label = "Expected\ncoverage\n20x", size = 7) 	#Adding expected coverage 


#Put the five plots togheter
ggarrange(a, b, c, d, e, labels = c("a", "b", "c", "d", "e"),font.label = list(size = 30), ncol = 2 , nrow = 3)


    
