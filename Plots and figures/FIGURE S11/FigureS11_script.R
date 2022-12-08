library(ggplot2)
library(readxl)
library(gridExtra)
library(grid)
library(scales)
library(dplyr)
library(ggpubr)

rm(list = ls())



png("macro.png", width = 3500, height = 2500, res=300)


##macrochromosomes
#import datasets
#egy
macro_egy<-read_excel("Egy_Smith/macro_Egy.xlsx")

macro_egy <- macro_egy %>% 
  mutate(midpoint = midpoint / 1000)


#ame
macro_ame<-read_excel("Ame_Smith/macro_Ame.xlsx")

macro_ame <- macro_ame %>% 
  mutate(midpoint = midpoint / 1000)


#gutt
macro_gutt<-read_excel("macro_gutturalis.xlsx")

macro_gutt <- macro_gutt %>% 
  mutate(midpoint = midpoint / 1000)


#rustica
macro_rustica<-read_excel("macro_rustica.xlsx")

macro_rustica <- macro_rustica %>% 
  mutate(midpoint = midpoint / 1000)


##plot 
#egy
egy <- ggplot (macro_egy, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
egy <- egy + theme_bw()

#adjust axis titles
egy <- egy + labs(x = "distance (kbp)")
egy <- egy + labs(y = "LD")


#mark more x vaues
egy <- egy + scale_x_continuous(breaks=seq(0,80000,10000))


#ame
ame <- ggplot (macro_ame, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
ame <- ame + theme_bw()

#adjust axis titles
ame <- ame + labs(x = "distance (kbp)")
ame <- ame + labs(y = "LD")



#mark more x vaues
ame <- ame + scale_x_continuous(breaks=seq(0,80000,10000))


#gutt
gutt <- ggplot (macro_gutt, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
gutt <- gutt + theme_bw()

#adjust axis titles
gutt <- gutt + labs(x = "distance (kbp)")
gutt <- gutt + labs(y = "LD")


#mark more x vaues
gutt <- gutt + scale_x_continuous(breaks=seq(0,80000,10000))


#rustica
rustica <- ggplot (macro_rustica, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
rustica <- rustica + theme_bw()

#adjust axis titles
rustica <- rustica + labs(x = "distance (kbp)")
rustica <- rustica + labs(y = "LD")


#mark more x vaues
rustica <- rustica + scale_x_continuous(breaks=seq(0,80000,10000))


#Plot everything together


figure1 <- ggarrange(ame, egy, gutt, rustica,
          labels = c("a", "b", "c","d"),
          ncol=2, nrow=2)


annotate_figure(figure1, 
                top = text_grob("Macrochromosomes", color = "black", face = "bold", size = 14))

                                

dev.off()


png("int.png", width = 3500, height = 2500, res=300)


##intermediate
#import datasets
#egy
int_egy<-read_excel("Egy_Smith/int_Egy.xlsx")

int_egy <- int_egy %>% 
  mutate(midpoint = midpoint / 1000)


#ame
int_ame<-read_excel("Ame_Smith/int_Ame.xlsx")


int_ame <- int_ame %>% 
  mutate(midpoint = midpoint / 1000)


#gutt
int_gutt<-read_excel("int_gutturalis.xlsx")

int_gutt <- int_gutt %>% 
  mutate(midpoint = midpoint / 1000)


#rustica
int_rustica<-read_excel("int_rustica.xlsx")

int_rustica <- int_rustica %>% 
  mutate(midpoint = midpoint / 1000)


##plot 
#egy
egy <- ggplot (int_egy, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
egy <- egy + theme_bw()

#adjust axis titles
egy <- egy + labs(x = "distance (kbp)")
egy <- egy + labs(y = "LD")


#ame
ame <- ggplot (int_ame, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
ame <- ame + theme_bw()

#adjust axis titles
ame <- ame + labs(x = "distance (kbp)")
ame <- ame + labs(y = "LD")


#gutt
gutt <- ggplot (int_gutt, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
gutt <- gutt + theme_bw()

#adjust axis titles
gutt <- gutt + labs(x = "distance (kbp)")
gutt <- gutt + labs(y = "LD")


#rustica
rustica <- ggplot (int_rustica, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
rustica <- rustica + theme_bw()

#adjust axis titles
rustica <- rustica + labs(x = "distance (kbp)")
rustica <- rustica + labs(y = "LD")




#Plot everything together
figure2 <- ggarrange(ame, egy, gutt, rustica,
                     labels = c("e", "f", "g","h"),
                     ncol=2, nrow=2)


annotate_figure(figure2, 
                top = text_grob("Intermediate chromosomes", color = "black", face = "bold", size = 14))



dev.off()



png("micro.png", width = 3500, height = 2500, res=300)


##microchromosomes
#import datasets
#egy
micro_egy<-read_excel("Egy_Smith/micro_Egy.xlsx")

micro_egy <- micro_egy %>% 
  mutate(midpoint = midpoint / 1000)


#ame
micro_ame<-read_excel("Ame_Smith/micro_Ame.xlsx")


micro_ame <- micro_ame %>% 
  mutate(midpoint = midpoint / 1000)


#gutt
micro_gutt<-read_excel("micro_gutturalis.xlsx")


micro_gutt <- micro_gutt %>% 
  mutate(midpoint = midpoint / 1000)


#rustica
micro_rustica <- read_excel("micro_rustica.xlsx")

micro_rustica <- micro_rustica %>% 
  mutate(midpoint = midpoint / 1000)


##plot 
#egy
egy <- ggplot (micro_egy, aes(x=midpoint, y=R2)) + geom_point(size=.8)


#adjust background
egy <- egy + theme_bw()

#adjust axis titles
egy <- egy + labs(x = "distance (kbp)")
egy <- egy + labs(y = "LD")


#ame
ame <- ggplot (micro_ame, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
ame <- ame + theme_bw()

#adjust axis titles
ame <- ame + labs(x = "distance (kbp)")
ame <- ame + labs(y = "LD")


#gutt
gutt <- ggplot (micro_gutt, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
gutt <- gutt + theme_bw()

#adjust axis titles
gutt <- gutt + labs(x = "distance (kbp)")
gutt <- gutt + labs(y = "LD")


#rustica
rustica <- ggplot (micro_rustica, aes(x=midpoint, y=R2)) + geom_point(size=.8)

#adjust background
rustica <- rustica + theme_bw()

#adjust axis titles
rustica <- rustica + labs(x = "distance (kbp)")
rustica <- rustica + labs(y = "LD")




#Plot everything together
figure3 <- ggarrange(ame, egy, gutt, rustica,
                     labels = c("i", "l", "m","n"),
                     ncol=2, nrow=2)


annotate_figure(figure3, 
                top = text_grob("Microchromosomes", color = "black", face = "bold", size = 14))



dev.off()


