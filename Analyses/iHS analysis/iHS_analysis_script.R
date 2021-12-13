library(rehh)
library(qvalue)
library(writexl)
library(readxl)

rm(list = ls())


##H. r. savignii population
#read in vcf file
egy_wgs_chr6 <- data2haplohh(hap_file = "egy/SUPER_6_egy_WGS_phased_nomiss.vcf.gz",
                         polarize_vcf = FALSE,
                         vcf_reader = "data.table",
                         min_perc_geno.mrk = 90)

#filter for maf and genotype missingness 
egy_wgs_chr6_filtered <- subset(egy_wgs_chr6, 
                                min_perc_geno.mrk = 90,
                                min_maf = 0.05)

#calculate iES statistics by using scan_hh()
#this function calculates these values for all markers in the haplohh object
egy_wgs_chr6_scan <- scan_hh(egy_wgs_chr6_filtered, polarized = FALSE)

#Use the output of this scan to calculate iHS with the ihh2ihs function
egy_wgs_chr6_ihs <- ihh2ihs(egy_wgs_chr6_scan, freqbin = 1)

#convert to data frame
egy_wgs_chr6_ihs<-as.data.frame(egy_wgs_chr6_ihs)


##FDR correction
#convert logpvalues to pvalues
Pvalues <- data.frame(10^-(abs(egy_wgs_chr6_ihs$ihs.LOGPVALUE)))

#calculate qvalues
qobj <- qvalue(p = Pvalues)

qvalues <- qobj$qvalues

colnames(qvalues) <- "qvalue"

#FDR corrected dataset
egy_wgs_chr6_ihs_FDR <- data.frame(egy_wgs_chr6_ihs$ihs.CHR, egy_wgs_chr6_ihs$ihs.POSITION, egy_wgs_chr6_ihs$ihs.IHS, qvalues$qvalue)

#set 5% as significative threshold
egy_wgs_chr6_ihs_FDR_0.05 <- subset.data.frame(egy_wgs_chr6_ihs_FDR, qvalues.qvalue < 0.05)

#write excel file containing markers with significative qvalues
write_xlsx(egy_wgs_chr6_ihs_FDR_0.05, "egy_wgs_chr6_ihs_FDR_0.05.xlsx")


##plot
#import significant values after FDR correction (present in the file "egy_wgs_chr6_ihs_FDR_0.05_dataframe.xlsx")
egy_FDR_markers<-read_excel ("egy_wgs_chr6_ihs_FDR_0.05_dataframe.xlsx")
egy_FDR_markers<-as.data.frame(egy_FDR_markers)

#import coordinates of the high-LD region containing BDNF gene
bdnf_region<-read_excel("chr6_bdnf_region.xlsx")
bdnf_region<-as.data.frame(bdnf_region)


#manhattanplot highlighting significant values after FDR correction and BDNF region
manhattanplot(egy_wgs_chr6_ihs, mrk = egy_FDR_markers, cex=0.3, mrk.lab.cex =0, mrk.pch=1, mrk.col="red", mrk.cex=0.5, cr=bdnf_region, cr.col="green4", cr.opacity=.5, cex.axis=2, cex.lab=1.5, yaxt="n")
axis(2, seq(-8,8,2), las=2, cex.axis=2)


##H. r. erythrogaster
#read in vcf file
ame_wgs_chr6 <- data2haplohh(hap_file = "ame/SUPER_6_ame_WGS_phased_nomiss.vcf.gz",
                             polarize_vcf = FALSE,
                             vcf_reader = "data.table",
                             min_perc_geno.mrk = 90)


#filter for maf and genotype missingness 
ame_wgs_chr6_filtered <- subset(ame_wgs_chr6,
                                min_perc_geno.mrk = 90,
                                min_maf = 0.05)


#calculate iES statistics by using scan_hh()
#this function calculates these values for all markers in the haplohh object
ame_wgs_chr6_scan <- scan_hh(ame_wgs_chr6_filtered, polarized = FALSE)

#Use the output of this scan to calculate iHS with the ihh2ihs function
ame_wgs_chr6_ihs <- ihh2ihs(ame_wgs_chr6_scan, freqbin = 1)

#convert to data frame
ame_wgs_chr6_ihs<-as.data.frame(ame_wgs_chr6_ihs)


##FDR correction
#convert logpvalues to pvalues
Pvalues <- data.frame(10^-(abs(ame_wgs_chr6_ihs$ihs.LOGPVALUE)))

#calculate qvalues
qobj <- qvalue(p = Pvalues)

qvalues <- qobj$qvalues

colnames(qvalues) <- "qvalue"

#FDR corrected dataset
ame_wgs_chr6_ihs_FDR <- data.frame(ame_wgs_chr6_ihs$ihs.CHR, ame_wgs_chr6_ihs$ihs.POSITION, ame_wgs_chr6_ihs$ihs.IHS, qvalues$qvalue)

#set 5% as significative threshold
ame_wgs_chr6_ihs_FDR_0.05 <- subset.data.frame(ame_wgs_chr6_ihs_FDR, qvalues.qvalue < 0.05)


#write excel file containing markers with significative qvalues
write_xlsx(ame_wgs_chr6_ihs_FDR_0.05, "ame_wgs_chr6_ihs_FDR_0.05.xlsx")

##plot
#import significant values after FDR correction (present in the file "ame_wgs_chr6_ihs_FDR_0.05_dataframe.xlsx")
ame_FDR_markers<-read_excel("ame_wgs_chr6_ihs_FDR_0.05_dataframe.xlsx")
ame_FDR_markers<-as.data.frame(ame_FDR_markers)

#manhattanplot highlighting significant values after FDR correction and BDNF region
manhattanplot(ame_wgs_chr6_ihs, mrk = ame_FDR_markers, cex=0.3, mrk.lab.cex =0, mrk.pch=1, mrk.col="red", mrk.cex=0.5, cr=bdnf_region, cr.col="green4", cr.opacity=.5, cex.axis=2, cex.lab=1.5, yaxt="n")
axis(2, seq(-8,8,2), las=2, cex.axis=2)


##H. r. gutturalis
#read in vcf file
gutturalis_chr6 <- data2haplohh(hap_file = "gutturalis/gutturalis_SUPER_6_phased_nomiss.vcf.gz",
                             polarize_vcf = FALSE,
                             vcf_reader = "data.table", 
                             min_perc_geno.mrk = 90)

#filter for maf and genotype missingness 
gutturalis_chr6_filtered <- subset(gutturalis_chr6,
                                   min_perc_geno.mrk = 90,
                                   min_maf = 0.05
                                   )

#calculate iES statistics by using scan_hh()
#this function calculates these values for all markers in the haplohh object
gutturalis_chr6_scan <- scan_hh(gutturalis_chr6_filtered, polarized = FALSE)

#Use the output of this scan to calculate iHS with the ihh2ihs function
gutturalis_chr6_ihs <- ihh2ihs(gutturalis_chr6_scan, freqbin = 1)

#convert to data frame
gutturalis_chr6_ihs<-as.data.frame(gutturalis_chr6_ihs)

##FDR correction
#convert logpvalues to pvalues
Pvalues <- data.frame(10^-(abs(gutturalis_chr6_ihs$ihs.LOGPVALUE)))

#calculate qvalues
qobj <- qvalue(p = Pvalues)

qvalues <- qobj$qvalues

colnames(qvalues) <- "qvalue"

#FDR corrected dataset
gutturalis_chr6_ihs_FDR <- data.frame(gutturalis_chr6_ihs$ihs.CHR, gutturalis_chr6_ihs$ihs.POSITION, gutturalis_chr6_ihs$ihs.IHS, qvalues$qvalue)

#set 5% as significative threshold
gutturalis_chr6_ihs_FDR_0.05 <- subset.data.frame(gutturalis_chr6_ihs_FDR, qvalues.qvalue < 0.05)

#write excel file containing markers with significative qvalues
write_xlsx(gutturalis_chr6_ihs_FDR_0.05, "gutturalis_wgs_chr6_ihs_FDR_0.05.xlsx")


##plot
#import significant values after FDR correction (present in the file "gutturalis_wgs_chr6_ihs_FDR_0.05_dataframe.xlsx")
gutt_FDR_markers<-read_excel ("gutturalis_wgs_chr6_ihs_FDR_0.05_dataframe.xlsx")
gutt_FDR_markers<-as.data.frame(gutt_FDR_markers)

#manhattanplot highlighting significant values after FDR correction and BDNF region
manhattanplot(gutturalis_chr6_ihs, mrk = gutt_FDR_markers, cex=0.3, mrk.lab.cex =0, mrk.pch=1, mrk.col="red", mrk.cex=0.5, cr=bdnf_region, cr.col="green4", cr.opacity=.5, cex.axis=2, cex.lab=1.5, yaxt="n")
axis(2, seq(-8,8,2), las=2, cex.axis=2)


##H. r. rustica
#read in vcf file
rustica_chr6 <- data2haplohh(hap_file = "rustica_rustica/rustica_SUPER_6_phased_nomiss.vcf.gz",
                                polarize_vcf = FALSE,
                                vcf_reader = "data.table",
                                min_perc_geno.mrk = 90)


#filter for maf and genotype missingness 
rustica_chr6_filtered <- subset(rustica_chr6,
                                min_perc_geno.mrk = 90,
                                min_maf = 0.05)


#calculate iES statistics by using scan_hh()
#this function calculates these values for all markers in the haplohh object
rustica_chr6_scan <- scan_hh(rustica_chr6_filtered, polarized = FALSE)

#Use the output of this scan to calculate iHS with the ihh2ihs function
rustica_chr6_ihs <- ihh2ihs(rustica_chr6_scan, freqbin = 1)

#convert to data frame
rustica_chr6_ihs<-as.data.frame(rustica_chr6_ihs)


##FDR correction
#convert logpvalues to pvalues
Pvalues <- data.frame(10^-(abs(rustica_chr6_ihs$ihs.LOGPVALUE)))

#calculate qvalues
qobj <- qvalue(p = Pvalues)

qvalues <- qobj$qvalues

colnames(qvalues) <- "qvalue"

#FDR corrected dataset
rustica_chr6_ihs_FDR <- data.frame(rustica_chr6_ihs$ihs.CHR, rustica_chr6_ihs$ihs.POSITION, rustica_chr6_ihs$ihs.IHS, qvalues$qvalue)

#set 5% as significative threshold
rustica_chr6_ihs_FDR_0.05 <- subset.data.frame(rustica_chr6_ihs_FDR, qvalues.qvalue < 0.05)

#write excel file containing markers with significative qvalues
write_xlsx(rustica_chr6_ihs_FDR_0.05, "rustica_chr6_ihs_FDR_0.05.xlsx")


##plot
#import significant values after FDR correction (present in the file "rustica_chr6_ihs_FDR_0.05_dataframe.xlsx")
rustica_FDR_markers<-read_excel("rustica_chr6_ihs_FDR_0.05_dataframe.xlsx")
rustica_FDR_markers<-as.data.frame(rustica_FDR_markers)

#manhattanplot highlighting significant values after FDR correction and BDNF region
manhattanplot(rustica_chr6_ihs, mrk = rustica_FDR_markers, cex=0.3, mrk.lab.cex=0, mrk.pch=1, mrk.col="red", mrk.cex=0.5, cr=bdnf_region, cr.col="green4", cr.opacity=.5, cex.axis=2, cex.lab=1.5, yaxt="n")
axis(2, seq(-8,8,2), las=2, cex.axis=2)