#FDR correction on phyloP basewise

setwd("/path/to/file")
bed_basewise <- read.table("PhyloP_basewise.bed") #read the file
colnames(bed_basewise) <- c("chr", "start", "end", "id", "logP")
nrow(bed_basewise) #988249638
range(bed_basewise$V5) #-7.030  0.597

Pvalues <- data.frame(10^-(abs(bed_basewise$logP))) #logPvalues to Pvalues
library(qvalue) #lead the library
qobj <- qvalue(p = Pvalues) #calculate qvalues
bed_FDR_basewise <- data.frame(bed_basewise, qobj$qvalues) #add qvalues to the data frame 
colnames(bed_FDR_basewise) <- c("chr", "start", "end", "id", "logP", "qvalue")

bed_0.05_basewise <- subset.data.frame(bed_FDR_basewise, qvalue < 0.05) 
nrow(bed_0.05_basewise) #0

#FDR correction on phyloP 10bp windows

bed <- read.table("phyloP_10bp.bed")
colnames(bed) <- c("chr", "start", "end", "logP")
nrow(bed) #98829932

bed$start <- (bed$start-1) #change start from 1-based to 0-based
bed$length <- (bed$start-bed$end) #calculate length
bed1 <- subset.data.frame(bed, length>=10) #remove windows smaller that 10bp
nrow(bedd1) #98820052 (removed 9880 windows)
range(bedd$logP) #20.000 3.508 

Pvalues <- data.frame(10^-(abs(bed1$logP))) #logPvalues in Pvalues
library(qvalue) #lead library
qobj <- qvalue(p = Pvalues) #calculate qvalues
qvalues <- qobj$qvalues #extract qvalue column
colnames(qvalues) <- "qvalue"

bed_FDR <- data.frame(bed1$chr,bed1$start, bed1$end, bed1$logP, qvalues$qvalue)
colnames(bed_FDR) <- c("chr", "start", "end", "logP", "qvalue")
bed_0.05 <- subset.data.frame(bedd_FDR, qvalues.qvalue < 0.05) 

nrow(bed_0.05) #3967050 number of FDR corrected windows

acc <- subset.data.frame(bedd_0.05 , logP<0) #extract accelerated windows
cons <- subset.data.frame(bedd_0.05, logP>0) #extract conserved windows
nrow(acc) #1036529
nrow(cons) #2930521
range(acc$logP) #--20.000  -2.698 
range(cons$logP) #2.698 3.508

#Bonferroni
soglia <- -log((0.05/98820052),10) #9.295875 

Bonferroni_acc <- subset.data.frame(acc, logP<(-9.295836))
Bonferroni_cons <- subset.data.frame(cons, logP>9.295836)
nrow(Bonferroni_acc) #6373 
nrow(Bonferroni_cons) #0

write.table(bed_0.05, file="PhyloP_CHR_10bp_FDR_4fold.bed", sep = "\t", col.names=FALSE, row.names=FALSE)
write.table(acc, file="PhyloP_CHR_10bp_FDR_4fold_ACC.bed", sep = "\t", col.names=FALSE, row.names=FALSE)
write.table(cons, file="PhyloP_CHR_10bp_FDR_4fold_CONS.bed", sep = "\t", col.names=FALSE, row.names=FALSE)



