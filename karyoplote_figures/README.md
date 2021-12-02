# karyoploteR figures

The R files present in the folders machrocromosomes and microchromosomes contain the commands used to generate Fig. 3 and Extended Data Fig. 7.
For Suppl. Fig. 9, the script is the same as the one used for machrocromosomes, except that a downsampled (5x) set of SNPs has been used to plot HiFi data.
The "chrom_list" file contains chromosome sizes and is requested by the karyoploteR package in order to plot ideograms using custom genomes.
Files required to plot the additional tracks of the figure (assembly gaps, repeats, GC content, Pacbio reads coverage) have been provided for chromosome 1 as example files.
The "chr_1_snp_pos" files contain SNP positions for all sequencing technologies analyzed, extracted directly from the relative vcf files with the bcftools query command.