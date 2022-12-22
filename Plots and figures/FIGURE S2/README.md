## Figure S2 - correlations

This folder contains the R script to generate all panels of Supplementary figure 4. Some modifications in the final figure were done manually. </br>

The repeats input was not uploaded gived its big size.</br>

GC content was calculated on the entire chromosomes like this:</br>
`bedtools nuc -fi bHirRus1_primary.fasta -bed chr_coords.bed > GC_CONTENT_1.bed`</br>
`awk '{print $1, $2, $3, $5*100}' GC_CONTENT_1.bed > GC_CONTENT.bed`</br>
Add back tabs.

CpG islands were downloaded from the [UCSC browser CpG island track](https://hgdownload.soe.ucsc.edu/hubs/GCF/015/227/805/GCF_015227805.1/bbi/GCF_015227805.1_bHirRus1.pri.v2.cpgIslandExt.bb)</br>
Only chromosomes coordinates were mantained.</br>

Genes coordinates were taken from [GenomicFeatures analysis](https://github.com/SwallowGenomics/BarnSwallow/tree/main/Analyses/GenomicFeatures) and merged with `bedtools merge`</br>

Repeats coordinates were generated as in [FIGURE1 panel c](https://github.com/SwallowGenomics/BarnSwallow/tree/main/Plots%20and%20figures/FIGURE1/panel_C)



