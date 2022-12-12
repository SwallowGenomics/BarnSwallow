## karyoploteR figures

The R files present in this folder contain the commands used to generate SNP density plots with the [karyoploteR](https://bernatgel.github.io/karyoploter_tutorial/) package. Macrochromosomes and intermediate chromosomes are plotted together, while microchromosomes are plotted separately (see Figure S5). <br />
The [chrom_list](https://github.com/SwallowGenomics/BarnSwallow/blob/main/Plots%20and%20figures/FIGURE%203/panel_B/macro_interm_chromosomes/chrom_list.ls) files contain chromosome sizes and are requested by the karyoploteR package in order to plot ideograms using custom genomes. <br />
Files required to plot the additional tracks of the figure (assembly gaps, repeats, GC content, Pacbio reads coverage) have been provided for chromosome 1 and 14 as example files. <br />
The `chr_snp_pos` files contain SNP positions for all sequencing technologies analyzed, extracted directly from the relative vcf files with the bcftools query command.