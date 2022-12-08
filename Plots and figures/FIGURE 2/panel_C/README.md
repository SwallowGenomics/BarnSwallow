## FIGURE1 panel C

This folder contains the R scripts and example files for Figure 2 panels c.

The package used for plotting is [circlize](https://github.com/jokergoo/circlize)

Tracks used for the figure are provided, exept for the repeats density track wich was too big.

Tracks were obtained as following:

1. generate file with 200 kbp windows from the genome

`bedtools makewindows -b CHR_coords.bed -w 200000 > binned_genome_200k.bed`

2. pacBio coverage

`minimap2 -t 32 -ax map-pb -L bHirRus1_primary.fasta PacBio.subreads.fastq | samtools view -b -S -@ 32 -o bHirRus1_primary_Pacbio_alignment.bam`</br>
`samtools sort -@ 32 bHirRus1_primary_Pacbio_alignment.bam -o bHirRus1_primary_Pacbio_alignment_sort.bam` </br>
`samtools index -@ 32 bHirRus1_primary_Pacbio_alignment_sort.bam`</br>
`mosdepth --no-per-base -by 200000 Pacbio_cov_200k bHirRus1_primary_Pacbio_alignment_sort.bam`</br>
Only chromosomes coordinates were mantained.

3. GC content

`bedtools nuc -fi bHirRus1_primary.fasta -bed binned_genome_200k.bed  > GC_CONTENT_200k_bedtools_nuc.bed`</br>
Only chromosomes coordinates were mantained.</br>
`awk '{print $1, $2, $3, $5*100}' GC_CONTENT_200k_bedtools_nuc.bed > GC_CONTENT_200k.bed`</br>
Replace spaces with tabs and remove first row (header).


4. CpG islands

CpG islands were downloaded from the [UCSC browser CpG island track](https://hgdownload.soe.ucsc.edu/hubs/GCF/015/227/805/GCF_015227805.1/bbi/GCF_015227805.1_bHirRus1.pri.v2.cpgIslandExt.bb)</br>
Only chromosomes coordinates were mantained.</br>

`bedtools intersect -c -a binned_genome_200k.bed -b CpG_ISLANDS_CHR.bed > CpG_ISLANDS_density_200k_n.bed` 

5. Repeats

Use [generate_masked_ranges_mod.py](https://github.com/SwallowGenomics/BarnSwallow/blob/main/Plots%20and%20figures/FIGURE1/panel_C/generate_masked_ranges_mod.py) modified from [generate_masked_ranges_mod.py](https://www.danielecook.com/generate-a-bedfile-of-masked-ranges-a-fasta-file/) to obtain soft-masked coordinates from the genomes. </br>
`python2 generate_masked_ranges_mod.py bHirRus1_primary_masked.fasta > REPEATS.bed` </br>
`bedtools intersect -wao -a binned_genome_200k_3839.bed -b REPEATS.bed > REPEATS_200k_1.bed` </br>
`awk '{print $1, $2, $3, $7}' REPEATS_200k_1.bed > REPEATS_200k.bed` </br>
Replace spaces with tabs.

6. Genes

Use the genes from GenomicFeatures with merged coordinates (`bedtools merge`).</br>
Only chromosomes coordinates were mantained.</br>

`bedtools intersect -c -a binned_genome_200k.bed -b genes_CHRs_merged.bed > genes_200k.bed`

7. phyloP accelerated

`bedtools intersect -c -a binned_genome_200k.bed -b PhyloP_10bp_ACC_NO_GAPS.bed > PhyloP_10bp_FDR_ACC_200k.bed`

8. phyloP conserved

`bedtools intersect -c -a binned_genome_200k.bed -b PhyloP_10bp_CONS_NO_GAPS.bed > PhyloP_10bp_FDR_CONS_200k.bed`

9. phastCons CEs

`bedtools intersect -c -a binned_genome_200k.bed -b most-conserved_NO_GAPS_FINAL.bed > phastCons_CEs_200k.bed`

10. HAL cov 

The coverage was calculated from the [Cactus alignment](https://github.com/SwallowGenomics/BarnSwallow/tree/main/Analyses/Cactus_alignment) for barn swallow chromosome separately. Example on chr.1:</br>
`halAlignmentDepth --noAncestors --step 200000 --targetGenomes Camarhynchus_parvulus,Gallus_gallus,Lonchura_striata,Molothrus_ater,Motacilla_alba,Taeniopygia_guttata,Passer_domesticus --outWiggle coverage_10_genomes_chr1_200k.wig --refSequence chr1 10_genomes.hal Hirundo_rustica`</br>
`wig2bed < coverage_10_genomes_chr1_200k.wig > coverage_10_genomes_chr1_200k.bed`    </br>
`cat *bed > coverage_10_genomes_200k.bed`</br>