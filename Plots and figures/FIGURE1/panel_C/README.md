## FIGURE1 panel C

This folder contains the R scripts and example files for Figure 1 panels c.

The package used for plotting is [circlize](https://github.com/jokergoo/circlize)

Tracks used for the figure are provided, exept for the repeats density track wich was too big.

Tracks were obtained as following:

1. generate file with 200 kbp windows from the genome

`bedtools makewindows -b CHR_coords.bed -w 200000 > binned_genome_200k.bed`

1. pacBio coverage

`minimap2 -t 32 -ax map-pb -L bHirRus1_primary.fasta PacBio.subreads.fastq | samtools view -b -S -@ 32 -o bHirRus1_primary_Pacbio_alignment.bam`
`samtools sort -@ 32 bHirRus1_primary_Pacbio_alignment.bam -o bHirRus1_primary_Pacbio_alignment_sort.bam` 
`samtools index -@ 32 bHirRus1_primary_Pacbio_alignment_sort.bam`
`mosdepth --no-per-base -by 200000 Pacbio_cov_200k bHirRus1_primary_Pacbio_alignment_sort.bam`
Only chromosomes coordinates were mantained.

2. 

3. CpG islands

CpG islands were downloaded from the [UCSC browser CpG island track](https://hgdownload.soe.ucsc.edu/hubs/GCF/015/227/805/GCF_015227805.1/bbi/GCF_015227805.1_bHirRus1.pri.v2.cpgIslandExt.bb)
Only chromosomes coordinates were mantained.

`bedtools intersect -c -a binned_genome_200k.bed -b CpG_ISLANDS_CHR.bed > CpG_ISLANDS_density_200k_n.bed` 

4. Repeats

Use [generate_masked_ranges_mod.py]() modified from [generate_masked_ranges_mod.py](https://www.danielecook.com/generate-a-bedfile-of-masked-ranges-a-fasta-file/) to obtain soft-masked coordinates from the genomes.
`python2 generate_masked_ranges_mod.py bHirRus1_primary_masked.fasta > REPEATS.bed`


5.

6.

7.
