## Extract genomic features from gff3 annotation file

This folder contains the scripts used to obtain genes, 5'UTRs, 3'UTRs, exons, introns, CDS and promoter regions from a gff3 file. <br />

[pre-processing.txt](pre-processing.txt) was used to remove tRNAs, pseudogenes and C/V_gene_segments from the annotation. In the gff3 file from NCBI, transcript ids are not assigned to these gene type.

[create_regions_from_gencode_modified.R](create_regions_from_gencode_modified.R) is the script used to generate diffrenet bed files from the annotation. <br />
The script was modified from: [create_regions_from_gencode.R](https://github.com/saketkc/gencode_regions/blob/master/create_regions_from_gencode.R)

[intergenic.txt](intergenic.txt) Commands used to extract intergenic coordinates.