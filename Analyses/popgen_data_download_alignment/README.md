## Pipeline for downloading and processing genomic data from NCBI

This pipeline can be used to download raw genomic data directly from NCBI public database, perform quality control and align them to the reference genome.  <br />
Use the command `bash popgen_pipeline.sh Acclist.txt` to run this pipeline.  <br />
`Acclist.txt` is a text file with all Accession numbers (SRR...) of the genomic data to download from the NCBI public archive. <br />
This pipeline uses the following bioinformatic softwares: [fasterq-dump](https://github.com/ncbi/sra-tools), [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [BBduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/), [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), [Samtools](https://github.com/samtools/samtools), [picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-), [Bam clipOverlap](https://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap)







