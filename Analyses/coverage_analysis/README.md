# coverage_analysis

This short script will first split the "per-base" genome coverage output from bedtools genomecov into files relative to the different chromosomes; next, it will compute the average coverage of each chromosome every 500 bp window and generate a bed file reporting this information.  <br />
first input: file with the list of chromosomes name as reported in the reference fasta.fai index  <br />
second input: "per-base" genome coverage output from bedtools genomecov   <br />
usage: bash avg_coverage.sh chrs per_base_genomecov_output 