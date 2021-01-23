#!/bin/bash

input=$1

readarray -t a < "${input}"

#download all data

for f in "${a[@]}"; do
	
	echo "processing: ${f}"
	fasterq-dump ${f} -O ${f} -e 12

done

#find all files
find ./* -name "*.fastq" > files.ls

#find dir of all files
find ./* -name '*.fastq' -exec dirname {} \; > dirs.ls

#generate QC dirs
cat dirs.ls | awk -F' ' '{a["mkdir "$1"/QC"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 1

#generate QC
paste files.ls dirs.ls | awk -F' ' '{a["fastqc "$1" -o "$2"/QC"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 10  &> fastqc.out &

#generate multiQC
current_path=$(pwd)
multiqc ${current_path}/*/QC/*_fastqc.zip

mv multiqc_data multiqc_data_beforeQC
mv multiqc_report.html multiqc_report_beforeQC.html

#generate list of fw and rv reads
find ./* -name "*_1.fastq.gz" > fw.list
find ./* -name "*_2.fastq.gz" > rv.list

#trim adapters from fw and rv reads
paste fw.list rv.list | awk -F' ' '{a["bbduk.sh in1="$1" in2="$2" out1="$1"_clean out2="$2"_clean ref=/home/gformenti/bin/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 10 &> bbdukPE.out &

#trim adapters from unique files
find ./* -name "SRR[0-9][0-9][0-9][0-9][0-9][0-9][0-9].fastq" | awk -F' ' '{a["bbduk.sh in="$1" out="$1"_clean ref=/home/gformenti/bin/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 10 &> bbdukS.out &

#change file names
find ./* -name '*fastq_clean' -exec bash -c 'echo mv $0 ${0/\.fastq_clean/\_clean\.fastq}' {} \; | parallel --gnu -j 1

#remove untrimmed fastq
rm S*/*[0-9].fastq

#find all cleaned files
find ./* -name "*_clean.fastq" > files_clean.ls

#find dir of all cleaned files
find ./* -name '*_clean.fastq' -exec dirname {} \; > dirs_clean.ls

#generate QC dirs
cat dirs_clean.ls | awk -F' ' '{a["mkdir "$1"/QC2"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 1

#generate QC
paste files_clean.ls dirs_clean.ls | awk -F' ' '{a["fastqc "$1" -o "$2"/QC2"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 10 &> fastqc2.out &

#generate multiQC
current_path=$(pwd)
multiqc ${current_path}/*/QC2/*_fastqc.zip

mv multiqc_data multiqc_data_afterQC
mv multiqc_report.html multiqc_report_afterQC.html

#find all cleaned unique files
find ./* -name "*[0-9]_clean.fastq" > files_clean_u.ls

#find all cleaned fw files
find ./* -name "*_1_clean.fastq" > files_clean_fw.ls

#find dir of all cleaned rv files
find ./* -name "*_2_clean.fastq" > files_clean_rv.ls

paste files_clean_fw.ls files_clean_rv.ls > files_clean_split.ls

arr_fw=( $(awk '{print $1}' files_clean_split.ls) )
arr_rv=( $(awk '{print $2}' files_clean_split.ls) )

#export bowtie2 index location
export BOWTIE2_INDEXES=/home/gformenti/hirundo/

while [  $COUNTER -lt ${#arr_fw[@]} ]
do
	NAME=$(basename ${arr_fw[$COUNTER]})
	FW=${arr_fw[$COUNTER]}
	RV=${arr_rv[$COUNTER]}
	#generate alignments directory
	mkdir ${NAME}/alignment
	#align cleaned files
	bowtie2 -x ../../../hirundoN -1 ${FW} -2 ${RV} -p 12 | samtools view -bSF4 - > "${NAME}/alignment/${NAME}.bam"

done

#find all bam
ls ./*/*/*bam > all_bam.ls

#sort all bam
cat all_bam.ls | awk -F' ' '{a["samtools sort "$1" -o "$1"_sorted; rm "$1]}END{for(i in a) print i a[i]}' | parallel --gnu -j 9

#change file names
find ./* -name '*bam_sorted' -exec bash -c 'echo mv $0 ${0/\.bam_sorted/\_sorted\.bam}' {} \; | parallel --gnu -j 1

#find all bam
ls ./*/*/*bam > all_bam.ls

#if WGS
#remove dup reads
cat all_bam.ls | awk -F' ' '{a["java -jar /home/gformenti/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I="$1" O="$1"_rmdup M="$1"_metrics"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 8 &> MarkDuplicates.out &

#if WGS
#rm bam with dups
cat all_bam.ls | rm

#if WGS
#change file names
find ./* -name '*.bam_rmdup' -exec bash -c 'echo mv $0 ${0/\.bam_rmdup/\_rmdup\.bam}' {} \; | parallel --gnu -j 1
find ./* -name '*.bam_metrics' -exec bash -c 'echo mv $0 ${0/\.bam_metrics/\_metrics\.bam}' {} \; | parallel --gnu -j 1

#find all bam without dups
ls ./*/*/*_rmdup.bam > all_bam_nodup.ls

#index all bam
cat all_bam_nodup.ls | awk -F' ' '{a["samtools index "$1]}END{for(i in a) print i a[i]}' | parallel --gnu -j 8 &> indexing.out &

#get bam names
find ./*/*/ -name '*.bam' -exec basename {} \; | grep -oe "SRR[0-9]*" > bam_names.ls

#if not already added with bowtie2 add RG tag SM and RG_ID
paste all_bam.ls bam_names.ls | awk -F' ' '{a["bamaddrg -b "$1" -s "$2" -r "$2" > "$1"_RG"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 3

#change names
find ./* -name '*bam_sorted' -exec bash -c 'echo mv $0 ${0/\.bam_sorted/\_sorted\.bam}' {} \; | parallel --gnu -j 1

#find all sorted bam
ls ./*/*/*_sorted.bam > all_sorted_bam.ls

#find all sorted bam filenames
cat all_sorted.bam | grep -oE "SRR[0-9]{7}" | uniq > all_sorted_files.ls

#add read group (SM) if not added with bowtie2
paste all_sorted.bam all_sorted_files.ls | awk -F' ' '{a["bamaddrg -b "$1" -s "$2" -r "$2" > "$1"sorted_RG"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 10

#find dir of all bam files
find ./* -name '*_sorted.bam' -exec dirname {} \; > dirs.ls

#get file names without ext
find ./* -name '*_rmdup.bam' -exec basename {} .bam \; > all_bam_nodup_noext.ls

#trim overlaps in reads
paste dirs.ls all_bam_nodup_noext.ls | awk -F' ' '{a["bam clipOverlap --in "$1"/"$2".bam --out "$1"/"$2"_clipped.bam --stats"]}END{for(i in a) print i a[i]}' | parallel --gnu -j 8 &> clipOverlap.out &

#find all clipped bam
find ./* -name '*clipped.bam' > all_clipped_bam.ls

#index all clipped bam
cat all_clipped_bam.ls | awk -F' ' '{a["samtools index "$1]}END{for(i in a) print i a[i]}' | parallel --gnu -j 12

#get masked regions
python generate_masked_ranges.py hirundo_Nmasked.fasta > masked.bed

#get chr length
cat hirundo_Nmasked.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'

#get unmasked regions
bedtools complement -i masked.bed -g chr.len

#get bam lists
ls ../data/american/smith2018/ddRAD/*/alignment/*bam > american_bam.ls
ls ../data/egyptian/smith2018/ddRAD/*/alignment/*bam > egyptian_bam.ls
ls ../data/european/von_ronn2016/*/alignment/*_clipped.bam > european_bam.ls

#concatenate
cat american_bam.ls egyptian_bam.ls european_bam.ls > all_bam.ls

#get sample name lists
ls ../data/american/smith2018/ddRAD/*/alignment/*bam | grep -oe "SRR[0-9]*" | uniq > american.ls
ls ../data/egyptian/smith2018/ddRAD/*/alignment/*bam | grep -oe "SRR[0-9]*" | uniq > egyptian.ls
ls ../data/european/von_ronn2016/*/alignment/*_clipped.bam | grep -oe "SRR[0-9]*" | uniq > european.ls

#add population information to lists
sed 's/$/\tame/' american.ls > american_pl.ls
sed 's/$/\tegy/' egyptian.ls > egyptian_pl.ls
sed 's/$/\teur/' european.ls > european_pl.ls

#concatenate
cat american_pl.ls egyptian_pl.ls european_pl.ls > sample_names_pop.ls

#freebayes vcf
nohup freebayes -f ../reference/hirundo_Nmasked.fasta -L all_bam.ls -t ../reference/unmasked.bed --min-mapping-quality 10 --min-base-quality 20 --min-alternate-count 10 --min-coverage 10 -v ame_egy_eur.vcf --populations sample_names_pop.ls &

#extract europeans
bcftools view -S european.ls ame_egy_eur.vcf > eur.vcf

#filter sites available in >90% of individuals
vcftools --recode --vcf eur.vcf --max-missing 0.9 --out all_notmissing

#filter only SNPs
vcftools --vcf eur_notmissing.recode.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only --out eurSNP

#automatic filter
dDocent_filters eurSNP.vcf _filtered

#convert to plink format
vcftools --vcf eurSNP.vcf --out eurSNP_plink --plink

#add chr indication
cut -f2 eurSNP_plink.map | grep -oE "[^:]*" | awk 'NR % 2 == 1' > chrs.ls

awk '{print $2,$3,$4,$5}' eurSNP_plink.map

nohup freebayes -f ../reference/hirundo_Nmasked.fasta -L all_bam.ls -t ../reference/unmasked.bed --min-mapping-quality 10 --min-base-quality 20 --min-alternate-count 5 --min-coverage 20 -v ame_egy_eur2.vcf --populations sample_names_pop.ls &
nohup freebayes -f ../reference/hirundo_Nmasked.fasta -L all_bam.ls -t ../reference/unmasked.bed --min-mapping-quality 10 --min-base-quality 20 --min-coverage 20 -v ame_egy_eur3.vcf --populations sample_names_pop.ls &

cat eurSNP.recode.vcf | grep -oE "([0-2]/[0-2])" | sort | uniq -c
