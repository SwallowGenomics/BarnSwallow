#!/bin/bash


input=$1

readarray -t a < "${input}"

#download all data from NCBI

for f in "${a[@]}" ; do

        echo "\nprocessing: ${f}"

    if ! [[ -e ${f}/${f}.fastq ]] ; then   

		fasterq-dump ${f} -O ${f} -e 12
	
	else
	
		echo "${f}/${f}.fastq already exists"
	
	fi 
	
#perform quality control

	if ! [[ -e ${f}/QC ]]; then
	
		mkdir ${f}/QC
    
	else
	
		echo "${f}/QC already exists"
		
	fi

	if ! [[ -e ${f}/QC/${f}_1_fastqc.html ]] ; then
	
		fastqc ${f}/${f}_1.fastq -o ${f}/QC -t 32

    else
	
		echo "${f}/QC/${f}_1_fastqc.html already exists"
		
	fi
	
	if ! [[ -e ${f}/QC/${f}_2_fastqc.html ]]; then
	
		fastqc ${f}/${f}_2.fastq -o ${f}/QC -t 32
	
	else 
	
		echo "${f}/QC/${f}_2_fastqc.html already exists"
		
	fi
	
#If adapters contamination was detected after quality control, perform adapters trimming	
	
	if ! [[ -e ${f}/trimmed_${f}_1.fastq ]] || ! [[ -e ${f}/trimmed_${f}_2.fastq ]] ; then
	
		bbduk.sh in1=${f}/${f}_1.fastq in2=${f}/${f}_2.fastq out1=${f}/trimmed_${f}_1.fastq out2=${f}/trimmed_${f}_2.fastq ref=.../bin/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo stats="${f}/${f}_bbduk.out"
		
	else
	
	
		echo "either ${f}/trimmed_${f}_1.fastq or ${f}/trimmed_${f}_2.fastq already exists"
	
	fi

	if  [[ -e ${f}/trimmed_${f}_1.fastq ]] || [[ -e ${f}/trimmed_${f}_2.fastq ]] ; then
	
	
	rm ${f}/${f}*.fastq
		
	else 	

		echo "could not detect either ${f}/trimmed_${f}_1.fastq or ${f}/trimmed_${f}_2.fastq"
	
	fi

#re-perform quality control after adapters trimming		
		
	if ! [[ -e ${f}/QC2 ]] ; then	
	
		mkdir ${f}/QC2
		
	else	
	
		echo "${f}/QC2 already exists"
		
	fi
	
	if ! [[ -e ${f}/QC2/trimmed_${f}_1_fastqc.html ]]; then
	
		fastqc ${f}/trimmed_${f}_1.fastq -o ${f}/QC2 -t 32
		
	else 
		
		echo "${f}/QC2/trimmed_${f}_1_fastqc.html already exists"
	
	fi	
	
	if ! [[ -e ${f}/QC2/trimmed_${f}_2_fastqc.html ]]; then
	
		fastqc ${f}/trimmed_${f}_2.fastq -o ${f}/QC2 -t 32
	
	else	
	
		echo "${f}/QC2/trimmed_${f}_2_fastqc.html already exists"
		
	fi
	
#perform reads alignment
	
	if ! [[ -e ${f}/alignment ]]; then

		mkdir ${f}/alignment

	else
	
		echo "${f}/alignment already exists"
		
	fi

	if ! [[ -e ${f}/alignment/${f}.bam ]]; then
	
	
		bowtie2 -x ref_index -1 ${f}/trimmed_${f}_1.fastq -2 ${f}/trimmed_${f}_2.fastq -p 2 --rg "SM:${f}" --rg-id "${f}" | samtools view -@ 16 -bS - > ${f}/alignment/${f}.bam
	
	else
	
		echo "${f}/alignment/${f}.bam already exists"
		
	fi

	if [[ -s ${f}/alignment/${f}.bam ]]; then
	
		rm ${f}/trimmed_${f}*.fastq
		
	else 
	
		echo "${f}/alignment/${f}.bam does not exist"
		
	fi

#alignment sorting and stats; remove unmapped reads

	if ! [[ -e ${f}/alignment/${f}_sorted.bam ]]; then

		samtools sort ${f}/alignment/${f}.bam -o ${f}/alignment/${f}_sorted.bam -@ 2
		
		
		samtools stats ${f}/alignment/${f}_sorted.bam 1> ${f}/alignment/${f}_stats.out


		if [[ -s ${f}/alignment/${f}_sorted.bam ]]; then
		
			rm ${f}/alignment/${f}.bam
			
		else 

			echo "${f}/alignment/${f}_sorted.bam does not exist"
		
		fi

		samtools view -bF4 ${f}/alignment/${f}_sorted.bam > ${f}/alignment/${f}_sorted_mapped.bam

	rm ${f}/alignment/${f}_sorted.bam
	
	mv ${f}/alignment/${f}_sorted_mapped.bam ${f}/alignment/${f}_sorted.bam
	
	
	else 
	
		echo "${f}/alignment/${f}_sorted.bam already exists"
		
	fi


#index bam file

if ! [[ -e ${f}/alignment/${f}_sorted.bam.bai ]]; then
		
        samtools index ${f}/alignment/${f}_sorted.bam
          
    else
	
		echo "${f}/alignment/${f}_sorted.bam.bai already exists"
		
	fi

#Remove duplicated reads (this step was not not applied to ddRAD data)

    if ! [[ -e ${f}/alignment/${f}_sorted_rmdup.bam ]]; then
		

       picard MarkDuplicates -REMOVE_DUPLICATES true -I ${f}/alignment/${f}_sorted.bam -O ${f}/alignment/${f}_sorted_rmdup.bam -M ${f}/alignment/${f}_sorted_rmdup_metrics.txt

    else
	
		echo "${f}/alignment/${f}_sorted_rmdup.bam already exists"
		
	fi

#Clip overlaps between paired reads (this step was applied only to paired-end data)

    if ! [[ -e ${f}/alignment/${f}_sorted_rmdup_clipped.bam ]]; then
	
	bam clipOverlap --in ${f}/alignment/${f}_sorted_rmdup.bam --out ${f}/alignment/${f}_sorted_rmdup_clipped.bam --stats 2> ${f}/alignment/${f}_clipstats.out
		
	else
	
		echo "${f}/alignment/${f}_sorted_rmdup_clipped.bam already exists"
		
	fi

#remove intermediate files


	if [[ -s ${f}/alignment/${f}_sorted_rmdup_clipped.bam ]]; then
	
		rm ${f}/alignment/${f}_sorted.bam
		rm ${f}/alignment/${f}_sorted_rmdup.bam

	else
	
		echo "${f}/alignment/${f}_sorted_rmdup_clipped.bam does not exist"
	
	fi


#re-index bam file

    if ! [[ -e ${f}/alignment/${f}_sorted_rmdup_clipped.bam.bai ]]; then
		

        samtools index ${f}/alignment/${f}_sorted_rmdup_clipped.bam
  
    else
	
		echo "${f}/alignment/${f}_sorted_rmdup_clipped.bam.bai already exists"
		
	fi


done