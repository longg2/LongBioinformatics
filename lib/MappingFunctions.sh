#! /usr/bin/env bash

alnMapping(){ # BWA aln Mapping.  Automatically determines if merged or paired.
	if [ "$merged" != "NA" ]; then # If I found a merged file
	#if [ -v merged ]; then # If I found a merged file
		bwa aln -o 2 -n 0.01 -l 16500 $ref $merged -t $ncores > tmp.sai
	
		bwa samse $ref tmp.sai $merged |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmpMBad.bam |\
			samtools sort -	> tmpM.bam  
	
		# Now to figure out if Multimappers are to be extracted
		if [ "$multi" != "TRUE" ]; then
			samtools fastq tmpMBad.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz

		else
			# Apr 20 2021 -- Want to separate the poor hits from the reads with multiple hits
			samtools fastq -f 4 tmpMBad.bam | gzip > tmpUnmapped.fastq.gz #${out}UnmappedReads/${sample}.fastq.gz # All the unmapped Reads
			samtools view -F 4 tmpMBad.bam | awk '{if($5 > 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip -c | cat tmpUnmapped.fastq.gz -  > ${out}UnmappedReads/${sample}.fastq.gz # The reads which were poor matches 
			samtools view -m $len -F 4 tmpMBad.bam | awk '{if($5 == 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip > ${out}Qual0Reads/${sample}.fastq.gz # The Multimappers
		fi
	fi
	
	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If I found a merged file
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r1 -t $ncores > tmp.sai
		bwa samse $ref tmp.sai $r1 |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmpSBad.bam |\
			samtools sort -	> tmpS.bam  
	
		samtools fastq tmpS.bam | gzip > ${out}UnmappedReads/${sample}_Single.fastq.gz

		if [ "$multi" != "TRUE" ]; then
			samtools fastq tmpSBad.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz

		else
			# Apr 20 2021 -- Want to separate the poor hits from the reads with multiple hits
			samtools fastq -f 4 tmpSBad.bam | gzip > tmpUnmapped.fastq.gz #${out}UnmappedReads/${sample}.fastq.gz # All the unmapped Reads
			samtools view -F 4 tmpSBad.bam | awk '{if($5 > 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip -c | cat tmpUnmapped.fastq.gz -  > ${out}UnmappedReads/${sample}_r1.fastq.gz # The reads which were poor matches 
			samtools view -m $len -F 4 tmpSBad.bam | awk '{if($5 == 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip > ${out}Qual0Reads/${sample}_r1.fastq.gz # The Multimappers
		fi
	fi

	if [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then # If I found a Paired set
	#if [ -v r1 ]; then # If we have paired reads
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r1 -t $ncores  > tmp1.sai # First index
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r2 -t $ncores  > tmp2.sai # Second Index
	
		bwa sampe $ref tmp1.sai tmp2.sai $r1 $r2  |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmpPBad.bam |\
		       	samtools sort - > tmpP.bam 
		
		if [ "$multi" != "TRUE" ]; then
			samtools fastq -c 6 tmpPBad.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null 

		else
			# Apr 20 2021 -- Want to separate the poor hits from the reads with multiple hits
			samtools fastq -f 4 -c 6 tmpPBad.bam -1 tmpUnmapped_r1.fastq.gz -2 tmpUnmapped_r2.fastq.gz -s /dev/null # All the unmapped Reads
			# R1 reads
			samtools view -F 4 tmpPBad.bam | ~/Scripts/V3Folder/SamQualInttoBin.awk | awk '{if($5 != 0 && substr($2, length($2) - 7, length($2) - 7) == 0){print "@"$1"\n"$10"\n+\n"$11}}'| gzip -c | cat tmpUnmapped_r1.fastq.gz - > ${out}UnmappedReads/${sample}_r1.fastq.gz # The reads which were poor matches
			samtools view -m $len -F 4 tmpPBad.bam | ~/Scripts/V3Folder/SamQualInttoBin.awk | awk '{if($5 == 0 && substr($2, length($2) - 7, length($2) - 7) == 0){print "@"$1"\n"$10"\n+\n"$11}}'| gzip > ${out}Qual0Reads/${sample}_r1.fastq.gz # The Multimappers
			## R2 reads
			samtools view -F 4 tmpPBad.bam | ~/Scripts/V3Folder/SamQualInttoBin.awk | awk '{if($5 != 0 && substr($2, length($2) - 7, length($2) - 7) == 1){print "@"$1"\n"$10"\n+\n"$11}}'| gzip -c | cat tmpUnmapped_r2.fastq.gz - > ${out}UnmappedReads/${sample}_r2.fastq.gz
			samtools view -m $len -F 4 tmpPBad.bam | ~/Scripts/V3Folder/SamQualInttoBin.awk | awk '{if($5 == 0 && substr($2, length($2) - 7, length($2) - 7) == 1){print "@"$1"\n"$10"\n+\n"$11}}'| gzip > ${out}Qual0Reads/${sample}_r2.fastq.gz # The Multimappers
		fi
	
	fi
	
	# Merging the files if needed.  Otherwise, we can ignore and simply mv it
	#if [ -v merged ] && [ -v r1 ]; then
	if [ "$merged" != "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
		#echo "MERGED and PAIRED"
		samtools merge -f ${out}MappedReads/$sample.bam tmpM.bam tmpP.bam
		samtools merge -f ${out}UnmappedBam/$sample.bam tmpMBad.bam tmpPBad.bam
		rm tmpP* tmpM*
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		#echo "MERGED"
		mv tmpM.bam ${out}MappedReads/$sample.bam 
		mv tmpMBad.bam ${out}UnmappedBam/$sample.bam 
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		mv tmpS.bam  ${out}MappedReads/$sample.bam 
		mv tmpSBad.bam ${out}UnmappedBam/$sample.bam 
	else
		#echo "PAIRED"
		mv tmpP.bam ${out}MappedReads/$sample.bam
		mv tmpPBad.bam ${out}UnmappedBam/$sample.bam 
	fi

	# Deleting temporary files
	rm tmp.bam tmp*.sai

}

memMapping(){ # BWA mem Mapping.  Automatically determines if merged or paired.
	if [ "$merged" != "NA" ]; then # If I found a merged file
	#if [ -v merged ]; then # If I found a merged file
		bwa mem $ref $merged -t $ncores |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\
			samtools sort -	> tmpM.bam 
	
		samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz
	fi

	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If I found a merged file
	#if [ -v merged ]; then # If I found a merged file
		bwa mem $ref $r1 -t $ncores |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\
			samtools sort -	> tmpSingle.bam 
	
		samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${sample}_Single.fastq.gz
	fi
	
	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If I found a merged file
	#if [ -v r1 ]; then # If we have paired reads
		bwa mem $ref $r1 $r2 -t $ncores | samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\
	       		samtools sort - > tmpP.bam
	
		samtools fastq -c 6 tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null
	
	fi
	
	# Merging the files if needed.  Otherwise, we can ignore and simply mv it
	#if [ -v merged ] && [ -v r1 ]; then
	if [ "$merged" != "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
		samtools merge -f ${out}MappedReads/$sample.bam tmpM.bam tmpP.bam 
		rm tmpP.bam tmpM.bam 
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		mv tmpSingle.bam  ${out}MappedReads/$sample.bam 
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		mv tmpM.bam ${out}MappedReads/$sample.bam 
	else
		mv tmpP.bam ${out}MappedReads/$sample.bam
	fi

	# Deleting temporary files
	rm tmp.bam

}
