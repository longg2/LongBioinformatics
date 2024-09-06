#! /usr/bin/env bash

STARMapping(){ # STAR Mapping.  Automatically determines if merged or paired.
	if [ "$merged" != "NA" ]; then # If I found a merged file
	#if [ -v merged ]; then # If I found a merged file
		STAR --runThreadN $ncores --genomeDir $ref --readFilesIn $merged \
			--readFilesCommand gunzip -c --outFileNamePrefix tmp/tmpM

		samtools view -b -h -F 4 -m $len -q $qual -U tmpMBad.bam tmp/tmpMAligned.out.sam |\
			samtools sort -	> tmpM.bam  
		
		samtools fastq tmpMBad.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz
	fi
	
	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If I found a single read
		STAR --runThreadN $ncores --genomeDir $ref --readFilesIn $r1 \
			--readFilesCommand gunzip -c --outFileNamePrefix tmp/tmpS

		samtools view -b -h -F 4 -m $len -q $qual -U tmpSBad.bam tmp/tmpSAligned.out.sam |\
			samtools sort -	> tmpS.bam  
	
		samtools fastq tmpS.bam | gzip > ${out}UnmappedReads/${sample}_Single.fastq.gz

	fi

	if [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then # If I found a Paired set
	#if [ -v r1 ]; then # If we have paired reads
		STAR --runThreadN $ncores --genomeDir $ref --readFilesIn $r1 $r2 \
			--readFilesCommand gunzip -c --outFileNamePrefix tmp/tmpP
	
		samtools view -b -h -F 4 -m $len -q $qual -U tmpPBad.bam tmp/tmpPAligned.out.sam |\
		samtools sort - > tmpP.bam 
		
		samtools fastq -c 6 tmpPBad.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null 

	
	fi
	
	# Merging the files if needed.  Otherwise, we can ignore and simply mv it
	#if [ -v merged ] && [ -v r1 ]; then
	if [ "$merged" != "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
		#echo "MERGED and PAIRED"
		samtools merge -r -f tmpMerged.bam tmpM.bam tmpP.bam 
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpMerged.bam
		samtools merge -f ${out}UnmappedBam/$sample.bam tmpMBad.bam tmpPBad.bam
		#samtools merge -f ${out}MQ0Bam/$sample.bam tmpMQual0.bam tmpPQual0.bam
		rm tmpP* tmpM*
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		#echo "MERGED"
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpM.bam
		#mv tmpM.bam ${out}MappedReads/$sample.bam 
		mv tmpMBad.bam ${out}UnmappedBam/$sample.bam 
		#mv tmpMQual0.bam ${out}MQ0Bam/$sample.bam 
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpS.bam
		#mv tmpS.bam  ${out}MappedReads/$sample.bam 
		mv tmpSBad.bam ${out}UnmappedBam/$sample.bam 
		#mv tmpSQual0.bam ${out}MQ0Bam/$sample.bam 
	else
		#echo "PAIRED"
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpP.bam
		#mv tmpP.bam ${out}MappedReads/$sample.bam
		mv tmpPBad.bam ${out}UnmappedBam/$sample.bam 
		#mv tmpPQual0.bam ${out}MQ0Bam/$sample.bam 
	fi

	# Deleting temporary files
	rm tmp*

}

alnMapping(){ # BWA aln Mapping.  Automatically determines if merged or paired.
	if [ "$merged" != "NA" ]; then # If I found a merged file
	#if [ -v merged ]; then # If I found a merged file
		bwa aln -o 2 -n 0.01 -l 16500 $ref $merged -t $ncores > tmp.sai
	
		bwa samse -n 10 $ref tmp.sai $merged |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmpMBad.bam |\
			samtools sort -	> tmpM.bam  
	
		# Now to figure out if Multimappers are to be extracted
		if [ "$multi" != "TRUE" ]; then
			samtools fastq tmpMBad.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz

		else
			# Apr 20 2021 -- Want to separate the poor hits from the reads with multiple hits
			samtools view -h -F 4 -m $len -q 0 tmpMBad.bam > tmpMQual0.bam 
		#	samtools fastq -f 4 tmpMBad.bam | gzip > tmpUnmapped.fastq.gz #${out}UnmappedReads/${sample}.fastq.gz # All the unmapped Reads
		#	samtools view -F 4 tmpMBad.bam | awk '{if($5 > 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip -c | cat tmpUnmapped.fastq.gz -  > ${out}UnmappedReads/${sample}.fastq.gz # The reads which were poor matches 
		#	samtools view -m $len -F 4 tmpMBad.bam | awk '{if($5 == 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip > ${out}Qual0Reads/${sample}.fastq.gz # The Multimappers
		fi
	fi
	
	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If I found a merged file
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r1 -t $ncores > tmp.sai
		bwa samse -n 10 $ref tmp.sai $r1 |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmpSBad.bam |\
			samtools sort -	> tmpS.bam  
	
		samtools fastq tmpS.bam | gzip > ${out}UnmappedReads/${sample}_Single.fastq.gz

		if [ "$multi" != "TRUE" ]; then
			samtools fastq tmpSBad.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz

		else
			# Apr 20 2021 -- Want to separate the poor hits from the reads with multiple hits
			samtools view -h -F 4 -m $len -q 0 tmpSBad.bam > tmpSQual0.bam 
			#samtools fastq -f 4 tmpSBad.bam | gzip > tmpUnmapped.fastq.gz #${out}UnmappedReads/${sample}.fastq.gz # All the unmapped Reads
			#samtools view -F 4 tmpSBad.bam | awk '{if($5 > 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip -c | cat tmpUnmapped.fastq.gz -  > ${out}UnmappedReads/${sample}_r1.fastq.gz # The reads which were poor matches 
			#samtools view -m $len -F 4 tmpSBad.bam | awk '{if($5 == 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip > ${out}Qual0Reads/${sample}_r1.fastq.gz # The Multimappers
		fi
	fi

	if [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then # If I found a Paired set
	#if [ -v r1 ]; then # If we have paired reads
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r1 -t $ncores  > tmp1.sai # First index
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r2 -t $ncores  > tmp2.sai # Second Index
	
		bwa sampe -n 10 $ref tmp1.sai tmp2.sai $r1 $r2  |\
		       	samtools view -b -h -m $len -q $qual -U tmpPBad.bam |\
			samtools sort -n - | samtools fixmate -r -m - - |\
			samtools view -b -h -F 2048 |\
	       		samtools sort - > tmpP.bam

			# 2024-09-06: Old way of handling paired reads. Wasn't really accounting for PEs the right way. The above fixes it since I've included fixmate
			#samtools view -b -h -F 4 -m $len -q $qual -U tmpPBad.bam |\
		       	#samtools sort - > tmpP.bam 
		
		if [ "$multi" != "TRUE" ]; then
			samtools fastq -c 6 tmpPBad.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null 

		else
			echo "PANIC"
			samtools view -h -F 4 -m $len -q 0 tmpPBad.bam > tmpPQual0.bam 
		fi
	
	fi
	
	# Merging the files if needed.  Otherwise, we can ignore and simply mv it
	#if [ -v merged ] && [ -v r1 ]; then
	if [ "$merged" != "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
		#echo "MERGED and PAIRED"
		samtools merge -r -f tmpMerged.bam tmpM.bam tmpP.bam 
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpMerged.bam
		samtools merge -f ${out}UnmappedBam/$sample.bam tmpMBad.bam tmpPBad.bam
		if [ "$multi" == "TRUE" ]; then
			samtools merge -f ${out}MQ0Bam/$sample.bam tmpMQual0.bam tmpPQual0.bam
		fi
		rm tmpP* tmpM*
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		#echo "MERGED"
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpM.bam
		#mv tmpM.bam ${out}MappedReads/$sample.bam 
		mv tmpMBad.bam ${out}UnmappedBam/$sample.bam 

		if [ "$multi" == "TRUE" ]; then
			mv tmpMQual0.bam ${out}MQ0Bam/$sample.bam 
		fi
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpS.bam
		#mv tmpS.bam  ${out}MappedReads/$sample.bam 
		mv tmpSBad.bam ${out}UnmappedBam/$sample.bam 
		mv tmpSQual0.bam ${out}MQ0Bam/$sample.bam 
	else
		#echo "PAIRED"
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpP.bam
		#mv tmpP.bam ${out}MappedReads/$sample.bam
		mv tmpPBad.bam ${out}UnmappedBam/$sample.bam 
		if [ "$multi" == "TRUE" ]; then
			mv tmpPQual0.bam ${out}MQ0Bam/$sample.bam 
		fi
	fi

	# Deleting temporary files
	rm tmp*
}

bowtie2Mapping(){ # Bowtie2 Mapping. Using local mapping with the very sensitive flag....
	if [ "$merged" != "NA" ]; then # If I found a merged file
	#if [ -v merged ]; then # If I found a merged file
		bowtie2 --local --very-sensitive-local -U $merged --threads $ncores -x $ref |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmpMBad.bam |\
		       	samtools sort - > tmpM.bam 
	
		# Now to figure out if Multimappers are to be extracted
		samtools fastq tmpMBad.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz
	fi
	
	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If I found a merged file
		bowtie2 --local --very-sensitive-local -U $r1 --threads $ncores -x $ref |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmpSBad.bam |\
		       	samtools sort - > tmpS.bam 
	
		samtools fastq tmpS.bam | gzip > ${out}UnmappedReads/${sample}_Single.fastq.gz

		samtools fastq tmpSBad.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz
	fi

	if [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then # If I found a Paired set
	#if [ -v r1 ]; then # If we have paired reads
		bowtie2 --local --very-sensitive-local -1 $r1 -2 $r2 --threads $ncores -x $ref |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmpPBad.bam |\
		       	samtools sort - > tmpP.bam 
	
		samtools fastq -c 6 tmpPBad.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null 
	fi
	
	# Merging the files if needed.  Otherwise, we can ignore and simply mv it
	#if [ -v merged ] && [ -v r1 ]; then
	if [ "$merged" != "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
		#echo "MERGED and PAIRED"
		samtools merge -r -f tmpMerged.bam tmpM.bam tmpP.bam 
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpMerged.bam
		samtools merge -f ${out}UnmappedBam/$sample.bam tmpMBad.bam tmpPBad.bam
		rm tmpP* tmpM*
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		#echo "MERGED"
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpM.bam
		#mv tmpM.bam ${out}MappedReads/$sample.bam 
		mv tmpMBad.bam ${out}UnmappedBam/$sample.bam 
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpS.bam
		#mv tmpS.bam  ${out}MappedReads/$sample.bam 
		mv tmpSBad.bam ${out}UnmappedBam/$sample.bam 
	else
		#echo "PAIRED"
		samtools addreplacerg -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpP.bam
		#mv tmpP.bam ${out}MappedReads/$sample.bam
		mv tmpPBad.bam ${out}UnmappedBam/$sample.bam 
	fi

	# Deleting temporary files
	rm tmp.bam tmp*.sai

}

samtoolsDeduplication(){ # Deduplication suing samtools markdup
	# Setting up the functions
	local inFile=$1
	local outFile=$2
	local statsFolder=$3
	local sample=$(basename $inFile .bam)

	# We need to resort the reads by name and run fix mate first
	samtools sort -n $inFile | samtools fixmate -m - - |\
		samtools sort - |\
		samtools markdup -r -S -f $statsFolder/$sample.tab - $outFile

	# The other two lines deduplicate the reads *and* resort bam file by co-ordinates
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
		bwa mem $ref $r1 $r2 -t $ncores |\
		       	samtools view -b -h -m $len -q $qual -U tmp.bam |\
			samtools sort -n - | samtools fixmate -r -m - - |\
			samtools view -b -h -F 2048 |\
	       		samtools sort - > tmpP.bam
			# 2024-09-06: Old way of handling paired reads. Wasn't really accounting for PEs the right way. The above fixes it since I've included fixmate
		#bwa mem $ref $r1 $r2 -t $ncores | samtools view -b -h -f 3 -m $len -q $qual -U tmp.bam |\
	       	#	samtools sort - > tmpP.bam
	
		samtools fastq -c 6 tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null
	
	fi
	
	# Merging the files if needed.  Otherwise, we can ignore and simply mv it
	#if [ -v merged ] && [ -v r1 ]; then
	if [ "$merged" != "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
		#samtools merge -r -f ${out}MappedReads/$sample.bam tmpM.bam tmpP.bam 
		samtools merge -r -f tmpMerged.bam tmpM.bam tmpP.bam 
		samtools addreplacerg -@ $((ncores - 1)) -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpMerged.bam
		#rm tmpP.bam tmpM.bam tmpMerged.bam
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
	#	mv tmpSingle.bam  ${out}MappedReads/$sample.bam 
		samtools addreplacerg -@ $((ncores - 1)) -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpSingle.bam
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		#mv tmpM.bam ${out}MappedReads/$sample.bam 
		samtools addreplacerg -@ $((ncores - 1)) -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpM.bam
	else
		#mv tmpP.bam ${out}MappedReads/$sample.bam
		samtools addreplacerg -@ $((ncores - 1)) -r ID:$sample -r SM:$sample -o ${out}MappedReads/$sample.bam tmpP.bam
	fi

	# Deleting temporary files
	rm tmp*.bam

}
vgMapping(){
	# Merged Reads
	if [ "$merged" != "NA" ]; then # If I found a merged file
		vg map -t $ncores -d $ref -f $merged -k 15 -w 1024 --surject-to bam | samtools view -b -h -F 4 -q $qual -U tmp.bam |\
			samtools sort - > tmpM.bam
		samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz
	fi

	# Paired Reads
	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then
	#if [ -v r1 ]; then # If we have paired reads
		vg map -t $ncores -d $ref -f $r1 -f $r2 -k 15 -w 1024 --surject-to bam | samtools view -b -h -F 4 -q $qual -U tmp.bam |\
			samtools sort - > tmpP.bam
	
		samtools fastq -c 6 tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null
	
	fi

	# Single End Sequencing
	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then
		vg map -t $ncores -d $ref -f $r1 -k 15 -w 1024 --surject-to bam | samtools view -b -h -F 4 -q $qual -U tmp.bam |\
			samtools sort - > tmpM.bam
		samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz
	
	fi



}
