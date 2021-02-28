#! /usr/bin/env sh

# BEGIN FUNCTION (MIGHT NOT BE NEEDED!)
# Testing if they're gzipped (Should be a function unto itself)
#if [[ $file =~ \.t?gz$ ]];
#then
#	echo "$sample was gzipped"
#	gunzip -c $in/$sample.fastq > merge.fastq
#	gunzip -c $file > r1.fastq
#	gunzip -c $r2 > r2.fastq
#else
#	cp $in/$sample.fastq merge.fastq
#	cp $file r1.fastq
#	cp $r2 r2.fastq
#fi
# END FUNCTION

#############################
# Performing the alignments #
#############################
# Should be a function itself
# BEGIN FUNCTION
alnMapping(){
	if [ -v merged ]; then # If I found a merged file
		bwa aln -o 2 -n 0.01 -l 16500 $ref $merged -t $ncores 2>/dev/null > tmp.sai
	
		bwa samse $ref tmp.sai $merged 2>/dev/null |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam|\
			samtools sort -	> tmpM.bam
	
		samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${sample}.fastq
	fi
	
	if [ -v r1 ]; then # If we have paired reads
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r1 -t $ncores 2>/dev/null > tmp1.sai # First index
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r2 -t $ncores 2>/dev/null > tmp2.sai # Second Index
	
		bwa sampe $ref tmp1.sai tmp2.sai $r1 $r2 2>/dev/null |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\
		       	samtools sort - > tmpP.bam
		
		samtools fastq -c 6 tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null
	
	fi
	
	# Merging the files if needed.  Otherwise, we can ignore and simply mv it
	if [ -v merged ] && [ -v r1 ]; then
		samtools merge -f ${out}MappedReads/$sample.bam tmpM.bam tmpP.bam
	elif [ -v merged ] && [ -z ${r1+x} ]; then
		mv tmpM.bam ${out}MappedReads/$sample.bam
	else
		mv tmpP.bam ${out}MappedReads/$sample.bam
	fi

	# Deleting temporary files
	rm tmpP.bam tmpM.bam

}
# END FUNCTION

# QoL progress
printf "$sample has been filtered ($count/$total)\n"
count=$(echo "$count + 1" | bc)

