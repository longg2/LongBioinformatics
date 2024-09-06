#! /usr/bin/env bash

snippyParallel(){ # Want this to figure out what kind of Snippy run its doing

	# Reading in the variables
	local seqFile=$1
	local mincov=$2
	local reference=$3
	local out=$4

	# Now I need to figure out which version of the command below I'm running (i.e. ctgs/se/contigs)
	local fileExtension=$(echo $seqFile | sed -e 's/\.gz//g' -e 's/.*\.//g')
	local sampleName=$(basename $seqFile | sed -e 's/\.fa.*//I' -e 's/\.bam//I' -e 's/\.fna.*//I' -e 's/\.fq.*//I')

	#printf "SeqFile:\t${seqFile}\nExtension:\t${fileExtension}\nSampleName:\t${sampleName}\nReference:\t${reference}\n"
	# This is what will run the correct version of the snippy command
	if echo $fileExtension | grep -P -i -q "bam"; then
	#	echo "BAM"
		snippy --bam $seqFile --reference $reference --outdir ${out}/$sampleName --prefix $sampleName --mincov $mincov --mapqual 30 --cpus 8 --force 
	elif echo $fileExtension | grep -P -i -q "fastq|fq"; then
#		if echo $seqFile | grep -P -i -q "_R1.f|_1.f"; then
#	#		echo "Paired End"
#			snippy --R1 $seqFile --R2 ${seqFile/1./2.} --reference $reference --outdir ${out}/$sampleName --prefix $sampleName --mincov $mincov --mapqual 30 --cpus 8 --force --quiet
		#elif echo $seqFile | grep -v -P -i -q "_R2.f|_R1.f|_1.f|_2.f"; then
	#		echo "Single Ended (Merged?)"
		snippy --se $seqFile --reference $reference --outdir ${out}/$sampleName --prefix $sampleName --mincov $mincov --mapqual 30 --cpus 8 --force --quiet
		#fi
	elif echo $fileExtension | grep -P -i -q "fasta|fna|fa"; then
	#	echo "Contigs/Genomes"
		snippy --ctgs $seqFile --reference $reference --outdir ${out}/$sampleName --prefix $sampleName --mincov $mincov --mapqual 30 --cpus 8 --force --quiet
	else
		echo "I don't know how to deal with this!"
		exit 1
	fi

}
