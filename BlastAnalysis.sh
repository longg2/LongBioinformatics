#! /usr/bin/env sh
######################################
### Functions that I'll be calling ###
######################################

GzipDetection(){ # Need to ID Gzipped files and decompress if needed
	local file=$1
	local folderV2=$2
	local name=${file/%.gz*}

	if file $folderV2/$file | grep -1 "compressed"; then
		gunzip -c $folderV2/$file > IntGzip/$name
	else
		cp $folderV2/$file IntGzip/$name
	fi
}
FastaorFastq(){ # Figuring out if fasta or fastq and converting to fasta.
	# All local as I'm avoiding accidental overwrites
	local tmp=$1
	local outFolder=$2
	local firstChar=$(head -c 1 $tmp)
	local name=$(basename $tmp | sed -e "s/\.fastq//" -e "s/\.fq//")
	
	if [ "$firstChar" == ">" ];then
		mv $tmp $outFolder/$(basename $tmp)
	elif [ "$firstChar" == "@" ];then # Need to convert the file to fasta
		sed -n '1~4s/^@/>/p;2~4p' $tmp > $outFolder/$name.fasta
	else
		echo "Skipping $tmp as it's not fasta/fastq" >> FastaFastq.log
	fi
}
export -f FastaorFastq # important as I'm running this in parallel
export -f GzipDetection # important as I'm running this in parallel


#blastCMD() {
#	local in=$1
#
#	if [ "$blast" == "blastn" ]; then
#		blastn -db $db -query $in -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores -perc_identity $Pident -task blastn &
#		pid=$! # Getting the PID of the blast run
#		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well
#
#		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
#
#			# Get the current position of the blast run
#			local curquery=$(tail -1 $blast | cut -f 1)
#			local curline=$(fgrep -n $curquery $query |  cut -f 1 -d ':')
#			local nblines=$(wc -l ${out}/$sample.tab | cut -f 1 -d " ")
#
#
#
#		done
#
#
#	elif [ "$blast" == "blastp" ]; then
#		blastp -db $db -query $in -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores
#	elif [ "$blast" == "blastx" ]; then
#		blastx -db $db -query_gencode 11 -query $in -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores 
#	else
#		echo "Please choose either blastn or blastp"
#		return 1
#	fi
#
#	return 0
#}


declare -r folder=$1
declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
mkdir -p IntGzip

##################
### The Workup ###
##################

# We need to determine if the file is gzipped
echo "Decompressing the files"
parallel -j $ncores --bar "GzipDetection {} $folder" ::: "${files[@]}"

# Now to detect if the file in IntGzip is fastq, fasta, or other and conver the file
mkdir -p FastaOnly
echo "Converting FastQ to Fasta"
parallel -j $ncores --bar "FastaorFastq {} ${out}FastaOnly" ::: IntGzip/*
rm -rf IntGzip

# Now to string deduplicate the files as I'd like to speed up the blast runs
echo "String Deduplciation with prinseq"
parallel -j $ncores --bar "perl /home/sam/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fasta {} -out_good ${out}StringDedup/{/.} -min_len $len -derep 14 -log prinseqLog/${sample}Merged.log" ::: ${out}FastaOnly/*

# The final step is to test if StringDeduplication will occur
if [[ $subsample != 0 ]]; then
	mkdir -p ${out}Subsample
	parallel -j $ncores --bar "seqkit sample {} -n $subsample > ${out}Subsample/{/}" ::: ${out}StringDedup/*
fi

#####################################
### Now to run the Blast commands ###
#####################################

if [[ $subsample != 0 ]]; then
	for final in StringDedup/*; do
		blastCMD ${out}Subsample/$final	
	done
else
	for final in StringDedup/*; do
		blastCMD ${out}StringDedup/$final	
	done

fi



