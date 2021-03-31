#! /usr/bin/env sh
# These are the functions that actually do the work.  Mean to make parallelization easier to do
FileIdentification(){ # Extract Files from an array using results from another array
	local sample=$1 # Basename of the file
	local arrayFiles=$files

	# Finding the indices which have the same samplename
	sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))
#	local numFiles=$(echo "$sampleFiles" | wc -l )
#	#echo "$numFiles"
#
#	if [ $numFiles = 1 ]; then
#		echo "$sample only has one file.  Assuming all Merged"
#		sampleFiles=($( echo $sampleFiles | tr '\012' ' ' ))
#		#FastaorFastq "$folder/$sampleFiles"	
#		## What file is it?
#		#if [ $? == "0" ]; then
#		#	echo "It's a fastq file.  Proceed to mapping."
#		#elif [ $? == "1" ]; then
#		#	echo "It's a fasta file.  Won't be mapping it."
#		#else
#		#	echo "Unknown file.  Halting."
#		#	exit 1
#		#fi
#
#		# Need to test if it's fasta or fastq
#	elif [ $numFiles = 2 ]; then
#		echo "$sample has two files.  Assuming only Paired"
#		sampleFiles=($( echo $sampleFiles | tr '\012' ' ' ))
#
#		#FastaorFastq "$folder/${sampleFiles[0]}"	
#		## What file is it?
#		#if [ $? == "0" ]; then
#		#	echo "It's a fastq file.  Proceed to mapping."
#		#elif [ $? == "1" ]; then
#		#	echo "It's a fasta file.  Won't be mapping it."
#		#else
#		#	echo "Unknown file.  Halting."
#		#	exit 1
#		#fi
#	elif [ $numFiles = 3 ]; then
#		echo "$sample has three files.  Assuming both paired and merged reads"
#		sampleFiles=($( echo $sampleFiles | tr '\012' ' ' ))
#
#		#echo "$folder/${sampleFiles[0]}"
#		#FastaorFastq "$folder/${sampleFiles[0]}"
#		## What file is it?
#		#if [ $? == "0" ]; then
#		#	echo "It's a fastq file.  Proceed to mapping."
#		#elif [ $? == "1" ]; then
#		#	echo "It's a fasta file.  Won't be mapping it."
#		#else
#		#	echo "Unknown file.  Halting."
#		#	exit 1
#		#fi
#	else
#		echo "$sample has more than three files.  Will be skipped."
#	fi

}

FileExtraction(){ # Assign files to their variables.  Assumes that $sampleFiles and $sample exists
	# Unsetting variables in case they're already defined from a previous run 
	# unset merged
	# unset r1
	# unset r2

	# This here is a fix that's needed if -v doesn't
	merged="NA"
	r1="NA"
	r2="NA"

	# Creating a hidden text file of file names
	printf '%s\n' "${sampleFiles[@]}" > .hiddenlist.list

	# Identifying the files
	if grep -P -i -q "r1" .hiddenlist.list; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r1')
	        r1="$folder/$fileName"
	fi

	if grep -P -i -q "r2" .hiddenlist.list; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r2')
	        r2="$folder/$fileName"
	fi

	# Two cases for the merged.  Want to control for shenanigans
	if grep -q -i "$sample\.f.*" .hiddenlist.list; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -v -q "_\.f*") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -i "$sample\.f.*")
	        merged="$folder/$fileName"
	fi

	if grep -P -i -q "merged" .hiddenlist.list; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'merged')
	        merged="$folder/$fileName"
	fi

	rm .hiddenlist.list
}

Trimming(){ # Performing the trimming
	# Trimming
	local r1=$1
	local r2=$1
	local sample=$3
	local out=$4
        fastp -i $r1 -I $r2 --merge \
        --merged_out ${out}Trimmed/${sample}_merged.fastq.gz \
	--out1 ${out}Trimmed/${sample}_r1.fastq.gz \
	--out2 ${out}Trimmed/${sample}_r2.fastq.gz \
        --adapter_fasta /usr/local/trimmomatic/adapters/TruSeq3-PE-2.fa --correction \
        --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
        --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
        --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
        --overlap_diff_limit 100 --n_base_limit 0\
        --overlap_len_require 15 --length_required $len\
        --html ${out}FastpLogs/${sample}.html \
        --json ${out}FastpLogs/${sample}.json -R $sample --thread 16 -R $sample \
        --unpaired1 ${out}Trimmed/${sample}_u1.fastq.gz \
        --unpaired2 ${out}Trimmed/${sample}_u2.fastq.gz\
        --failed_out ${out}FailedQC/${sample}_failed.fastq.gz;
}

FastpWrapper(){ # Convenient Wrapper for parallelization
	FileIdentificationInFunction $1 $folder
	FileExtraction $folder
	Trimming
}

#StringDeduplication(){ # Convenient Wrapper for parallelization
#	FileIdentification $1
#	FileExtraction $2
#	
#	# Now for the actual deduplication
#
#	if [ "$merged" != "NA" ]; then # If I found a merged file
#		zcat $merged | perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq stdin -out_bad null -out_good stdout -min_len $len -derep 14 -log ${out}PrinseqLog/${sample}Merged.log 2> /dev/null | gzip > ${out}PooledLanes/${sample}_Merged.fastq.gz
#	fi
#
#	if [ "$r1" != "NA" ]; then # If I found a paired file
#		perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq <(zcat $r1) -fastq2 <(zcat $r2) -out_bad null -out_good TMP/${sample}_r -min_len $len -derep 14 -log ${out}PrinseqLog/${sample}Paired.log 2> /dev/null
#
#	# Because of how prinseq is coded, I'll need to compress separately	
#		gzip -c ${sample}_r1.fastq > ${out}PooledLanes/${sample}_r1.fastq.gz
#		gzip -c ${sample}_r2.fastq > ${out}PooledLanes/${sample}_r2.fastq.gz
#
#		rm -rf TMP/*
#	fi
#}

#export -f FileIdentificationInFunction
#export -f FileExtraction
export -f Trimming
#export -f FastpWrapper


DeduplicateArray(){ # Deduplicating the sample names
	local array=("$@") # This feels weird, but, it'll allow you pass an array to it
	local sampleNames=$(for name in ${array[@]}; do
		tmp="$(echo $name |sed -e 's/_r1.*//I' -e 's/_r2.*//I' -e 's/_merged.*//I' -e 's/\.fa.*//I')";
		echo $tmp;
       	done) # Getting only the sample names
	samples=( $(echo ${sampleNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') )

}

# These are the files and variables that will be needed
usage() { printf 'Modern QC Script V1
	This is a minimal script.  All it will do is trim, merge,
	and Pool the reads.  Assumes you have fastp.
        -i\tThe folder that contains the Raw sequencing files
	-o\tThe output prefix (Default: QC)
	-n\tNumber of CPU Threads to be used (Default: 16)
	-l\tLog File Name (Default: $date)
	-k\tMinimum Read Length (Default: 30)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Mapping settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

#Default Values
out="QC"
ncores=8
export len=30
log="$(date +'%Y%m%d').log"

while getopts "i:k:n:o:l:hs" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
			export "${folder}"
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                o)
                        out=${OPTARG}
			export $out
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        #echo "Using $ncores"
                        ;;
                l)
                        log=${OPTARG}
                        #echo "Settings are being outputted to $log"
                        ;;
		s)
			dedupMethod="String"
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done
logBWA="${log/%.*}BWA.log"

# Writing the Log File
log | tee $log

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples
#echo ${samples[@]}
#exit 0
####################################
###Trimming the reads with leeHom###
####################################
mkdir -p ${out}Trimmed
mkdir -p ${out}FastpLogs
mkdir -p ${out}FailedQC
echo "Trimming and Merging Reads"

njobs=$(echo "scale=0;var1=$ncores/16;var1"|bc) # Will round down!!!

for sample in ${samples[@]}; do
	FileIdentification $sample
	FileExtraction
	sem -j $njobs "Trimming $r1 $r2 $sample $out" > /dev/null 2> /dev/null
	
done

#parallel --bar -j $njobs "FastPWrapper {}" ::: "${samples[@]}"

###################
###Pooling Lanes###
###################

mkdir -p ${out}PooledLanes
echo "Pooling Lanes"

parallel --bar -j $ncores "cat ${out}Trimmed/{}_*merged.fastq.gz > ${out}PooledLanes/{}.fastq.gz;
	cat ${out}Trimmed/{}_*r1.fastq.gz > ${out}PooledLanes/{}_r1.fastq.gz;
	cat ${out}Trimmed/{}_*r2.fastq.gz > ${out}PooledLanes/{}_r2.fastq.gz;" ::: "${samples[@]}"

####################################
### Running String Deduplication ###
####################################
#if [ $dedupMethod = "String" ]; then
#
#	# The setup
#	echo "Performing String Deduplication"
#	mkdir -p ${out}PrinseqLog
#	mv ${out}PooledLanes ${out}PooledLanesDups # Want to keep the non deduplicated section in case there's stuff of interest there
#	mkdir -p ${out}PooledLanes
#	mkdir -p TMP
#
#	parallel --bar -j $ncores "StringDeduplication {} ${out}PooledLanesDups" ::: "${samples[@]}"
#
#	rm -rf TMP
#
#fi
