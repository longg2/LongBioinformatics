#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
usage() { printf 'BWA aln Mapping Split Test

	-i\tThe folder the fasta files (REQUIRED)
	-o\tOutput folder of the BAM files (Default: BWAMappingScript)
	-r\tThe BWA index.  Will be mapping reads against it (REQUIRED)
	-k\tMinimum Fragment Length (Default: 30)
	-q\tMinimum Mapping Quality (Default: 30)
	-d\tDeduplicate the bam files? (Default: FALSE)
	-m\tSeparate the MultiMappers? (Default: FALSE)
	-n\tNumber of CPU Threads to be used (Default: 8)
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

log() {	printf "BWA aln settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Reference:\t${ref}
	CPU Threads:\t${ncores}
	-------------------------------------
	Min Quality:\t${qual}
	Min Length:\t${len}
	Deduplicated:\t${dedup}
	Multimappers Separated:\t${multi}
	-------------------------------------\n"; exit 0;
}

FileIdentification(){ # Extract Files from an array using results from another array
	local sample=$1 # Basename of the file
	local arrayFiles=$files

	# Finding the indices which have the same samplename
	sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))

}
DeduplicateArray(){ # Deduplicating the sample names
	local array=("$@") # This feels weird, but, it'll allow you pass an array to it
	#echo "$array"
	#local sampleNames=$(for name in ${array[@]}; do tmp=$(echo ${name/%%_*}); echo ${tmp/%.f*}; done) # Getting only the sample names
	local sampleNames=$(for name in ${array[@]}; do
		tmp="$(echo $name |sed -e 's/_r1.*//I' -e 's/_r2.*//I' -e 's/_merged.*//I' -e 's/\.f.*//I')";
		echo $tmp;
       	done) # Getting only the sample names
	#local sampleNames=$(for tmp in ${array[@]}; do tmp=$(echo ${tmp/_r1}); tmp=$(echo ${tmp/_r2}); tmp=$(echo ${tmp/_merged});echo ${tmp/%.f*}; done) # Getting only the sample tmps
	samples=( $(echo ${sampleNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') )

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
	if grep -P -i -q "_r1" .hiddenlist.list; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i '_r1')
	        r1="$folder/$fileName"
	fi

	if grep -P -i -q "_r2" .hiddenlist.list; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i '_r2')
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

SplitAlnMapping(){ # BWA aln Mapping.  Split ALN mapping is a different beast. Will be run in parallel, but, that requires a rewrite in certain parts.  Only going
		   # to do deal with SE for now.
	local merged=$1
	local ncores=$2
	local jobNum=$3

	bwa aln -o 2 -n 0.01 -l 16500 $ref $merged -t $ncores > tmp_$jobNumb.sai

	bwa samse $ref tmp.sai $merged |\
		samtools view -b -h -F 4 -m $len -q $qual -U tmpMBad_$jobNumb.bam |\
		samtools sort -	> TMPBAM/tmpM_$jobNumb.bam  

	# Now to figure out if Multimappers are to be extracted
	if [ "$multi" != "TRUE" ]; then
		samtools fastq tmpMBad_$jobNumb.bam | gzip > TMPUNMAPPED/${sample}_$jobNumb.fastq.gz

	else
		# Apr 20 2021 -- Want to separate the poor hits from the reads with multiple hits
		samtools fastq -f 4 tmpMBad_$jobNumb.bam | gzip > tmpUnmapped_$jobNumb.fastq.gz #${out}UnmappedReads/${sample}.fastq.gz # All the unmapped Reads
		samtools view -F 4 tmpMBad_$jobNumb.bam | awk '{if($5 > 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip -c | cat tmpUnmapped_$jobNumb.fastq.gz -  > TMPUNMAPPED/${sample}_$jobNumb.fastq.gz # The reads which were poor matches 
		samtools view -m $len -F 4 tmpMBad.bam | awk '{if($5 == 0){print "@"$1"\n"$10"\n+\n"$11}}' | gzip > TMPQUAL0/${sample}_$jobNumb.fastq.gz # The Multimappers
	fi
	
}

ProgressBar() { # From github.com/fearside/ProgressBar
	# Process data
		let _progress=(${1}*100/${2}*100)/100
		let _done=(${_progress}*4)/10
		let _left=40-$_done
	# Build progressbar string lengths
	_done=$(printf "%${_done}s")
	_left=$(printf "%${_left}s")
	printf "\rProgress : [${_done// />}${_left// /-}] ${_progress}%%"

}	

export -f SplitAlnMapping

######################
### Default values ###
######################

declare -i qual=30
declare -i len=30
out="BWAMappingScript"
declare -i ncores=8
log="$(date +'%Y%m%d').log"
dedup="FALSE"
multi="FALSE"

##############
### The UI ###
##############

while getopts "i:o:q:r:l:k:n:hdm" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find -L $folder/* -type f -printf "%f\n") # Making an array of files
                        #echo "The raw sequencing files are located in $in"
                        ;;
                o)
                        out=${OPTARG}
                        #echo "$out is the output folder"
                        ;;
		q)
			declare -i qual=${OPTARG}
			#echo "Requiring a minimum mapping quality of $qual"
                        ;;
                r)
                        declare -r ref=${OPTARG}
                        #echo "We're mapping the QCed reads against $ref"
                        ;;
                l)
                        log=${OPTARG}
                        #echo "Settings are being outputted to $log"
                        ;;
                k)
                        declare -i len=${OPTARG}
                        #echo "Using a minimum length of ${len}bp"
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        #echo "Using $ncores CPU Threads"
                        ;;
		m)
			multi="TRUE"
			;;
		d)
			dedup="TRUE"
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Testing if files options are missing
if [ -z ${folder+x} ] || [ -z ${ref+x} ]; then
	echo "You are missing either the Input folder or the reference"
	usage
	exit 1
fi

#################
### The Setup ###
#################

log | tee $log # The inital log file

# The Folders
mkdir -p ${out}MappedReads
mkdir -p ${out}UnmappedReads 
mkdir -p ${out}UnmappedBam 
[ "${multi}" == "TRUE" ] && mkdir -p ${out}Qual0Reads 
mkdir -p ${out}BWALogs
[ "${dedup}" == "TRUE" ] && mkdir -p ${out}DeduplicatedMappings

# Exporting needed functions
export ref

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples
#echo ${samples[@]}

# The actual loop
#echo "${samples[@]}"
total=${#samples[@]}
count=0

ProgressBar $count $total
printf "\n"
for sample in ${samples[@]}; do # Iterating over an array of Samples
	FileIdentification $sample # Extracting the file names.  Will be saved as $sampleFiles
	FileExtraction

#	printf "\n$sample"
#	printf "\nMERGED:$merged\nR1:$r1\nR2:$r2\n" #| tee -a $log # Debugging only
	
	# Now to split the file into 1M reads
	#seqkit split2 -s 1000000 -j $ncores -O SplitFiles -1 $merged  #2> /dev/null # Only informing that the split is occuring
	export parThreads=$(( $ncores/3 ))
	export sample

	# Making Temporary folders
	mkdir -p TMPBAM
	mkdir -p TMPUNMAPPED
	[ "${multi}" == "TRUE" ] && mkdir -p TMPQUAL0

	# Now to do the mapping
	parallel -j 3 --bar "SplitAlnMapping {} $parThreads {#} 2> ${out}BWALogs/{/.}.log" ::: SplitFiles/*

	#Merging the files together
	samtools merge ${out}MappedReads/$sample.bam TMPBAM/*
	cat TMPUNMAPPED/*.gz > ${out}UnmappedReads/$sample.fastq.gz
	[ "${multi}" == "TRUE" ] && cat TMPQUAL0/*.gz ${out}Qual0Reads/$sample.fastq.gz

	# A simple counter
	#rm -rf SplitFiles TMPBAM TMPUNMAPPED
	[ "${multi}" == "TRUE" ] && rm -rf TMPQUAL0
	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
	printf "\n"
done

printf "\nMapping Complete\n"

# Only if requested
if [ "${dedup}" == "TRUE" ]; then
	echo "Deduplication Requested"
       parallel -j $ncores --bar "/usr/local/biohazard/bin/bam-rmdup -c -o ${out}DeduplicatedMappings/{/} {} > /dev/null 2> /dev/null" ::: ${out}MappedReads/*bam # Removes Duplicates
fi
