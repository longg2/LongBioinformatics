#! /usr/bin/env sh
######################################
### Functions that I'll be calling ###
######################################
usage() { printf 'BWA aln Mapping V1
	Mapping using aln is more complicated than it warrants,
	so here is a useful wrapper that will hopefully make life
	easier.
	-i\tThe folder the fasta files (REQUIRED)
	-o\tOutput folder of the BAM files (Default: BWAMappingScript)
	-r\tThe BWA index.  Will be mapping reads against it (REQUIRED)
	-k\tMinimum Fragment Length (Default: 30)
	-q\tMinimum Mapping Quality (Default: 30)
	-d\tDeduplicate the bam files? (Default: FALSE)
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
	-------------------------------------\n"; exit 0;
}

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

DeduplicateArray(){ # Deduplicating the sample names
	local array=("$@") # This feels weird, but, it'll allow you pass an array to it
	#echo "$array"
	local sampleNames=$(for name in ${array[@]}; do tmp=$(echo ${name/_*}); echo ${tmp/%.f*}; done) # Getting only the sample names
	samples=( $(echo ${sampleNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') )

}

FastaorFastq(){ # Figuring out if fasta or fastq
	local tmp=$1
	if file $tmp | grep -q "compressed"; then 
		local firstChar=$(zcat $1 | head -c 1)
	else 
		local firstChar=$(cat $1 | head -c 1)
	fi

	if [ "$firstChar" == "@" ];then
		exit 0
	elif [ "$firstChar" == ">" ];then
		exit 1
	else
		exit 2
	fi
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

alnMapping(){ # BWA aln Mapping.  Automatically determines if merged or paired.
	if [ "$merged" != "NA" ]; then # If I found a merged file
	#if [ -v merged ]; then # If I found a merged file
		bwa aln -o 2 -n 0.01 -l 16500 $ref $merged -t $ncores 2>/dev/null > tmp.sai
	
		bwa samse $ref tmp.sai $merged 2>/dev/null |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam|\
			samtools sort -	> tmpM.bam
	
		samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${sample}.fastq.gz
	fi
	
	if [ "$r1" != "NA" ]; then # If I found a merged file
	#if [ -v r1 ]; then # If we have paired reads
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r1 -t $ncores 2>/dev/null > tmp1.sai # First index
		bwa aln -o 2 -n 0.01 -l 16500 $ref $r2 -t $ncores 2>/dev/null > tmp2.sai # Second Index
	
		bwa sampe $ref tmp1.sai tmp2.sai $r1 $r2 2>/dev/null |\
			samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\
		       	samtools sort - > tmpP.bam
		
		samtools fastq -c 6 tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null
	
	fi
	
	# Merging the files if needed.  Otherwise, we can ignore and simply mv it
	#if [ -v merged ] && [ -v r1 ]; then
	if [ "$merged" != "NA" ] && [ "r1" != "NA" ]; then
		samtools merge -f ${out}MappedReads/$sample.bam tmpM.bam tmpP.bam
	elif [ "$merged" != "NA" ] && [ "r1" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		mv tmpM.bam ${out}MappedReads/$sample.bam
	else
		mv tmpP.bam ${out}MappedReads/$sample.bam
	fi

	# Deleting temporary files
	rm tmpP.bam tmpM.bam tmp.bam tmp*.sai

}

######################
### Default values ###
######################

declare -i qual=30
declare -i len=30
out="BWAMappingScript"
declare -i ncores=8
log="$(date +'%Y%m%d').log"
dedup="FALSE"

##############
### The UI ###
##############

while getopts "i:o:q:r:l:k:n:hd" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
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
	exit 1
fi

#################
### The Setup ###
#################

log | tee $log # The inital log file

# The Folders
mkdir -p ${out}MappedReads
mkdir -p ${out}UnmappedReads 
[ "${dedup}" == "TRUE" ] && mkdir -p ${out}DeduplicatedMappings

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

# The actual loop
total=${#samples[@]}
count="1"

for sample in ${samples[@]}; do # Iterating over an array of Samples
	FileIdentification $sample # Extracting the file names.  Will be saved as $sampleFiles
	FileExtraction

	#printf "$sample\n"
	#printf "MERGED:$merged\nR1:$r1\nR2:$r2\n" | tee -a $log # Debugging only

	# This here is to prevent odd scenarios where I only have r1 or Merged + r2
	if [ "$merged" != "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
	#if [ -v $merged ] && [ -v $r1 ] && [ -v $r2 ]; then
		#printf "$sample will run the whole shebang\n--------\n"
		alnMapping
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v $merged ] && [ -z ${r1+x} ] && [ -z ${r2+x} ]; then
		#printf "$sample will run only the merged file\n--------\n"
		alnMapping
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
	#elif [ -z ${merged+x} ] && [ -v $r1 ] && [ -v $r2 ]; then
		#printf "$sample will only run the paired file\n--------\n"
		alnMapping
	else
		printf "$sample has an odd combination.  It has been skipped\n" >> $log
	fi

	# A simple counter
	printf "$sample has been filtered ($count/$total)\n--------\n"
	count=$(echo "$count + 1" | bc)

done

# Only if requested
if [ "${dedup}" == "TRUE" ]; then
	echo "Deduplication Requested"
       parallel -j $ncores --bar "/usr/local/biohazard/bin/bam-rmdup -c -o ${out}DeduplicatedMappings/{/} {} 2> /dev/null" ::: ${out}MappedReads/*bam # Removes Duplicates
fi
