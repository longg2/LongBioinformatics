#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
usage() { printf 'BWA aln Mapping V1.1
	Mapping using aln is more complicated than it warrants,
	so here is a useful wrapper that will hopefully make life
	easier.
	V1.1 -> Multimappers Separated

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
	#local sampleNames=$(for name in ${array[@]}; do tmp=$(echo ${name/%%_*}); echo ${tmp/%.f*}; done) # Getting only the sample names
	local sampleNames=$(for name in ${array[@]}; do
		tmp="$(echo $name |sed -e 's/_r1.*//I' -e 's/_r2.*//I' -e "s/_merged.*//I" -e "s/\.f.*//I")";
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

	if [ "$r1" != "NA" ] && ["$r2" != "NA" ]; then # If I found a Paired set
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
		mv tmpMBad.bam ${out}UnmappedReads/$sample.bam 
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v merged ] && [ -z ${r1+x} ]; then
		mv tmpS.bam  ${out}MappedReads/$sample.bam 
		mv tmpSBad.bam ${out}UnmappedReads/$sample.bam 
	else
		#echo "PAIRED"
		mv tmpP.bam ${out}MappedReads/$sample.bam
		mv tmpPBad.bam ${out}UnmappedReads/$sample.bam 
	fi

	# Deleting temporary files
	rm tmp.bam tmp*.sai

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

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples
#echo ${samples[@]}

# The actual loop
#echo "${samples[@]}"
total=${#samples[@]}
count=0

ProgressBar $count $total
for sample in ${samples[@]}; do # Iterating over an array of Samples
	FileIdentification $sample # Extracting the file names.  Will be saved as $sampleFiles
	FileExtraction

#	printf "\n$sample"
#	printf "\nMERGED:$merged\nR1:$r1\nR2:$r2\n" #| tee -a $log # Debugging only

	# This here is to prevent odd scenarios where I only have r1 or Merged + r2
	if [ "$merged" != "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
	#if [ -v $merged ] && [ -v $r1 ] && [ -v $r2 ]; then
		#printf "$sample will run the whole shebang\n--------\n"
		alnMapping 2> ${out}BWALogs/$sample.log
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v $merged ] && [ -z ${r1+x} ] && [ -z ${r2+x} ]; then
		#printf "$sample will run only the merged file\n--------\n"
		alnMapping 2> ${out}BWALogs/$sample.log
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
	#elif [ -z ${merged+x} ] && [ -v $r1 ] && [ -v $r2 ]; then
		#printf "$sample will only run the paired file\n--------\n"
		alnMapping 2> ${out}BWALogs/$sample.log
	else
		printf "$sample has an odd combination.  It has been skipped\n" 
	fi

	# A simple counter
	#printf "$sample has been filtered ($count/$total)\n--------\n"
	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done

printf "\nMapping Complete\n"

# Only if requested
if [ "${dedup}" == "TRUE" ]; then
	echo "Deduplication Requested"
       parallel -j $ncores --bar "/usr/local/biohazard/bin/bam-rmdup -c -o ${out}DeduplicatedMappings/{/} {} > /dev/null 2> /dev/null" ::: ${out}MappedReads/*bam # Removes Duplicates
fi
