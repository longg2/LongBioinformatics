#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh || { echo 'Missing BasicCommands.sh'; exit 1;} # This loads the basic things I need.
source $script_full_path/lib/MappingFunctions.sh  || { echo 'Missing MappingFunctions.sh'; exit 1;} # The mapping functions
source $script_full_path/lib/QCFunctions.sh || { echo 'Missing QCFunctions.sh'; exit 1;} # The mapping functions

usage() { printf 'vg Mapping V0.1
	Implementing the use of variation graphs for new aDNA projects.
	See https://doi.org/10.1186/s13059-020-02160-7 and
	https://doi.org/10.1186/s13059-020-1941-7 for rationale.

	NOTES: 
		• Only works on info2020 due to old kernels on 113-115.
		• Cannot run bam-rmdup on info2020.  Need to talk to Brian...

	-i\tThe folder the fasta files <--- REQUIRED
	-o\tOutput folder (Default: BWAMappingScript)
	-r\tThe basename of the vg index.  <--- REQUIRED
	-k\tMinimum Fragment Length (Default: 30)
	-q\tMinimum Mapping Quality (Default: 50)
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

if [ ${HOSTNAME} == "info2020" ]; then
	echo "This script needs to be run on info2020"
	exit 1
elif ! command -v vg &> /dev/null ; then
	printf "To run this script, the vg program is needed. Please create a conda environment or talk to Sam.\n-------------\n"
	usage
	exit 1
elif [ -z ${folder+x} ] || [ -z ${ref+x} ]; then
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

#########################
### Mapping The Reads ###
#########################

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
		vgMapping 2> ${out}vgLogs/$sample.log
	elif [ "$merged" != "NA" ] && [ "$r1" == "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -v $merged ] && [ -z ${r1+x} ] && [ -z ${r2+x} ]; then
		#printf "$sample will run only the merged file\n--------\n"
		vgMapping 2> ${out}vgLogs/$sample.log
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" != "NA" ]; then
	#elif [ -z ${merged+x} ] && [ -v $r1 ] && [ -v $r2 ]; then
		#printf "$sample will only run the paired file\n--------\n"
		vgMapping 2> ${out}vgLogs/$sample.log
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
	#elif [ -z ${merged+x} ] && [ -v $r1 ] && [ -v $r2 ]; then
		#printf "$sample will only run the paired file\n--------\n"
		vgMapping 2> ${out}vgLogs/$sample.log
	else
		printf "$sample has an odd combination.  It has been skipped\n" 
	fi

	# A simple counter
	#printf "$sample has been filtered ($count/$total)\n--------\n"
	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done

printf "\nMapping Complete\n"

################################################
### Now, what if I want to get the variants? ###
################################################

# Commands:
# vg pack
# vg call
# How do I combine gam files -> gfakluge is a possibility
