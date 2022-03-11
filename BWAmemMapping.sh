#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/MappingFunctions.sh # Contains mem and aln mapping code

usage() { printf 'BWA mem Mapping V1
	BWA MEM mapping.  Much simpler than with BWA aln, but, a wrapper
	to do the workup is still nice to have
	-i\tThe folder the fasta files (REQUIRED)
	-o\tOutput folder of the BAM files (Default: BWAMappingScript)
	-r\tThe BWA index.  Will be mapping reads against it (REQUIRED)
	-k\tMinimum Fragment Length (Default: 30)
	-q\tMinimum Mapping Quality (Default: 30)
	-d\tDeduplicate the bam files? (Default: FALSE)
	-n\tNumber of CPU Threads to be used (Default: 8)
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

log() {	printf "BWA mem settings for $(date):
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
mkdir -p ${out}BWALogs
[ "${dedup}" == "TRUE" ] && mkdir -p ${out}DeduplicatedMappings

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples
#echo ${samples[@]}

# The actual loop
total=${#samples[@]}
count=0

ProgressBar $count $total
for sample in ${samples[@]}; do # Iterating over an array of Samples

	FileIdentification $sample # Extracting the file names.  Will be saved as $sampleFiles
	FileExtraction

	#printf "$sample\n"
	#printf "\nMERGED:$merged\nR1:$r1\nR2:$r2\n" | tee -a $log # Debugging only
	memMapping 2> ${out}BWALogs/$sample.log

	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total

done

mkdir ${out}Depths/
# Now we split based on if we want to deduplicate
if [ "${dedup}" == "TRUE" ]; then
	echo "Deduplication Requested"
	parallel -j $ncores --bar "/usr/local/biohazard/bin/bam-rmdup -c -o ${out}DeduplicatedMappings/{/} {} > /dev/null 2> /dev/null" ::: ${out}MappedReads/*bam # Removes Duplicates
	parallel -j $ncores --bar "samtools depth -aa {} > ${out}Depths/{/.}.tab" ::: ${out}DeduplicatedMappings/*bam # Getting the read depths. Something that I end up doing often anyways
	#parallel -j $ncores --bar "$script_full_path/DepthStatistics.awk {}" ::: ${out}Depths/*bam > DepthStatistics.tab # This is the script that will calculate the depths.
	gzip ${out}Depths/*
else
	parallel -j $ncores --bar "samtools depth -aa {} > ${out}Depths/{/.}.tab" ::: ${out}MappedReads/*bam # Getting the read depths. Something that I end up doing often anyways
#	parallel -j $ncores --bar "$script_full_path/DepthStatistics.awk {}" ::: ${out}Depths/*bam > DepthStatistics.tab # This is the script that will calculate the depths.
	gzip ${out}Depths/*
fi

printf "\nMapping Complete\n"

