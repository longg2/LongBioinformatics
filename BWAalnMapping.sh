#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)
export script_full_path

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/MappingFunctions.sh # The mapping functions
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
mkdir -p ${out}BWALogs
mkdir -p ${out}Depths
[ "${multi}" == "TRUE" ] && mkdir -p ${out}MQ0Bam 
#[ "${multi}" == "TRUE" ] && mkdir -p ${out}Qual0Reads 
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
	elif [ "$merged" == "NA" ] && [ "$r1" != "NA" ] && [ "$r2" == "NA" ]; then
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

# Now we split based on if we want to deduplicate

# This is taking into account cases where the slash isn't being given
if [ "$out" != "" ]; then
	outFolder="${out}/"
else
	outFolder=$out
fi

if [ "${dedup}" == "TRUE" ]; then
	printf "\nDeduplication Requested\n"
	parallel -j $ncores --bar "/usr/local/biohazard/bin/bam-rmdup -c -o ${out}DeduplicatedMappings/{/} {} > /dev/null 2> /dev/null" ::: ${out}MappedReads/*bam # Removes Duplicates
	printf "Getting Read Depths\n"
	parallel -j $ncores --bar "samtools depth -a {} > ${out}Depths/{/.}.tab" ::: ${out}DeduplicatedMappings/*bam # Getting the read depths. Something that I end up doing often anyways
	printf "Getting Mean and SD of Depths\n"
	parallel -j $ncores --bar "$script_full_path/DepthStatistics.awk {}" ::: ${out}Depths/*tab > ${outFolder}DepthStatistics.tab # This is the script that will calculate the depths.
	gzip ${out}Depths/*
else
	printf "\nGetting Read Depths\n"
	parallel -j $ncores --bar "samtools depth -a {} > ${out}Depths/{/.}.tab" ::: ${out}MappedReads/*bam # Getting the read depths. Something that I end up doing often anyways
	printf "Getting Mean and SD of Depths\n"
	parallel -j $ncores --bar "$script_full_path/DepthStatistics.awk {}" ::: ${out}Depths/*tab > ${outFolder}DepthStatistics.tab # This is the script that will calculate the depths.
	gzip ${out}Depths/*
fi

printf "Mapping Complete\n"

