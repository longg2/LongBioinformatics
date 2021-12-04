#! /usr/bin/env bash

# Getting the commands I neeed
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/QCFunctions.sh # This loads the QC commands
usage() { printf 'Ancient QC Script V1
	Very simple.  Want this to focus primarily the workup
	prior to analysis.
        -i\tThe folder that contains the Raw sequencing files
	-o\tThe Output Prefix
        -n\tNumber of CPU Threads to be used
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Ancient Mapping settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	CPU Threads:\t${ncores}\n"; exit 0;
}

##################################
#Default Values
out="AncientQC/"
declare -i ncores=8
log="$(date +'%Y%m%d').log"

while getopts "i:o:n:l:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
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
                o)
                        out="${OPTARG}"
                        #echo "Settings are being outputted to $log"
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done
###################################################################
# Writing the Log File
log | tee $log

mkdir -p ${out}Trimmed
mkdir -p ${out}Failed
mkdir -p ${out}PooledLanes
###################################################################

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

echo "Trimming the ${#samples[@]} samples found in ${folder}.  This may take a while"
total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in ${samples[@]}; do
	FileIdentification $sample
	FileExtraction
	AncientTrimming
	
	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done

###################
###Pooling Lanes###
###################

pooledNames=$(for name in ${samples[@]}; do # This was lifted from DeduplicateArray().
	tmp="$(echo $name |sed -e 's/_L00.*//I')";
		echo $tmp;
       	done) # Getting only the sample names
pooledNames=( $(echo ${pooledNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') ) # This was lifted from Deduplicate Array

mv ${out}Trimmed/*fail* ${out}Failed # Hiding the Failed before anything happens....
printf "\nPooling the lanes together\n"
parallel --bar -j $ncores "cat ${out}Trimmed/{}_L00{1,2}.fq.gz > ${out}PooledLanes/{}.fastq.gz 2> /dev/null;
	cat ${out}Trimmed/{}*r1.fq.gz > ${out}PooledLanes/{}_r1.fastq.gz;
	cat ${out}Trimmed/{}*r2.fq.gz > ${out}PooledLanes/{}_r2.fastq.gz;" ::: "${pooledNames[@]}" # The 2> is included as we won't have two lanes all the time.  Want to hide these errors if that is the case since it's working as intended.  Need a better way to do it....
