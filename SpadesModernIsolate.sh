#! /usr/bin/env bash
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.

SPAdesFunction(){
	spades --isolate -1 $r1 -2 $r2 --merged $merged -t $ncores -o ${out}/$sample
}
# These are the files and variables that will be needed
usage() { printf 'SPAdes Modern Assembly Script V0.5
	Absolute minimal script to get SPAdes running
        -i\tThe folder that contains the Raw sequencing files
	-o\tThe output prefix (Default: ModernIsolate)
	-n\tNumber of CPU Threads to be used (Default: 16)
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "SPAdes Isolate Settings (Modern) for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

out="ModernIsolate"
ncores=16
log="$(date +'%Y%m%d').log"

while getopts "i:n:o:l:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
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
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

if [ -z ${folder} ]; then # Testing if I have an input folder.  That's all I need to get running
	printf "The input was not set"
	usage
	exit 1
fi

# Writing the Log File
log | tee $log

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

mkdir -p ${out}
mkdir -p ${out}/SPAdesLogs
#########################
###Making the Assembly###
#########################

total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in ${samples[@]}; do
	FileIdentification $sample
	FileExtraction
	SPAdesFunction > ${out}/SPAdesLogs/$sample.log

	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done
