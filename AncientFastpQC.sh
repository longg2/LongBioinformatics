#! /usr/bin/env bash
# These are the functions that actually do the work.
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/QCFunctions.sh # This loads the QC commands
# These are the files and variables that will be needed
usage() { printf 'Ancient QC Script V2
	This is minimal script. Assumes that all you want to do
	is trim, merge, and pool lanes together. Will detect adapters
	automatically. Requires fastp to run.
        -i\tThe folder that contains the Raw sequencing files
	-o\tThe output prefix (Default: QC)
	-n\tNumber of CPU Threads to be used (Default: 16)
	-l\tLog File Name (Default: $date)
	-k\tMinimum Read Length (Default: 30)
	-a\tDisable Adapter Detection (Default: FALSE)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Mapping settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

# Exporting the required functions
export -f Trimming
export -f FileIdentificationInFunction
export -f FileExtractionInFunction

#Default Values
out="ModernQC"
ncores=8
export len=30
log="$(date +'%Y%m%d').log"
#dedupMethod="NONE"
adapter="TRUE"

while getopts "i:k:n:o:l:ha" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find -L $folder/* -type f -printf "%f\n") # Making an array of files
			export folder
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                o)
                        out=${OPTARG}
			export out
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        #echo "Using $ncores"
                        ;;
                k)
                        len=${OPTARG}
			export len
                        #echo "Settings are being outputted to $log"
                        ;;
                l)
                        log=${OPTARG}
                        #echo "Settings are being outputted to $log"
                        ;;
		a)
			adapter="FALSE"
			;;
#		d)
#			dedupMethod="String"
#			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Writing the Log File
log | tee $log

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples
#echo ${samples[@]}
#exit 0
###################################
###Trimming the reads with Fastp###
###################################
mkdir -p ${out}Trimmed
mkdir -p ${out}FastpLogs
mkdir -p ${out}FailedQC
mkdir -p ${out}FastpLogNorm
echo "Trimming and Merging Reads"

njobs=$(echo "scale=0;var1=$ncores/16;var1"|bc) # Will round down!!!

if [[ $njobs -eq 0 ]]; then
	njobs=1
	export njobs # This is to tell the fastp function how many threads we're working with
fi

if [ "$adapter" == "TRUE" ]; then
	export -f FastpWrapperAncient # Only want to export the functions I need
	export -f AncientTrimmingFastp
	parallel -j $njobs --bar "FastpWrapperAncient {}" ::: ${samples[@]}
else
	export -f FastpWrapperAncientNoAdapter # Only want to export the functions I need
	export -f AncientTrimmingFastpNoAdapter
	parallel -j $njobs --bar "FastpWrapperAncientNoAdapter {}" ::: ${samples[@]}
fi

###################
###Pooling Lanes###
###################

mkdir -p ${out}PooledLanes
pooledNames=$(for name in ${samples[@]}; do # This was lifted from DeduplicateArray().
	tmp="$(echo $name |sed -e 's/_L00.*//I')";
		echo $tmp;
       	done) # Getting only the sample names
pooledNames=( $(echo ${pooledNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') ) # This was lifted from Deduplicate Array
#echo $pooledNames

printf "\nPooling the lanes together\n"

parallel --bar -j $ncores "cat ${out}Trimmed/{}*merged.fastq.gz > ${out}PooledLanes/{}.fastq.gz;
	cat ${out}Trimmed/{}*r1.fastq.gz > ${out}PooledLanes/{}_r1.fastq.gz;
	cat ${out}Trimmed/{}*r2.fastq.gz > ${out}PooledLanes/{}_r2.fastq.gz;" ::: "${pooledNames[@]}" # "${samples[@]}"

find ${out}PooledLanes -type f -empty -exec rm {} \;
###############################################
### Now to get the FLDs of the Merged Reads ###
###############################################

mkdir -p ${out}FLD

parallel --bar -j $ncores "zcat ${out}PooledLanes/{}.fastq.gz | $script_full_path/lib/FLDFastq.awk | sort -n | uniq -c | sed -e 's/^ *//g' -e 's/ /\t/g' > ${out}FLD/{}.tab" ::: "${pooledNames[@]}"
