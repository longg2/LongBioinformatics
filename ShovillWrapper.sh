#! /usr/bin/env bash
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I use in every script
source $script_full_path/lib/AssemblyFunctions.sh # This is needed for SPAdes
#source $script_full_path/lib/QCFunctions.sh # This loads the all of the QC functions

#export -f AdapterRemovalHeavy
#export -f GzipDetection
export -f FileIdentificationInFunction
export -f FileExtractionInFunction

# These are the files and variables that will be needed
usage() { printf "Shovill Wrapper V0.1
	Should be seen as a better spades assembly script.
	Includes assembly correction with Pilon by default, however,
	we want to potentially run successive refinements
        -i\tThe folder that contains the Raw sequencing files
	-o\tThe output prefix (Default: Assembly)
	-n\tNumber of CPU Threads to be used (Default: 16)
	-p\tNumber of additional pilon runs
	-P\tOnly run the pilon steps. Assumes this script was run before!
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit
	-------------------\n" 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Denovo Assembly settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Additional Pilon:\t${pilonRuns}
	Pilon Only?:\t${pilonOnly}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

out="Assembly"
ncores=16
pilonRuns=0
pilonOnly="FALSE"
log="$(date +'%Y%m%d').log"


while getopts "i:n:o:l:Pp:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find -L $folder/* -type f -printf "%f\n") # Making an array of files
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                o)
                        out=${OPTARG}
                        #echo "The raw sequencing files are located in $raw"
                        ;;
		p)
			declare -i pilonRuns=${OPTARG}
			;;
		P)
			pilonOnly="TRUE"
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

# Exporting the required files
export folder
export out

# Need to make the ouput folder
mkdir -p ${out}
# Writing the Log File
log | tee $log

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

###############################################
###Let's filter by length and trim the reads###
###############################################
mkdir -p ${out}/Logs

njobs=$(echo "scale=0;var1=$ncores/8;var1"|bc) # Will round down!!!

if [[ $njobs -eq 0 ]]; then
	njobs=1
	#export njobs # This is to tell the fastp function how many threads we're working with
fi

if [[ $pilonOnly != "TRUE" ]]; then
	export -f shovillFunction
	printf "\nCreating the assembly\n"
	parallel -j $njobs --bar "shovillFunction {} 2> /dev/null" ::: ${samples[@]}
fi
# If pilon runs were called!
if [[ $pilonRuns != 0 ]]; then
	export -f pilonFunction
	printf "\n$pilonRuns Pilon run(s) requested\n"
	for ((i = 1; i <= $pilonRuns; i++)); do # This is the for loop for successive runs
		printf "Iteration $i\n"
		parallel -j $njobs --bar "pilonFunction {} $i" ::: ${samples[@]}
	done

fi
printf "\n"
