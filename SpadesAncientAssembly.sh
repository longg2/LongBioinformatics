#! /usr/bin/env bash
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I use in every script
source $script_full_path/lib/AssemblyFunctions.sh # This is needed for SPAdes
source $script_full_path/lib/QCFunctions.sh # This loads the all of the QC functions

export -f AdapterRemovalHeavy
export -f GzipDetection
export -f FileIdentificationInFunction
export -f FileExtractionInFunction
# These are the files and variables that will be needed
usage() { printf 'SPAdes Ancient Assembly Script V0.2
	Implements suggestions from Ana to improve de novo success rate.
	Assumes basic QC has been done (ie. leeHom & Pooling)
	No deduplication of the reads!!
        -i\tThe folder that contains the Raw sequencing files
	-o\tThe output prefix (Default: ML50CHOMP0)
	-n\tNumber of CPU Threads to be used (Default: 16)
	-k\tMin Length (Default: 50)
	-c\tHow many bp to trim off ends? (Default: 0)
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit
	-------------------	
	TODO:
		• Figure out how to automate the addition of orphan reads to the assembly\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Denovo Assembly settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	CPU Threads:\t${ncores}
	-------------------------------------
	Minimum Length:\t${len}
	Chomp size:\t${chomp}
	-------------------------------------\n"; exit 0;
}

out="ML50CHOMP0"
ncores=16
len=50
chomp=0
log="$(date +'%Y%m%d').log"

while getopts "i:n:o:l:c:k:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                o)
                        out=${OPTARG}
                        #echo "The raw sequencing files are located in $raw"
                        ;;
		c)
			chomp=${OPTARG}
			#export chomp
			;;
		k)
			len=${OPTARG}
			#export len
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

# Need to decompress the files as prinseq can't handle them *shakes fist*
echo "Decompressing the files"
mkdir -p IntGzip
parallel -j $ncores --bar "GzipDetection {} $folder" ::: "${files[@]}"
########################################################
### Let's trim the reads and chomp them if requested ###
########################################################
mkdir -p TMP
echo "Filtering out short sequences and chomping them if needed..."

total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in ${samples[@]}; do
	FileIdentificationInFunction $sample IntGzip
	FileExtractionInFunction IntGzip

	AssemblyPreparation 2> /dev/null

	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done
#parallel -j 4 --bar "AdapterRemovalHeavy {} $folder $out" ::: "${samples[@]}" # Hard coded because of all the pipes...
rm TMP/*
rm -rf IntGzip # Not needed anymore

###############################################
###Let's filter by length and trim the reads###
###############################################
printf "\nRemoving adaptemers and creating the assembly\n"
mkdir -p ${out}SPAdesLogs
total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in ${samples[@]}; do
	AdapterRemovalHeavy $sample ${out}LengthFiltered $out

	echo "Making the Assembly: ML $len & CHOMP $chomp"
	SPAdesAncientFunction $sample ${out}AdaptersFiltered > ${out}/SPAdesLogs/${sample}.log

	if [ $? == 0 ]; then
		printf "\nRunning Quast\n"
		if [ $HOSTNAME == "info2020" ]; then
			python2 /usr/local-centos6/quast/version_3.1/quast.py -o ${out}/$sample/quast ${out}/$sample/contigs.fasta 
		else
			python2 /usr/local/quast/version_3.1/quast.py -o ${out}/$sample/quast ${out}/$sample/contigs.fasta 
		fi
	else
		printf "\n$sample encountered an error...\n"
	fi

	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done
rm -rf TMP
printf "\n"
