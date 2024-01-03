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
usage() { printf 'SPAdes Ancient Assembly Script V0.5
	Second version of assembly script. Follows what I
	have found to be useful and successful
        -i\tThe folder that contains the Raw sequencing files
	-o\tThe output prefix (Default: Assembly)
	-n\tNumber of CPU Threads to be used (Default: 16)
	-m\tUse metaSpades
	-k\tMin Length (Default: 30)
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
	metaSpades:\t${meta}
	-------------------------------------\n"; exit 0;
}

out="Assembly"
ncores=16
len=30
log="$(date +'%Y%m%d').log"
meta=FALSE

while getopts "i:n:o:l:c:k:hm" arg; do
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
	#	c)
	#		chomp=${OPTARG}
	#		#export chomp
	#		;;
		k)
			len=${OPTARG}
			#export len
			;;
		m)
			meta=TRUE
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
printf "\nCreating the assembly\n"
mkdir -p ${out}/SPAdesLogs
total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in ${samples[@]}; do
	#printf "\n$sample\n"
	FileIdentification $sample
	FileExtraction
	#printf "Forward:\t$r1\nReverse:\t$r2\nMerged:\t$merged\n"
	if [ $meta == TRUE ]; then
		spades --meta -1 $r1 -2 $r2 --merged $merged -t $ncores -o ${out}/$sample -k 15,21,29 > $out/SPAdesLogs/${sample}.log  # The kmers are what let it run with the smaller read lengths
	else
		spades -1 $r1 -2 $r2 --merged $merged -t $ncores -o ${out}/$sample -k 15,21,29 > $out/SPAdesLogs/${sample}.log  # The kmers are what let it run with the smaller read lengths
	fi

	if [ $? == 0 ]; then
		printf "\nRunning Quast\n"
		if [ $HOSTNAME == "info2020" ]; then
			python2 /usr/local-centos6/quast/version_3.1/quast.py --min-contig 200 -o ${out}/$sample/quast ${out}/$sample/contigs.fasta 
		else
			python2 /usr/local/quast/version_3.1/quast.py --min-contig 200 -o ${out}/$sample/quast ${out}/$sample/contigs.fasta 
		fi
	else
		printf "\n$sample encountered an error...\n"
	fi

	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done
rm -rf TMP
printf "\n"
