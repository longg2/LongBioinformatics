#! /usr/bin/env bash

######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)
export script_full_path

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/MappingFunctions.sh # The mapping functions

BCFtoolsWrapper() { # This is what's going to call my snps
	local BamFile=$1 # The bam file
	local Reference=$2 # The reference file
	local depth=$3
	local qual=$4
	local ploidy=$5

	bcftools mpileup -f $Reference -O u $BamFile |\
		bcftools call --ploidy $ploidy -O u -v -m |\ # Calling the SNPs
		bcftools filter -e "QUAL<${qual} || DP <${depth}" -O v |\ # Filtering out the low quality and depth SNPs
		bcftools norm -f $Reference # Normalizing the variants and outputting the results to STDOUT
}

usage() { printf "BCFtools SNP calling V0.1
	Standardizing the process in which I call SNPs
	Will do the following: mpileup -> call -> filter -> norm

	-i\tThe folder containing the BAM files (REQUIRED)
	-o\tOutput folder of the BAM files (Default: $out)
	-r\tGenome/fasta file we are calling SNPs against (REQUIRED)
	-d\tMinimum SNP depth (Default: 10)
	-q\tMinimum SNP Quality (Default: 100)
	-p\tPloidy (Default: $ploidy)
	-n\tNumber of CPU Threads to be used when parallelizing (Default: 8)
	-l\tLog File Name (Default: $log)
        -h\tShow this help message and exit\n" 1>&2; exit 1; }

log() {	printf "BCFtools SNP settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Reference:\t${ref}
	CPU Threads:\t${ncores}
	-------------------------------------
	Min Quality:\t${qual}
	Min Depth:\t${depth}
	Ploidy:\t${ploidy}
	-------------------------------------\n"; exit 0;
}

######################
### Default values ###
######################

declare -i qual=30
declare -i depth=30
declare -i ploidy=1
out="BCFToolsOut"
declare -i ncores=8
log="$(date +'%Y%m%d').log"

##############
### The UI ###
##############

while getopts "i:o:p:q:r:l:d:n:hm" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find -L $folder/* -printf "%f\n") # Making an array of files
                        ;;
                o)
                        out=${OPTARG}
                        ;;
		p)
			declare -i ploidy=${OPTARG}
                        ;;
		q)
			declare -i qual=${OPTARG}
                        ;;
                r)
                        declare -r ref=${OPTARG}
                        ;;
                l)
                        log=${OPTARG}
                        ;;
                d)
                        declare -i len=${OPTARG}
                        ;;
                n)
                        declare -i ncores=${OPTARG}
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

# The Folder
mkdir -p ${out}

# Need to check if the file has been indexed
if [[ -n $(find . -path "${ref}.fai") ]]; then # This is testing if we're not finding the fai file.
	echo "${ref}.fai has been found!"
else
	echo "${ref} has not been indexed. Running: samtools faidx ${ref}"
	samtools faidx ${ref} # Indexing the file!
fi


export -f BCFtoolsWrapper
# Now to actually run this
parallel --dry-run -j 1 --bar "BCFtoolsWrapper {} $ref $depth $qual $ploidy > ${out}/{/.}.vcf " ::: ${folder}/*
