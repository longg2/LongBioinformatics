#! /usr/bin/env bash

######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)
export script_full_path

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/MappingFunctions.sh # The mapping functions

usage() { printf "BCFtools SNP calling V0.1
	Standardizing the process in which I call SNPs
	Will do the following: mpileup -> call -> filter -> norm

	-i\tThe folder containing the BAM files (REQUIRED)
	-o\tOutput folder of the BAM files (Default: $out)
	-r\tGenome/fasta file we are calling SNPs against (REQUIRED)
	-d\tMinimum SNP depth (Default: 10)
	-q\tMinimum SNP Quality (Default: 100)
	-p\tPloidy (Default: $ploidy)
	-l\tLog File Name (Default: $log)
        -h\tShow this help message and exit\n" 1>&2; exit 1; }

log() {	printf "BCFtools SNP settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Reference:\t${ref}
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
log="$(date +'%Y%m%d').log"

##############
### The UI ###
##############

while getopts "i:o:p:q:r:l:d:hm" arg; do
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
mkdir -p ${out}Log

# Need to check if the file has been indexed
if [[ -n $(find $(dirname $ref) -name "$(basename $ref).fai") ]]; then # This is testing if we're not finding the fai file.
	echo "${ref}.fai has been found!"
else
	echo "${ref} has not been indexed. Running: samtools faidx ${ref}"
	samtools faidx ${ref} # Indexing the file!
fi


# Now to actually run this
total=$(ls -1q $folder | wc -l)
count=0

ProgressBar $count $total
for file in $folder/*; do # Iterating over an array of Samples

	#printf "BCFtoolsWrapper $file $ref $depth $qual $ploidy > ${out}/$(basename $file .bam).vcf\n"
	bcftools mpileup -f $ref -O u $file 2> /dev/null |\
		bcftools call --ploidy $ploidy -O u -v -m |\
		bcftools filter -e "QUAL<${qual} || DP <${depth}" -O v |\
		bcftools norm -f $ref > ${out}/$(basename $file .bam).vcf 2> ${out}Log/$(basename $file .bam).log
	#BCFtoolsWrapper $file $ref $depth $qual $ploidy > ${out}/$(basename $file .bam).vcf

	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done
printf "\nComplete!\n"
