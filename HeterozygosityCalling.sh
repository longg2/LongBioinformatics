#! /usr/bin/env bash
# TODO: PE reads need to be properly handled

script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
usage() { printf "Heterozygous SNP Calling V1.0
	Given a folder of mapped sequences, call SNPs and figure out which
	ones could be heterozygous. This script will allow you to choose
	either a preset ploidy or supply a custom ploidy file.
	
        -i\tInput folder
	-o\tOutput Prefix (Default = HetSNPs)
	-r\tReference Sequence
	-n\tNumber of CPU Threads to be used (Default = 10)
	-l\tLog File Name
	-s\tSeparate VCF Files?
	-o\tBecomes a folder name (Default = FALSE)
	-q\tMinimum SNP Quality calling (Default = 100)
	-p\tPloidy of the call. For possible options, see bcftools call --ploidy ? Overwritten by -f
	-f\tA ploidy file in the form of CHROM<TAB>FROM<TAB>TO<TAB>SEX<TAB>PLOIDY Overwrites -p
        -h\tShow this help message and exit\n" 1>&2; exit 0; }
log() {	printf "Heterozygous SNPs settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output:\t${out}
	Separate Files:\t${separate}
	Reference:\t${reference}
	Minimum Quality:\t${minqual}
	Ploidy:\t${ploidy}
	Ploidy File:\t${ploidyFile}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

log="$(date +'%Y%m%d').log"
declare -i ncores=10
export out="HetSNPs"
minqual=100
ploidy=2
separate=FALSE
while getopts "i:n:o:p:q:r:s:l:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
                        ;;
                r)
                        declare -r reference=${OPTARG}
			export reference
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        ;;
		o)
			export out=${OPTARG}
			;;
		p)
			ploidy=${OPTARG}
			;;
		f)
			ploidyFile=${OPTARG}
			unset ploidy
			;;
		l)
                        log=${OPTARG}
			;;
		s)
			separate=TRUE
			;;
		q)
			export minqual=${OPTARG}
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Testing if there's anything here
if [ -z ${folder+x} ] || [ -z ${reference+x} ]; then
	printf "You are missing either the Input folder or the Reference Genome\n\n"
	usage
	exit 1
fi


nthreads=$((($ncores - $ncores % 3) / 3)) # Want multithreading, but, since I'm piping need to be somewhat reasonable here

# Need to test if we want all the results in one file or separately

if [ $separate == "FALSE" ]; then
	if [ -z ${ploidy+x} ]; then # Testing if ploidy is unset
	
		bcftools mpileup -Ou -f $reference $folder/* --threads $nthreads \|
		       bcftools call -Ou -mv --ploidy-file $ploidyFile --threads $nthreads \|
		       bcftools filter -s LowQual -e "%QUAL<$minqual" -Oz6 --threads $nthreads > ${out}.vcf.gz
	
	else
		bcftools mpileup -Ou -f $reference $folder/* --threads $nthreads \|
			bcftools call -Ou -mv --ploidy $ploidy --threads $nthreads \|
		       	bcftools filter -s LowQual -e "%QUAL<$minqual" -Oz6 --threads $nthreads > ${out}.vcf.gz
	
	fi
else

	# Setting up the loop
	mkdir -p ${out}VCFFiles
	total=$(find $folder/ -type f | grep -cE "\.sam$|\.bam$")
	count=0
	
	ProgressBar $count $total
	if [ -z ${ploidy+x} ]; then # Testing if ploidy is unset
	
		for file in Input/*{bam,sam}; do

			# Preparing the output name
			sample=$(basename $file)
			sample=${sample%.*}
			# Actually running the command
			bcftools mpileup -Ou -f $reference $file --threads $nthreads \|
			       bcftools call -Ou -mv --ploidy-file $ploidyFile --threads $nthreads \|
			       bcftools filter -s LowQual -e "%QUAL<$minqual" -Oz6 --threads $nthreads > ${out}VCFFiles/${sample}.vcf.gz

			# Moving the progress bar
			count=$(echo "$count + 1" | bc)
			ProgressBar $count $total
		done
	else
		for file in Input/*{bam,sam}; do

			# Preparing the output name
			sample=$(basename $file)
			sample=${sample%.*}
			# Actually running the command
			bcftools mpileup -Ou -f $reference $file --threads $nthreads \|
			       bcftools call -Ou -mv --ploidy $ploidy --threads $nthreads \|
			       bcftools filter -s LowQual -e "%QUAL<$minqual" -Oz6 --threads $nthreads > ${out}VCFFiles/${sample}.vcf.gz

			# Moving the progress bar
			count=$(echo "$count + 1" | bc)
			ProgressBar $count $total
		done
	fi
fi
