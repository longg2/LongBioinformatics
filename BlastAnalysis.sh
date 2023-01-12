#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/BlastFunctions.sh # This loads the QC commands
usage() { printf "BlastN/P Wrapper Script V0.9
	Outputs tab deliminated BlastN/P report file in the form of std staxid.  Taxa counts
	for each step are also outputted.  Need seqkit for subsampling.	String Deduplication occurs.
	LCA Analysis requires taxonkit <https://github.com/shenwei356/taxonkit>
	TBD: Remove Seqkit Completely
        -i\tInput Folder
	-o\tOutput Folder (Default = BlastOut)
	-m\tUse Diamond for the analysis (Default = FALSE)
	-b\tblastn, blastp, or blastx? (Default = blastn)
	-d\tBlast database (Default = /1/scratch/blastdb/nt)
	-e\tMaximum Evalue (Default = 1e-5)
	-k\tRead Length (Default = 30)
	-p\tPercent Identity (Default = 90)
	-s\tReads to Subsample (Default = 0)
	-t\tLCA Analysis (Default = FALSE)
	-n\tNumber of CPU Threads to be used (Default = 8)
	-l\tLog File Name (Default is date+Ymd)
        -h\tShow this help message and exit\n" 1>&2; exit 1; }
log() {	printf "Blast settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Diamond:\t${diamond}
	Blast Type:\t${blast}
	Blast Database:\t${db}
	Computer:\t${HOSTNAME}
	-------------------------------------
	E-value:\t${eval}
	Percent Identity:\t${Pident}
	Subsampled Reads:\t${subsample}
	Min Read Length:\t${len}	
	LCA:\t${lca}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

# Not a fan of this as it limits it to me only....
source ~/miniconda3/etc/profile.d/conda.sh # Activating 
conda activate base
conda activate pangenome

export -f FastaorFastq # important as I'm running this in parallel
export -f GzipDetection # important as I'm running this in parallel
export -f LCA # important as I'm running this in parallel
export -f RandomFastaSelection # This is to remove my reliance on seqkit
#export -f ProgressBar

# Preset Variables
#export alias blastn="~/miniconda3/envs/pangenome/bin/blastn" # Will it work?
#echo $(which blastn) 
#exit 0
blast="blastn"
if [[ $HOSTNAME == "info113" ]]; then
	db="/1/scratch/blastdb/nt_v5"
else
	db="/1/scratch/blastdb/nt"
fi

declare -i Pident=90
declare -i len=30
declare -i subsample=0
declare -i ncores=8
eval="1e-5"
log="$(date +'%Y%m%d').log"
lca="FALSE"
diamond="FALSE"
out="BlastOut"

while getopts "i:k:e:d:n:o:b:p:l:s:hmt" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find -L $folder/* -type f -printf "%f\n") # Making an array of files
                        ;;
		e)
			eval=${OPTARG}
			;;
                d)
                        db=${OPTARG}
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        ;;
		o)
			out=${OPTARG}
			;;
		p)
			declare -i Pident=${OPTARG}
			;;
		b)
			blast=${OPTARG}
			;;
		s)
			declare -i subsample=${OPTARG}
			;;
		l)
			log=${OPTARG}
			;;
		k)
                        declare -i len=${OPTARG}
			;;
		t)
                        lca="TRUE"
			;;
		m)
			diamond="TRUE"
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

##################
### The Workup ###
##################
log | tee $log

# We need to determine if the file is gzipped
echo "Decompressing the files"
mkdir -p IntGzip
parallel -j $ncores --bar "GzipDetection {} $folder" ::: "${files[@]}"

# Now to detect if the file in IntGzip is fastq, fasta, or other and convert the file
mkdir -p ${out}FastaOnly
echo "Converting FastQ to Fasta"
parallel -j $ncores --bar "FastaorFastq {} ${out}FastaOnly" ::: IntGzip/*
rm -rf IntGzip

# Now to string deduplicate the files as I'd like to speed up the blast runs
mkdir -p ${out}StringDedup

if [[ "${blast}" == "blastp" ]]; then
	echo "BlastP requested.  String deduplication won't be performed"
	cp ${out}FastaOnly/* ${out}StringDedup/
else

	mkdir -p ${out}prinseqLog
	echo "String Deduplciation with prinseq"
	parallel -j $ncores --bar "perl /home/sam/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fasta {} -out_good ${out}StringDedup/{/.} -out_bad null -min_len $len -derep 14 -log ${out}prinseqLog/{/.}.log 2> /dev/null" ::: ${out}FastaOnly/*
fi

# The final step is to test if StringDeduplication will occur
if [[ $subsample != 0 ]]; then
	echo "Subsampling the files to ~$subsample reads"
	mkdir -p ${out}Subsample
	parallel -j $ncores --bar "RandomFastaSelection {} $subsample > ${out}Subsample/{/}" ::: ${out}StringDedup/*
	#parallel -j $ncores --bar "seqkit sample {} -n $subsample > ${out}Subsample/{/} 2> /dev/null" ::: ${out}StringDedup/*
fi

#####################################
### Now to run the Blast commands ###
#####################################

#rm -rf ${out} # Delete the folder if it already exists
mkdir -p ${out}BlastResults
if [[ $subsample != 0 ]]; then
	for final in ${out}Subsample/*; do
		if [[ $diamond == "TRUE" ]]; then
			diamondCMD $final	
		else
			blastCMD $final	
		fi
	done
else
	for final in ${out}StringDedup/*; do
		if [[ $diamond == "TRUE" ]]; then
			diamondCMD $final	
		else
			blastCMD $final	
		fi
	done
fi

#############################
### LCA Script if desired ###
#############################

if [[ $lca == "TRUE" ]]; then
	echo "Determining the LCA for each hit"
	mkdir -p ${out}LCA
	# Because taxonkit defaults to 4 cores, we don't want to accidentally swamp the machines
	let taxonKitcores=$ncores/8 # 8 Because I'm using two instances of taxonkit per script
	parallel -j $taxonKitcores --bar "LCA {}" ::: ${out}BlastResults/*
fi

###############################
### Compressing the Results ###
###############################
echo "Compressing the blast files"

parallel -j $ncores --bar "gzip {}" ::: ${out}BlastResults/*tab

echo "Blast is Finished!"
