#! /usr/bin/env bash
# TODO:
# • Bracken Option as well (CAN BE PARALLELIZED EASILY)
# Functions that do the heavy lifting.  Will also allow me to parallelize to my hearts content

script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/QCFunctions.sh # Contains StringDedup
source $script_full_path/lib/KrakenFunctions.sh # Contains the KrakenFunctions


usage() { printf 'Kraken2 Wrapper Script V0.75
	Will pipe Kraken2 to the correct folders and create Krona plots.
        -i\tWhere the sequences are
	-o\tOutput Prefix (Default = Kraken2)
	-r\tKraken2 Database
	-n\tNumber of CPU Threads to be used (Default = 8)
	-l\tLog File Name
	-k\tMinimum Fragment Length (Default = 30)
	-t\tTaxa Being Extracted (Comma Separated list)? (Default = NONE)
	-d\tString Deduplication? (Default = FALSE)
        -h\tShow this help message and exit\n' 1>&2; exit 0; }
log() {	printf "Kraken settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Database:\t${db}
	CPU Threads:\t${ncores}
	-------------------------------------
	Min Length:\t${len}
	Taxa Extraction:\t${taxa}
	-------------------------------------\n"; exit 0;
}

# Need to export some of the functions
export -f GzipDetection
export -f FileIdentificationInFunction
export -f FileExtractionInFunction
export -f FileExtraction
export -f StringDeduplication
export -f StringDeduplicationParallel
export -f KrakenAnalysis
export -f TaxaExtraction

declare -i len=30
log="$(date +'%Y%m%d').log"
declare -i ncores=8
export out="Kraken2"
taxa="NULL"
Dedup="FALSE"
while getopts "i:r:n:o:l:t:k:hd" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find -L $folder/* -type f -printf "%f\n") # Making an array of files
                        #echo "The raw sequencing files are located in $in"
                        ;;
                r)
                        declare -r db=${OPTARG}
			#export $db
                        #echo "Using the Kraken2 database at $db"
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        #echo "Using $ncores threads"
                        ;;
		o)
			export out=${OPTARG}
			#export $out
			#echo "Outputing to $outPrefix"
			;;
		l)
                        log=${OPTARG}
			;;
		k)
			export len=${OPTARG}
			;;
		t)
			taxa=${OPTARG}
			;;
		d)
			Dedup="TRUE"
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

if [ -z ${folder+x} ] || [ -z ${db+x} ]; then
	echo "You are missing either the Input folder or the Kraken Database"
	usage
	exit 1
fi

log | tee $log
printf "\n"

# Making the Directories
mkdir -p ${out}Reports
mkdir -p ${out}CombinedReports
mkdir -p ${out}UnClass
mkdir -p ${out}Class
mkdir -p ${out}Krona
mkdir -p ${out}Output
mkdir -p ${out}KrakenLog

## Getting the Sample Names
DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

# We need to determine if the file is gzipped
echo "Decompressing the files"
mkdir -p IntGzip
if [ $ncores > 30  ];
then
	#Preventing an accidental swamping of the cluster
	parallel -j 30 --bar "GzipDetection {} $folder" ::: "${files[@]}"
else
	parallel -j $ncores --bar "GzipDetection {} $folder" ::: "${files[@]}"
fi

# Now to string deduplicate the files as I'd like to speed up the blast runs
mkdir -p ${out}StringDedup
mkdir -p ${out}prinseqLog
mkdir -p TMP
echo "String Deduplciation with prinseq"

#total=${#samples[@]}
#count=0
#ProgressBar $count $total
#for sample in ${samples[@]}
#do
#	StringDeduplication $sample IntGzip $len
#	count=$(echo "$count + 1" | bc)
#	ProgressBar $count $total
#done
#printf "\n"
parallel -j $ncores --bar "StringDeduplicationParallel {} IntGzip $len 2> /dev/null" ::: "${samples[@]}" # DOESN'T WORK!!
rm -rf TMP
rm -rf IntGzip

# Now for the Kraken Loop
echo "Running Kraken2 on ${#samples[@]} samples"
total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in "${samples[@]}";
do
	KrakenAnalysis $sample ${out}StringDedup $ncores 2> ${out}KrakenLog/${sample}.log # This here will do both Paired and unpaired runs.  Can't parallelize it either so it's the rate limiting step as well
	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done

# Now to combine the reports and create the Krona Plot

echo "Combining the Merged and Paired Reports and preparing for a Krona plot"
parallel --bar -j $ncores "combine_kreports.py -r ${out}Reports/{}*tab --only-combined --no-header -o ${out}CombinedReports/{}.tab > /dev/null 2> /dev/null; kreport2krona.py -r ${out}CombinedReports/{}.tab -o ${out}Krona/{}.txt" ::: "${samples[@]}"
ktImportText -o KronaPlot.html ${outPrefix}Krona/* # making the KronaPlot
#
## If requested, we want to also pull out the taxa of interest
if [ $taxa != "NULL" ];then
	mkdir -p ${out}TaxaInterest
	export taxa2="$(echo $taxa | sed -e 's/,/ /g')"
	echo "Extracting sequences identified as taxons $taxa2"

	# Now to do the actual work
	parallel --bar -j $ncores "TaxaExtraction {} ${out}StringDedup > /dev/null 2> /dev/null" ::: "${samples[@]}"
fi
