#! /usr/bin/env bash
# TODO:
# • Bracken Option as well (CAN BE PARALLELIZED EASILY)
# Functions that do the heavy lifting.  Will also allow me to parallelize to my hearts content

script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/QCFunctions.sh # Contains StringDedup
source $script_full_path/lib/KrakenFunctions.sh # Contains the KrakenFunctions


#FileIdentificationInFunction(){ # Since I can't export Functions, this here is a workaround using find
#	local sample=$1 # Basename of the file
#	local location=$2
#
#	# Finding the indices which have the same samplename
#	sampleFiles=( $(find $location -name "${sample}*" -type f -exec basename {} \;) )
#	#sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))
#}
#FileIdentification(){ # Extract Files from an array using results from another array
#	local sample=$1 # Basename of the file
#	local arrayFiles=$files
#
#	# Finding the indices which have the same samplename
#	sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))
#}
#FileExtractionInFunction(){ # Assign files to their variables.  Assumes that $sampleFiles and $sample exists
#	# Unsetting variables in case they're already defined from a previous run 
#	# unset merged
#	# unset r1
#	# unset r2
#
#	local folderExtraction=$1
#	# This here is a fix that's needed if -v doesn't
#	merged="NA"
#	r1="NA"
#	r2="NA"
#
#	# Creating a hidden text file of file names
#	#printf '%s\n' "${sampleFiles[@]}" #> .hiddenlist.list
#
#	# Identifying the files
#	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r1"; then
#		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r1')
#	        r1="$folderExtraction/$fileName"
#	fi
#
#	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r2"; then
#		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r2')
#	        r2="$folderExtraction/$fileName"
#	fi
#
#	# Two cases for the merged.  Want to control for shenanigans
#	if printf '%s\n' "${sampleFiles[@]}" | grep -q -i "$sample\.f.*"; then
#	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -v -q "_\.f*") ]; then
#		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -i "$sample\.f.*")
#	        merged="$folderExtraction/$fileName"
#	fi
#
#	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged"; then
#	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged") ]; then
#		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'merged')
#	        merged="$folderExtraction/$fileName"
#	fi
#
#}
#FileExtraction(){ # Assign files to their variables.  Assumes that $sampleFiles and $sample exists
#	# Unsetting variables in case they're already defined from a previous run 
#	# unset merged
#	# unset r1
#	# unset r2
#
#	# This here is a fix that's needed if -v doesn't
#	merged="NA"
#	r1="NA"
#	r2="NA"
#
#	# Creating a hidden text file of file names
#	printf '%s\n' "${sampleFiles[@]}" #> .hiddenlist.list
#
#	# Identifying the files
#	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r1"; then
#		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r1')
#	        r1="$folder/$fileName"
#	fi
#
#	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r2"; then
#		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r2')
#	        r2="$folder/$fileName"
#	fi
#
#	# Two cases for the merged.  Want to control for shenanigans
#	if printf '%s\n' "${sampleFiles[@]}" | grep -q -i "$sample\.f.*"; then
#	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -v -q "_\.f*") ]; then
#		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -i "$sample\.f.*")
#	        merged="$folder/$fileName"
#	fi
#
#	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged"; then
#	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged") ]; then
#		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'merged')
#	        merged="$folder/$fileName"
#	fi
#
#}
#DeduplicateArray(){ # Deduplicating the sample names
#	local array=("$@") # This feels weird, but, it'll allow you pass an array to it
#
#	local sampleNames=$(for name in ${array[@]}; do
#		tmp="$(echo $name |sed -e 's/_r1.*//I' -e 's/_r2.*//I' -e 's/_merged.*//I' -e 's/\.fq.*//I' -e 's/\.fa.*//I')";
#		echo $tmp;
#       	done) # Getting only the sample names
#	samples=( $(echo ${sampleNames[@]} | tr ' ' '\n' | uniq | tr '\n' ' ') )
#}
#GzipDetection(){ # Need to ID Gzipped files and decompress if needed
#	local file=$1
#	local folderV2=$2
#	local name=${file/%.gz*}
#
#	if file $folderV2/$file | grep -q "compressed"; then
#		gunzip -c $folderV2/$file > IntGzip/$name
#	else
#		cp $folderV2/$file IntGzip/$name
#	fi
#}
#KrakenAnalysis(){ # This here will do the heavy lifting for the script
#
#	local sample=$1 # Need to get the sample name
#	local folderFunction=$2 # Here we'll get the StringDedup Folder
#	local ncores=$3 # How many cores will be used?
#	FileIdentificationInFunction $sample $folderFunction
#	FileExtractionInFunction $folderFunction
#
#	# Now to actually run Kraken
#	if [ "$merged" != "NA" ]; then # If I found a merged file
#		kraken2 --db $db --threads $ncores --report ${out}Reports/${sample}Merged.tab \
#			--classified-out ${out}Class/${sample}Merged.fastq --unclassified-out ${out}Unclass/${sample}Merged.fastq \
#		       --use-names $merged --output ${out}Output/${sample}Merged.out --confidence 0.1
#	fi
#
#	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If dealing with a single end library
#		kraken2 --db $db --threads $ncores --report ${out}Reports/${sample}Single.tab \
#			--classified-out ${out}Class/${sample}Single.fastq --unclassified-out ${out}Unclass/${sample}Single.fastq \
#		       --use-names $r1 --output ${out}Output/${sample}Single.out --confidence 0.1
#	fi 
#	
#	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
#		kraken2 --db $db --threads $ncores --report ${out}Reports/${sample}Paired.tab \
#			--classified-out ${out}Class/${sample}Paired.fastq --unclassified-out ${out}Unclass/${sample}Paired.fastq \
#		       --use-names $r1 --output ${out}Output/${sample}Paired.out --confidence 0.1
#	fi
#		
#}
#TaxaExtraction(){ # If requested, we need to get the reads which mapped against a particular taxa
#	local sample=$1 # Need to get the sample name
#	local folderFunction=$2 # Here we'll get the StringDedup Folder
#	FileIdentificationInFunction $sample $folderFunction
#	FileExtractionInFunction $folderFunction
#
#	if [ "$merged" != "NA" ]; then # If I found a merged file
#		extract_kraken_reads.py -k ${out}Output/${sample}Merged.out -r ${out}Reports/${sample}Merged.tab \
#			-o ${out}TaxaInterest/$sample.fastq -t $taxa2 --include-children --fastq-output \
#			-s $merged
#	fi
#
#	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If dealing with a single end library
#		extract_kraken_reads.py -k ${out}Output/${sample}Single.out -r ${out}Reports/${sample}Single.tab \
#			-o ${out}TaxaInterest/${sample}_single.fastq -t $taxa2 --include-children --fastq-output \
#			-s $r1
#	fi 
#
#	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
#			extract_kraken_reads.py -k ${out}Output/${sample}Paired.out -r ${out}Reports/${sample}Paired.tab \
#				-o ${out}TaxaInterest/${sample}_r1.fastq -o2 ${out}TaxaInterest/${sample}_r2.fastq \
#			       	-t $taxa2 --include-children --fastq-output \
#				-s $r1 -s2 $r2
#	fi
#
#}
#StringDeduplication(){ # Convenient Wrapper for parallelization
#	local sample=$1
#	local folderFunction=$2
#	local len=$3
#	FileIdentificationInFunction $sample $folderFunction
#	FileExtractionInFunction $folderFunction
#
#	# Now for the actual deduplication
#	if [ "$merged" != "NA" ]; then # If I found a merged file
#		if [ "$Dedup" == "TRUE" ]; then
#			prinseq -fastq $merged -out_bad null -out_good stdout -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20 | gzip > ${out}StringDedup/${sample}_Merged.fastq.gz
#		else
#			prinseq -fastq $merged -out_bad null -out_good stdout -min_len $len -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20 | gzip > ${out}StringDedup/${sample}_Merged.fastq.gz
#		fi
#	fi
#
#	if [[ $r1 != "NA" && $r2 == "NA" ]]; then # If dealing with a single end library
#		echo "Single"
#		if [ "$Dedup" == "TRUE" ]; then
#			prinseq -fastq $r1 -out_bad null -out_good stdout -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 | gzip > ${out}StringDedup/${sample}_r1.fastq.gz
#		else
#			prinseq -fastq $r1 -out_bad null -out_good stdout -min_len $len -log ${out}prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 | gzip > ${out}StringDedup/${sample}_r1.fastq.gz
#		fi
#	fi 
#
#	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
#		if [ "$Dedup" == "TRUE" ]; then
#			prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20
#		else
#			prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20
#
#		fi
#
#	# Because of how prinseq is coded, I'll need to compress separately	
#		gzip -c TMP/${sample}_1.fastq > ${out}StringDedup/${sample}_r1.fastq.gz
#		gzip -c TMP/${sample}_2.fastq > ${out}StringDedup/${sample}_r2.fastq.gz
#
#	fi
#}
#StringDeduplicationParallel(){ # Convenient Wrapper for parallelization
#	local sample=$1
#	local folderFunction=$2
#	local len=$3
#	FileIdentificationInFunction $sample $folderFunction
#	FileExtractionInFunction $folderFunction
#
#	# Now for the actual deduplication
#	if [ "$merged" != "NA" ]; then # If I found a merged file
#		if [ "$Dedup" == "TRUE" ]; then
#			prinseq -fastq $merged -out_bad null -out_good TMP/${sample} -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20 
#		else
#			prinseq -fastq $merged -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20
#		fi
#	fi
#
#	if [[ $r1 != "NA" && $r2 == "NA" ]]; then # If dealing with a single end library
#		if [ "$Dedup" == "TRUE" ]; then
#			prinseq -fastq $r1 -out_bad null -out_good TMP/${sample} -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 
#		else
#			prinseq -fastq $r1 -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 
#		fi
#	fi 
#
#	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
#		if [ "$Dedup" == "TRUE" ]; then
#			prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20
#		else
#			prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20
#
#		fi
#
#	# Because of how prinseq is coded, I'll need to compress separately	
#		gzip -c TMP/${sample}.fastq > ${out}StringDedup/${sample}.fastq.gz
#		gzip -c TMP/${sample}_1.fastq > ${out}StringDedup/${sample}_r1.fastq.gz
#		gzip -c TMP/${sample}_2.fastq > ${out}StringDedup/${sample}_r2.fastq.gz
#
#	fi
#}
#ProgressBar() { # From github.com/fearside/ProgressBar
#	# Process data
#		let _progress=(${1}*100/${2}*100)/100
#		let _done=(${_progress}*4)/10
#		let _left=40-$_done
#	# Build progressbar string lengths
#	_done=$(printf "%${_done}s")
#	_left=$(printf "%${_left}s")
#	printf "\rProgress : [${_done// />}${_left// /-}] ${_progress}%%"
#
#}	

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
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
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
                        declare -i len=${OPTARG}
			export $len
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
#echo "Decompressing the files"
#mkdir -p IntGzip
#if [ $ncores > 30  ];
#then
#	#Preventing an accidental swamping of the cluster
#	parallel -j 30 --bar "GzipDetection {} $folder" ::: "${files[@]}"
#else
#	parallel -j $ncores --bar "GzipDetection {} $folder" ::: "${files[@]}"
#fi

# Now to string deduplicate the files as I'd like to speed up the blast runs
mkdir -p ${out}StringDedup
mkdir -p ${out}prinseqLog
mkdir -p TMP
echo "String Deduplciation with prinseq"

total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in ${samples[@]}
do
	StringDeduplication $sample IntGzip $len
	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done
printf "\n"
#parallel -j $ncores --bar "StringDeduplicationParallel {} IntGzip $len" ::: "${samples[@]}" # DOESN'T WORK!!
rm -rf TMP
#rm -rf IntGzip

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
