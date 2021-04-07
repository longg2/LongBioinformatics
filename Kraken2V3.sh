# TODO:
# • Bracken Option as well (CAN BE PARALLELIZED EASILY)
# Functions that do the heavy lifting.  Will also allow me to parallelize to my hearts content
FileIdentificationInFunction(){ # Since I can't export Functions, this here is a workaround using find
	local sample=$1 # Basename of the file
	local location=$2

	# Finding the indices which have the same samplename
	sampleFiles=(find $location -name ${sample} -type f | basename)
	#sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))
}

FileIdentification(){ # Extract Files from an array using results from another array
	local sample=$1 # Basename of the file
	local arrayFiles=$files

	# Finding the indices which have the same samplename
	sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))
}

StringDeduplication(){ # Convenient Wrapper for parallelization
	local sample=$1
	local folder=$2
	FileIdentificationInFunction $sample $folder
	FileExtraction $folder
	
	# Now for the actual deduplication
	if [ "$merged" != "NA" ]; then # If I found a merged file
		zcat $merged | perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq stdin -out_bad null -out_good stdout -min_len $len -derep 14 -log ${out}PrinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20 2> /dev/null | gzip > ${out}PooledLanes/${sample}_Merged.fastq.gz
	fi

	if [ "$r1" != "NA" ]; then # If I found a paired file
		perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq <(zcat $r1) -fastq2 <(zcat $r2) -out_bad null -out_good TMP/${sample}_r -min_len $len -derep 14 -log ${out}PrinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20 2> /dev/null

	# Because of how prinseq is coded, I'll need to compress separately	
		gzip -c ${sample}_r1.fastq > ${out}PooledLanes/${sample}_r1.fastq.gz
		gzip -c ${sample}_r2.fastq > ${out}PooledLanes/${sample}_r2.fastq.gz

	fi
}

FileExtraction(){ # Assign files to their variables.  Assumes that $sampleFiles and $sample exists
	# Unsetting variables in case they're already defined from a previous run 
	# unset merged
	# unset r1
	# unset r2

	# This here is a fix that's needed if -v doesn't
	merged="NA"
	r1="NA"
	r2="NA"

	# Creating a hidden text file of file names
	printf '%s\n' "${sampleFiles[@]}" > .hiddenlist.list

	# Identifying the files
	if grep -P -i -q "r1" .hiddenlist.list; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r1')
	        r1="$folder/$fileName"
	fi

	if grep -P -i -q "r2" .hiddenlist.list; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r2')
	        r2="$folder/$fileName"
	fi

	# Two cases for the merged.  Want to control for shenanigans
	if grep -q -i "$sample\.f.*" .hiddenlist.list; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -v -q "_\.f*") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -i "$sample\.f.*")
	        merged="$folder/$fileName"
	fi

	if grep -P -i -q "merged" .hiddenlist.list; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'merged')
	        merged="$folder/$fileName"
	fi

	rm .hiddenlist.list
}
DeduplicateArray(){ # Deduplicating the sample names
	local array=("$@") # This feels weird, but, it'll allow you pass an array to it

	local sampleNames=$(for name in ${array[@]}; do
		tmp="$(echo $name |sed -e 's/_r1.*//I' -e 's/_r2.*//I' -e 's/_merged.*//I' -e 's/\.fa.*//I')";
		echo $tmp;
       	done) # Getting only the sample names
	samples=( $(echo ${sampleNames[@]} | tr ' ' '\n' | uniq | tr '\n' ' ') )
}

GzipDetection(){ # Need to ID Gzipped files and decompress if needed
	local file=$1
	local folderV2=$2
	local name=${file/%.gz*}

	if file $folderV2/$file | grep -q "compressed"; then
		gunzip -c $folderV2/$file > IntGzip/$name
	else
		cp $folderV2/$file IntGzip/$name
	fi
}

KrakenAnalysis(){ # This here will do the heavy lifting for the script

	local sample=$1 # Need to get the sample name
	local folder=$2 # Here we'll get the StringDedup Folder
	FileIdentificationInFunction $sample $folder
	FileExtraction $folder

	# Now to actually run Kraken
	if [ "$merged" != "NA" ]; then # If I found a merged file
		kraken2 --db $db --threads $ncores --report ${out}Reports/${sample}Merged.tab \
			--classified-out ${out}Class/${sample}Merged.fastq --unclassified-out ${out}Unclass/${sample}Merged.fastq \
		       --use-names $merged --output ${out}Output/${sample}Merged.out
	fi
	
	if [ "$r1" != "NA" ]; then # If I found a merged file
		kraken2 --db $db --threads $ncores --report ${out}Reports/${sample}Paired.tab \
			--classified-out ${out}Class/${sample}Paired.fastq --unclassified-out ${out}Unclass/${sample}Paired.fastq \
		       --use-names $r1 --output ${out}Output/${sample}Paired.out
	fi
		
}

TaxaExtraction(){ # If requested, we need to get the reads which mapped against a particular taxa
	local sample=$1 # Need to get the sample name
	local folder=$2 # Here we'll get the StringDedup Folder
	FileIdentificationInFunction $sample $folder # This is all being done as we need to ensure that we're doing both paired and merged reads
	FileExtraction $folder

	if [ "$merged" != "NA" ]; then # If I found a merged file
		extract_kraken_reads.py -k ${out}Output/${sample}Merged.out -r ${out}Reports/${sample}Merged.tab \
			-o ${out}TaxaInterest/$sample.fastq -t $taxa2 --include-children --fastq-output \
			-s $merged
	fi

	if [ "$r1" != "NA" ]; then # If I found a merged file
			extract_kraken_reads.py -k ${out}Output/${sample}Paired.out -r ${out}Reports/${sample}Paired.tab \
				-o ${out}TaxaInterest/${sample}_r1.fastq -o2 ${out}TaxaInterest/${sample}_r2.fastq \
			       	-t $taxa2 --include-children --fastq-output \
				-s $r1 -s2 $r2
	fi

}

usage() { printf 'Kraken2 Wrapper Script V0.5
	Will pipe Kraken2 to the correct folders and create Krona plots.
        -i\tWhere the sequences are
	-o\tOutput Prefix
	-d\tKraken2 Database
        -n\tNumber of CPU Threads to be used
	-l\tLog File Name
	-k\tMinimun Fragment Length
	-s\tString Deduplication? (Default = FALSE)
	-t\tTaxa Being Extracted (Comma Separated list)? (Default = NONE)
	-z\tAre the files gzipped? (Default = FALSE)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }
log() {	printf "Kraken settings for $(date):
	Log File:\t${log}
	Input folder:\t${in}
	Output folder:\t${outPrefix}
	Database:\t${db}
	CPU Threads:\t${ncores}
	Gzip:\t${zipped}
	-------------------------------------
	Min Length:\t${len}
	String Deduplication:\t${dedup}
	Taxa Extraction:\t${taxa}\n"; exit 0;
}

export -f GzipDetection
export -f FileIdentificationInFunction
export -f FileExtraction
export -f StringDeduplication
export -f KrakenAnalysis
export -f TaxaExtraction

declare -i len=30
log="$(date +'%Y%m%d').log"
declare -i ncores=8
export $ncores
out="Kraken2"
export $out
taxa="NULL"
while getopts "i:d:n:o:l:t:k:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
                        #echo "The raw sequencing files are located in $in"
                        ;;
                d)
                        declare -r db=${OPTARG}
			export $db
                        #echo "Using the Kraken2 database at $db"
                        ;;
                n)
                        declare -i ncores=${OPTARG}
			export $ncores
                        #echo "Using $ncores threads"
                        ;;
		o)
			out=${OPTARG}
			export $out
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
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done
log | tee $log

# Making the Directories
mkdir -p ${out}Reports
mkdir -p ${out}CombinedReports
mkdir -p ${out}UnClass
mkdir -p ${out}Class
mkdir -p ${out}Krona
mkdir -p ${out}Output
mkdir -p prinseqLog
#mkdir -p WorkFolder
#mkdir -p FilteredData


# Getting the Sample Names
DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

# We need to determine if the file is gzipped
#echo "Decompressing the files"
#mkdir -p IntGzip
#parallel -j $ncores --bar "GzipDetection {} $folder" ::: "${files[@]}"

# Now to string deduplicate the files as I'd like to speed up the blast runs
mkdir -p ${out}StringDedup
mkdir -p ${out}prinseqLog
echo "String Deduplciation with prinseq"
parallel -j $ncores --bar "StringDeduplication {} $folder" ::: "${samples[@]}"

# Now for the Kraken Loop
for sample in "${samples[@]}";
do
	KrakenAnalysis $sample ${out}StringDedup # This here will do both Paired and unpaired runs.  Can't parallelize it either so it's the rate limiting step as well
done

# Now to combine the reports and create the Krona Plot
parallel --bar -j $ncores "combine_kreports.py -r ${out}Reports/{} --only-combined --no-header -o ${out}CombinedReports/{}.tab > /dev/null 2> /dev/null;
			kreport2krona.py -r ${out}CombinedReports/{}.tab -o ${out}Krona/{}.txt" ::: "${samples[@]}"
ktImportText -o ${outPrefix}Krona.html ${outPrefix}Krona/* # making the KronaPlot

# If requested, we want to also pull out the taxa of interest
if [ $taxa != "NULL" ];then
	mkdir -p ${out}TaxaInterest
	taxa2="$(echo $taxa | sed -e 's/,/ /g')"
	export $taxa2

	# Now to do the actual work
	parallel --bar -j $ncores "TaxaInterest {} ${out}StringDedup" ::: "${samples[@]}"
fi
