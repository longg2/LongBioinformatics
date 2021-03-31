# TODO:
# • Make it function focussed
# • Parallelize as much as possible
# • Automate GUNZIP and Sequence Dedup
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
	FileIdentification $1 $2
	FileExtraction $2
	
	# Now for the actual deduplication

	if [ "$merged" != "NA" ]; then # If I found a merged file
		zcat $merged | perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq stdin -out_bad null -out_good stdout -min_len $len -derep 14 -log ${out}PrinseqLog/${sample}Merged.log 2> /dev/null | gzip > ${out}PooledLanes/${sample}_Merged.fastq.gz
	fi

	if [ "$r1" != "NA" ]; then # If I found a paired file
		perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq <(zcat $r1) -fastq2 <(zcat $r2) -out_bad null -out_good TMP/${sample}_r -min_len $len -derep 14 -log ${out}PrinseqLog/${sample}Paired.log 2> /dev/null

	# Because of how prinseq is coded, I'll need to compress separately	
		gzip -c ${sample}_r1.fastq > ${out}PooledLanes/${sample}_r1.fastq.gz
		gzip -c ${sample}_r2.fastq > ${out}PooledLanes/${sample}_r2.fastq.gz

		rm -rf TMP/*
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
		tmp="$(echo $name |sed -e 's/_r1.*//I' -e 's/_r2.*//I' -e "s/_merged.*//I")";
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

export -f GzipDetection
export -f FileIdentificationInFunction
export -f FileExtraction

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

len=30
log="$(date +'%Y%m%d').log"
ncores=8
outPrefix="Kraken2"
dedup="FALSE"
taxa="NULL"
zipped="FALSE"
while getopts "i:d:n:o:l:t:k:hsz" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
                        #echo "The raw sequencing files are located in $in"
                        ;;
                d)
                        db=${OPTARG}
                        #echo "Using the Kraken2 database at $db"
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        #echo "Using $ncores threads"
                        ;;
		o)
			outPrefix=${OPTARG}
			#echo "Outputing to $outPrefix"
			;;
		l)
                        log=${OPTARG}
			;;
		k)
                        len=${OPTARG}
			;;
		s)
			dedup="TRUE"
			;;
		t)
			taxa=${OPTARG}
			;;
		z)
			zipped="TRUE"
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
mkdir -p WorkFolder
mkdir -p FilteredData

# Preparing taxa names (if needed)
if [ $taxa != "NULL" ];then
	mkdir -p ${out}TaxaInterest
	taxa2="$(echo $taxa | sed -e 's/,/ /g')"
fi

# Getting the Sample Names
DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

# We need to determine if the file is gzipped
echo "Decompressing the files"
mkdir -p IntGzip
parallel -j $ncores --bar "GzipDetection {} $folder" ::: "${files[@]}"

# Now to string deduplicate the files as I'd like to speed up the blast runs
mkdir -p ${out}StringDedup/
mkdir -p ${out}prinseqLog
echo "String Deduplciation with prinseq"
parallel -j $ncores --bar "perl /home/sam/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fasta {} -out_good ${out}StringDedup/{/.} -out_bad null -min_len $len -derep 14 -log ${out}prinseqLog/{/.}.log 2> /dev/null" ::: ${out}FastaOnly/*

# Doing the 
for sample in "${samples[@]}"; # The loop to analyze the paired and unpaired files
do
	echo "$sample being analyzed"
	if [ $zipped == "TRUE" ]; then
		merged="$in/$sample.fastq.gz"
		paired="$in/${sample}_r1.fastq.gz"
	else
		merged="$in/$sample.fastq"
		paired="$in/${sample}_r1.fastq"
	fi
	#echo $paired
	if [ -f  "$merged" ]; then # This is for the merged samples
		echo "Merged file found"

		if [ $zipped == "TRUE" ] ; then
			zcat $merged > WorkFolder/${sample}.fastq
		else
			cat $merged >  WorkFolder/${sample}.fastq
		fi

		# This is affected by the option s
		if [ $dedup = "TRUE" ];
		then
			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq WorkFolder/${sample}.fastq -out_good FilteredData/${sample} -derep 14 -min_len $len -log prinseqLog/${sample}Merged.log # Need to dedup
		else
			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq WorkFolder/${sample}.fastq -out_good FilteredData/${sample} -min_len $len -log prinseqLog/${sample}Merged.log # Need to dedup
		fi

	else
		echo "No Merged File Found"
	fi
	
	if [ -f "$paired" ]; then # Are there PE files
		echo "Paired file found"

		if [ $zipped == "TRUE" ]; then
			zcat $paired > WorkFolder/${sample}_r1.fastq
			zcat ${paired/r1/r2} > WorkFolder/${sample}_r2.fastq
		else
			cat $paired > WorkFolder/${sample}_r1.fastq
			cat ${paired/r1/r2} > WorkFolder/${sample}_r2.fastq

		fi

		# This is affected by the option s
		if [ $dedup = "TRUE" ];then
			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq WorkFolder/${sample}_r1.fastq -fastq2 WorkFolder/${sample}_r2.fastq -out_good FilteredData/${sample} -derep 14 -min_len $len -log prinseqLog/${sample}Paired.log# Need to dedup
		else
			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq WorkFolder/${sample}_r1.fastq -fastq2 WorkFolder/${sample}_r2.fastq -out_good FilteredData/${sample} -min_len $len -log prinseqLog/${sample}Paired.log# Need to dedup
		fi

	else
		echo "No Paired File Found"
	fi

done

# Now for the Kraken Loop
for sample in "${samples[@]}"; # The loop to analyze the paired and unpaired files
do
	merged="FilteredData/$sample.fastq"
	paired="FilteredData/${sample}_1.fastq"

	if [ -f  "$merged" ]; then # This is for the merged samples
		echo "Paired Found"
		kraken2 --db $db --threads $ncores --report ${outPrefix}Reports/${sample}Merged.tab \
			--classified-out ${outPrefix}Class/${sample}Merged.fastq --unclassified-out ${outPrefix}Unclass/${sample}Merged.fastq \
		       --use-names $merged --output ${outPrefix}Output/${sample}Merged.out
	fi
	
	if [ -f "$paired" ]; then # Are there PE files
		echo "Paired Found"
		kraken2 --db $db --threads $ncores --report ${outPrefix}Reports/${sample}Paired.tab \
			--classified-out ${outPrefix}Class/${sample}Paired.fastq --unclassified-out ${outPrefix}Unclass/${sample}Paired.fastq \
		       --use-names $paired --output ${outPrefix}Output/${sample}Paired.out
	fi
done

# Now to extract the Reads matching taxa of interest (If wanted)
if [ $taxa != "NULL" ]; then
#	parallel -j $ncores --bar "extract_kraken_reads.py -k ${outPrefix}Output/{}Merged.out -r ${outPrefix}Reports/{}Merged.tab \
#			-o ${outPrefix}TaxaInterest/{}.fastq -t $taxa2 --include-children --fastq-output \
#			-s FilteredData/{}.fastq;
#	
#			extract_kraken_reads.py -k ${outPrefix}Output/{}Paired.out -r ${outPrefix}Reports/{}Paired.tab \
#				-o ${outPrefix}TaxaInterest/{}_r1.fastq -o2 ${outPrefix}TaxaInterest/{}_r2.fastq \
#			       	-t $taxa2 --include-children --fastq-output \
#				-s FilteredData/{}_r1.fastq -s2 FilteredData/{}_r2.fastq" ::: "${samples[@]}"
	for sample in "${samples[@]}"; # The loop to analyze the paired and unpaired files
	do
		# Merged
		extract_kraken_reads.py -k ${outPrefix}Output/${sample}Merged.out -r ${outPrefix}Reports/${sample}Merged.tab \
			-o ${outPrefix}TaxaInterest/$sample.fastq -t $taxa2 --include-children --fastq-output \
			-s FilteredData/$sample.fastq

		# Paired
			extract_kraken_reads.py -k ${outPrefix}Output/${sample}Paired.out -r ${outPrefix}Reports/${sample}Paired.tab \
				-o ${outPrefix}TaxaInterest/${sample}_r1.fastq -o2 ${outPrefix}TaxaInterest/${sample}_r2.fastq \
			       	-t $taxa2 --include-children --fastq-output \
				-s FilteredData/${sample}_1.fastq -s2 FilteredData/${sample}_2.fastq
	done
fi

# Combining the reports
for sample in "${samples[@]}"; do
	combine_kreports.py -r ${outPrefix}Reports/$sample* --only-combined --no-headers -o ${outPrefix}CombinedReports/$sample.tab # Combining the reports
	kreport2krona.py -r ${outPrefix}CombinedReports/$sample.tab -o ${outPrefix}Krona/$sample.txt # Preparing for a Krona Plot
done
#parallel -j $ncores --bar "combine_kreports.py -r ${outPrefix}Reports/{}* --only-combined --no-headers -o ${outPrefix}CombinedReports/{}.tab;
#			kreport2krona.py -r ${outPrefix}CombinedReports/{}.tab -o ${outPrefix}Krona/{}.txt" :::"${samples[@]}"

ktImportText -o ${outPrefix}Krona.html ${outPrefix}Krona/* # making the KronaPlot
