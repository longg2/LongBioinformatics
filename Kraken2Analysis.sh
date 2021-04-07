usage() { printf 'Kraken2 Wrapper Script V0.5
	Will pipe Kraken2 to the correct folders and create Krona plots.
        -i\tWhere the sequences are
	-o\tOutput Prefix
	-d\tKraken2 Database
	-s\tString Deduplication
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
                        in=${OPTARG}
                        #echo "The sequences are in $in"
                        ;;
                d)
                        db=${OPTARG}
                        #echo "Using the Kraken2 database at $db"
                        ;;
                n)
                        ncores=${OPTARG}
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
mkdir -p ${outPrefix}Reports
mkdir -p ${outPrefix}CombinedReports
mkdir -p ${outPrefix}UnClass
mkdir -p ${outPrefix}Class
mkdir -p ${outPrefix}Krona
mkdir -p ${outPrefix}TaxaInterest
mkdir -p ${outPrefix}Output
mkdir -p prinseqLog
mkdir -p WorkFolder
mkdir -p FilteredData

# Getting the Sample Names
#samples=($(ls -1 $in | sed -e "s/_r.*//g" -e "s/\.fastq\.gz//g" | sort | uniq))
samples=($(ls -1 $in | sed -e "s/_r.*//g" -e "s/\.fastq.*//g" | sort | uniq))

# Preparing taxa names (if needed)
if [ $taxa != "NULL" ];then
	taxa2="$(echo $taxa | sed -e 's/,/ /g')"
	echo $taxa2
fi

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
			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq WorkFolder/${sample}_r1.fastq -fastq2 WorkFolder/${sample}_r2.fastq -out_good FilteredData/${sample} -derep 14 -min_len $len -log prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20 # Need to dedup
		else
			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq WorkFolder/${sample}_r1.fastq -fastq2 WorkFolder/${sample}_r2.fastq -out_good FilteredData/${sample} -min_len $len -log prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20 # Need to dedup
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
		       --use-names $merged --output ${outPrefix}Output/${sample}Merged.out --confidence 0.1
	fi
	
	if [ -f "$paired" ]; then # Are there PE files
		echo "Paired Found"
		kraken2 --db $db --threads $ncores --report ${outPrefix}Reports/${sample}Paired.tab \
			--classified-out ${outPrefix}Class/${sample}Paired.fastq --unclassified-out ${outPrefix}Unclass/${sample}Paired.fastq \
		       --use-names $paired --output ${outPrefix}Output/${sample}Paired.out --confidence 0.1
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


#for sample in "${samples[@]}"; # The loop to analyze the paired and unpaired files
#do
#	echo "$sample being analyzed"
#	if [ $zipped == "TRUE" ]; then
#		merged="$in/$sample.fastq.gz"
#		paired="$in${sample}_r1.fastq.gz"
#	else
#		merged="$in/$sample.fastq"
#		paired="$in${sample}_r1.fastq"
#	fi
#	#echo $paired
#	if [ -f  "$merged" ]; then # This is for the merged samples
#		echo "Merged file found"
#
#		if [ $zipped == "TRUE" ] ; then
#			zcat $merged > tmp.fastq
#		else
#			cat $merged > tmp.fastq
#		fi
#
#		# This is affected by the option s
#		if [ $dedup = "TRUE" ];
#		then
#			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq tmp.fastq -out_good MergedDedup -derep 14 -min_len $len -log prinseqLog/${sample}Merged.log # Need to dedup
#		else
#			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq tmp.fastq -out_good MergedDedup -min_len $len -log prinseqLog/${sample}Merged.log # Need to dedup
#		fi
#
#		kraken2 --db $db --threads $ncores --report ${outPrefix}Reports/${sample}Merged.tab \
#			--classified-out ${outPrefix}Class/${sample}Merged.fastq --unclassified-out ${outPrefix}Unclass/${sample}Merged.fastq \
#		       --use-names MergedDedup.fastq --output ${outPrefix}Output/${sample}Merged.out
#
#		# Now to extract Reads matching taxa of interest (if applicable)
#		if [ $taxa != "NULL" ];then
#			extract_kraken_reads.py -k ${outPrefix}Output/${sample}Merged.out -r ${outPrefix}Reports/${sample}Merged.tab \
#				-o ${outPrefix}TaxaInterest/$sample.fastq -t $taxa2 --include-children --fastq-output \
#				-s MergedDedup.fastq
#		fi
#	else
#		echo "No Merged File Found"
#	fi
#	
#	if [ -f "$paired" ]; then # Are there PE files
#		echo "Paired file found"
#
#		if [ $zipped == "TRUE" ]; then
#			zcat $paired > Paired_r1.fastq
#			zcat ${paired/r1/r2} > Paired_r2.fastq
#		else
#			cat $paired > Paired_r1.fastq
#			cat ${paired/r1/r2} > Paired_r2.fastq
#
#		fi
#
#		# This is affected by the option s
#		if [ $dedup = "TRUE" ];then
#			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq Paired_r1.fastq -fastq2 Paired_r2.fastq -out_good PairedDedup -derep 14 -min_len $len -log prinseqLog/${sample}Paired.log# Need to dedup
#		else
#			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq Paired_r1.fastq -fastq2 Paired_r2.fastq -out_good PairedDedup -min_len $len -log prinseqLog/${sample}Paired.log# Need to dedup
#		fi
#
#		kraken2 --db $db --threads $ncores --report ${outPrefix}Reports/${sample}Paired.tab \
#			--classified-out ${outPrefix}Class/${sample}Paired.fastq --unclassified-out ${outPrefix}Unclass/${sample}Paired.fastq \
#		       --use-names PairedDedup_1.fastq --output ${outPrefix}Output/${sample}Paired.out
#
#		# Now to extract Reads matching taxa of interest (if applicable)
#		if [ $taxa != "NULL" ];then
#			extract_kraken_reads.py -k ${outPrefix}Output/${sample}Paired.out -r ${outPrefix}Reports/${sample}Paired.tab \
#				-o ${outPrefix}TaxaInterest/${sample}_r1.fastq -o2 ${outPrefix}TaxaInterest/${sample}_r2.fastq \
#			       	-t $taxa2 --include-children --fastq-output \
#				-s PairedDedup_1.fastq -s2 PairedDedup_2.fastq
#		fi
#	else
#		echo "No Paired File Found"
#	fi
#
#	combine_kreports.py -r ${outPrefix}Reports/$sample* --only-combined --no-headers -o ${outPrefix}CombinedReports/$sample.tab # Combining the reports
#	kreport2krona.py -r ${outPrefix}CombinedReports/$sample.tab -o ${outPrefix}Krona/$sample.txt # Preparing for a Krona Plot
#
#	rm -f Paired*fastq Merged*fastq tmp*fastq # Some Quick Cleanup of the data
#
#done
#
#ktImportText -o ${outPrefix}Krona.html ${outPrefix}Krona/* # making the KronaPlot

# The Kraken Runs
#for file in $(ls -1 $in);
#do
#	sampleWith=$(basename $file)
#	sample=${sampleWith%%.*}
#	kraken2 --db ~/Scratch/25merKraken --threads $ncores --report ${outPrefix}Reports/$sample.tab \
#		--classified-out ${outPrefix}Class/$sampleWith --unclassified-out ${outPrefix}UnClass/$sampleWith \
#	       	--use-names $file --output /dev/null
#	kreport2krona.py -r ${outPrefix}Reports/$sample.tab -o ${outPrefix}/$sample.txt
#done
