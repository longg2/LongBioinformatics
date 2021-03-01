# To be done:
# 	* Complete rewrite with more of a focus on functions
# 	* Remove reliance of seqkit
# 	* Progress bar for Blast (will require a trap and while loop)
# 	* Integration (with option for ignoring) the LCA script for blast

usage() { printf "BlastN/P Wrapper Script V0.6
	Outputs tab deliminated BlastN/P report file in the form of std staxid.  Taxa counts
	for each step are also outputted.  Need seqkit and taxonkit installed.
	String Deduplication occurs.
	IMPORTANT: Assumes that fasta/fastq files aren't compressed.  Will result in fireworks
	if compressed.  To be fixed.
        -i\tInput Folder
	-o\tOutput Folder (Default = BlastOut)
	-b\tblastn, blastp, or blastx? (Default = blastn)
	-d\tBlast database
	-e\tMaximum Evalue (Default = 1e-5)
	-k\tRead Length (Default = 30)
	-l\tLog File
	-p\tPercent Identity (Default = 90)
	-s\tReads to Subsample (Default = 0)
        -n\tNumber of CPU Threads to be used
        -h\tShow this help message and exit\n" 1>&2; exit 1; }

blastCMD() {
	if [ "$blast" == "blastn" ]; then
		blastn -db $db -query tmp.fasta -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores -perc_identity $Pident
	elif [ "$blast" == "blastp" ]; then
		blastp -db $db -query tmp.fasta -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores
	elif [ "$blast" == "blastx" ]; then
		blastx -db $db -query_gencode 11 -query tmp.fasta -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores 
	else
		echo "Please choose either blastn or blastp"
		return 1
	fi

	return 0
}

log() {	printf "Blast settings for $(date):
	Log File:\t${log}
	Input folder:\t${in}
	Output folder:\t${out}
	Blast Type:\t${blast}
	Blast Database:\t${db}
	Computer:\t${HOSTNAME}
	-------------------------------------
	E-value:\t${eval}
	Percent Identity:\t${Pident}
	Subsampled Reads:\t${subsample}
	Min Read Length:\t${len}	
	CPU Threads:\t${ncores}\n"; exit 0;
}

# Preset Variables
eval="1e-5"
Pident="90"
out="BlastOut"
blast="blastn"
subsample=0
len=30
log="$(date +'%Y%m%d').log"

while getopts "i:e:d:n:o:b:p:l:s:h" arg; do
        case $arg in
                i)
                        in=${OPTARG}
                        ;;
		e)
			eval=${OPTARG}
			;;
                d)
                        db=${OPTARG}
                        ;;
                n)
                        ncores=${OPTARG}
                        ;;
		o)
			out=${OPTARG}
			;;
		p)
			Pident=${OPTARG}
			;;
		b)
			blast=${OPTARG}
			;;
		s)
			subsample=${OPTARG}
			;;
		l)
			log=${OPTARG}
			;;
		k)
                        len=${OPTARG}
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Printing the Logfile
log | tee $log

# Making the Directories
mkdir -p ${out}
mkdir -p prinseqLog
if [ "$blast" == "blastn" ]; then
	mkdir -p ${out}/${out}Script
fi

for file in $in/*; do
	sample=$(basename $file)
	sample=${sample%%.*}
	echo $sample

	# First, we need to gunzip the files if applicable
#	if [[ $file =~ \.t?gz$ ]];
#	then
#		echo "$sample was gzipped"
#		gunzip -c $file > $file
#	else
#		cp $in/$sample.fastq merge.fastq
#		cp $file r1.fastq
#		cp $r2 r2.fastq
#	fi

	# Need to check the file format.  Easiest way will be to look at the first character 
	firstChar=$(cut -c1 $file | head -n 1)
	if [ "$firstChar" == ">" ];
	then
		cp $file tmp1.fasta
		perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fasta tmp1.fasta -out_good tmp -min_len $len -derep 14 -log prinseqLog/${sample}Merged.log # Remove Short Reads
	elif [ "$firstChar" == "@" ];
	then
		seqkit fq2fa -o tmp1.fasta $file
		perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fasta tmp1.fasta -out_good tmp -min_len $len -derep 14 -log prinseqLog/${sample}Merged.log # Remove Short Reads
	else
		echo "$file is not in fasta or fastq format. Exiting"
		exit 1
	fi

	# Now checking if Subsampling is being done on the file
	if [[ $subsample != 0 ]]; then
		# Saving the Subsampled Amounts
		mkdir -p Subsample${subsample}

		mv tmp.fasta tmpFull.fasta
		seqkit sample tmpFull.fasta -n $subsample > tmp.fasta

		cp tmp.fasta Subsample${subsample}/$sample.fasta

	fi
	
	blastCMD

	if [ $? -ne 0 ]; then
		printf "Either Blast panicked or you did not provide the correct blast command... Stopping the script on file $sample\n"
		exit 1
	fi

	# Now to run taxonkit
	if [ "$blast" == "blastn" ]; then
		cat ${out}/$sample.tab | cut -f 1,13 |\
			sort --compress-program gzip | uniq -c | sed -e "s/^ *//g" -e "s/ /\t/g" | sort -k3 > tmp_Counts.tab # Getting the Counts

		cut -f 3 tmp_Counts.tab | uniq | taxonkit lineage | taxonkit reformat -t -R "NA"| cut -f 1,4 > taxonKitBlast.out
		join -1 3 -2 1 tmp_Counts.tab taxonKitBlast.out | cut -f 2- -d " " | sed -e "s/ /\t/g" -e "s/;/\t/g" | sort -k 2 > ${out}/${out}Script/$sample.tab
	fi
done

rm -f tmp.fasta
rm -f tmpFull.fasta
rm -f tmp_Counts.tab
rm -f taxonKitBlast.out
