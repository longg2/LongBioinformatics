#! /usr/bin/env sh
# These are the files and variables that will be needed
usage() { printf 'BWA ALN Mapping V0.5
	Mapping using aln is more complicated than it warrants,
	so here is a useful wrapper that will hopefully make life
	easier.
	TBD: Dynamically detect whether single or paired reads are used.
	TBD: Gunzip when needed
        -i\tThe folder the fasta files
	-o\tOutput folder of the BAM files (Default: BWAMappingScript)
        -r\tThe BWA index.  Will be mapping reads against it
	-k\tMinimum Fragment Length (Default: 30)
	-q\tMinimum Mapping Quality (Default: 30)
        -n\tNumber of CPU Threads to be used
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Mapping settings for $(date):
	Log File:\t${log}
	Input folder:\t${in}
	Output folder:\t${out}
	Reference:\t${ref}
	CPU Threads:\t${ncores}
	-------------------------------------
	Min Quality:\t${qual}
	Min Length:\t${len}\n"; exit 0;
}


# Creating the functions that will do the heavy lifting
PairedalnMapping(){
	local sample=$(basename $file)
	local r2="${file/r1/r2}"
	local sample="${sample%_r*}"
	echo "$sample"
	#echo $r2
	
	# Testing if the file is gzipped
	if [[ $file =~ \.t?gz$ ]];
	then
		echo "$sample was gzipped"
		gunzip -c $in/$sample.fastq > merge.fastq
		gunzip -c $file > r1.fastq
		gunzip -c $r2 > r2.fastq
	else
		cp $in/$sample.fastq merge.fastq
		cp $file r1.fastq
		cp $r2 r2.fastq
	fi


	# Performing the mapping for the Merged Reads
	bwa aln -o 2 -n 0.01 -l 16500 $ref merge.fastq -t $ncores 2>/dev/null > tmp.sai
	bwa samse $ref tmp.sai merge.fastq 2>/dev/null |\
		samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam|\
		samtools sort -	> tmpM.bam

	samtools fastq tmp.bam > ${out}UnmappedReads/${sample}.fastq
	
	# Performing the mapping for the Paired Reads
	bwa aln -o 2 -n 0.01 -l 16500 $ref r1.fastq -t $ncores 2>/dev/null > tmp1.sai # First index
	bwa aln -o 2 -n 0.01 -l 16500 $ref r2.fastq -t $ncores 2>/dev/null > tmp2.sai # Second Index
	bwa sampe $ref tmp1.sai tmp2.sai r1.fastq r2.fastq 2>/dev/null |\
		samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\
	       	samtools sort - > tmpP.bam

	samtools fastq tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq -2 ${out}UnmappedReads/${sample}_r2.fastq -s /dev/null

	# Merging the bam files
	samtools merge -f ${out}MappedReads/$sample.bam tmpM.bam tmpP.bam

	# QoL progress
	printf "$sample has been filtered ($count/$total)\n"
	count=$(echo "$count + 1" | bc)
	rm -rf Merged.sam Forward.sam *tmp*

}

# Need to figure out if it's paired or not
PairedOrSingle(){
	local sample=$1
	local NumFiles=$(grep -c $1)
	if [ $NumFiles = 1 ]; then
		echo "$sample only has one file.  Assuming all Merged"
		MergedalnMapping
	elif [ $NumFiles = 2 ]; then
		echo "$sample has two files.  Assuming only Paired"
		PairedOnlyalnMapping
	elif [ $NumFiles = 3 ]; then
		echo "$sample"
}

#Default values
qual=30
len=30
out="BWAMappingScript"
ncores=8
log="$(date +'%Y%m%d').log"

while getopts "i:o:q:r:l:k:n:h" arg; do
        case $arg in
                i)
                        in=${OPTARG}
                        #echo "The raw sequencing files are located in $in"
                        ;;
                o)
                        out=${OPTARG}
                        #echo "$out is the output folder"
                        ;;
		q)
			qual=${OPTARG}
			#echo "Requiring a minimum mapping quality of $qual"
                        ;;
                r)
                        ref=${OPTARG}
                        #echo "We're mapping the QCed reads against $ref"
                        ;;
                l)
                        log=${OPTARG}
                        #echo "Settings are being outputted to $log"
                        ;;
                k)
                        len=${OPTARG}
                        #echo "Using a minimum length of ${len}bp"
                        ;;
                n)
                        ncores=${OPTARG}
                        #echo "Using $ncores CPU Threads"
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Writing the Log File
log | tee $log
# The actual loop

mkdir -p ${out}MappedReads
mkdir -p ${out}UnmappedReads 
mkdir -p ${out}DeduplicatedMappings
total=$(find $in/*r1* | wc -l)
count=1

for file in $in/*r1.fastq*; do
	# Getting set up for the two mappings needed here
	sample=$(basename $file)
	r2="${file/r1/r2}"
	sample="${sample%_r*}"
	echo "$sample"
	#echo $r2
	
	# Testing if the file is gzipped
	if [[ $file =~ \.t?gz$ ]];
	then
		echo "$sample was gzipped"
		gunzip -c $in/$sample.fastq > merge.fastq
		gunzip -c $file > r1.fastq
		gunzip -c $r2 > r2.fastq
	else
		cp $in/$sample.fastq merge.fastq
		cp $file r1.fastq
		cp $r2 r2.fastq
	fi


	# Performing the mapping for the Merged Reads
	#bwa aln -o 2 -n 0.01 -l 16500 $ref $in/$sample.fastq -t $ncores 2>/dev/null > tmp.sai
	#bwa samse $ref tmp.sai $in/$sample.fastq 2>/dev/null |\
	bwa aln -o 2 -n 0.01 -l 16500 $ref merge.fastq -t $ncores 2>/dev/null > tmp.sai
	bwa samse $ref tmp.sai merge.fastq 2>/dev/null |\
		samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam|\
		samtools sort -	> tmpM.bam

	samtools fastq tmp.bam > ${out}UnmappedReads/${sample}.fastq
	
	# Performing the mapping for the Paired Reads
	#bwa aln -o 2 -n 0.01 -l 16500 $ref $file -t $ncores 2>/dev/null > tmp1.sai # First index
	#bwa aln -o 2 -n 0.01 -l 16500 $ref $r2 -t $ncores 2>/dev/null > tmp2.sai # Second Index
	#bwa sampe $ref tmp1.sai tmp2.sai $file $r2 2>/dev/null |\
	bwa aln -o 2 -n 0.01 -l 16500 $ref r1.fastq -t $ncores 2>/dev/null > tmp1.sai # First index
	bwa aln -o 2 -n 0.01 -l 16500 $ref r2.fastq -t $ncores 2>/dev/null > tmp2.sai # Second Index
	bwa sampe $ref tmp1.sai tmp2.sai r1.fastq r2.fastq 2>/dev/null |\
		samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\
	       	samtools sort - > tmpP.bam

	samtools fastq tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq -2 ${out}UnmappedReads/${sample}_r2.fastq -s /dev/null
	# Merging the bam files
	samtools merge -f ${out}MappedReads/$sample.bam tmpM.bam tmpP.bam

	# QoL progress
	printf "$sample has been filtered ($count/$total)\n"
	count=$(echo "$count + 1" | bc)
	rm -rf Merged.sam Forward.sam *tmp*
done

rm -f merge.fastq r1.fastq r2.fastq

## Deduplicating the Reads
parallel -j $ncores --bar "/usr/local/biohazard/bin/bam-rmdup -c -o ${out}DeduplicatedMappings/{/} {} 2> /dev/null" ::: ${out}MappedReads/*bam # Removes Duplicates

# In case the depths are being sent elsewhere
#tar -czf ${out}Dedup.tar.gz ${out}Dedup

echo "Mapping is complete!"
