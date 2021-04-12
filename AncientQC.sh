#! /usr/bin/env sh
# These are the files and variables that will be needed
usage() { printf 'Ancient QC Script V0.5
	Note that if planning on assembling contigs, use <FILE> instead
        -s\tThe folder that contains the Raw sequencing files
        -r\tThe BWA index.  Will be mapping reads against it
	-k\tMinimum Fragment Length (Default: 30)
        -n\tNumber of CPU Threads to be used
	-q\tMinimum QUality (Default: 30)
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Mapping settings for $(date):
	Log File:\t${log}
	Input folder:\t${raw}
	Reference:\t${ref}
	CPU Threads:\t${ncores}
	-------------------------------------
	Min Quality:\t${qual}
	Min Length:\t${len}\n"; exit 0;
}

#Default Values
qual=30
len=30
out="BWAMappingScript"
ncores=8
log="$(date +'%Y%m%d').log"

while getopts "s:r:k:n:q:l:h" arg; do
        case $arg in
                s)
                        raw=${OPTARG}
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                r)
                        ref=${OPTARG}
                        #echo "We're mapping the QCed reads against $ref"
                        ;;
                k)
                        len=${OPTARG}
                        #echo "Using a minimum length of ${len}bp"
                        ;;
                n)
                        ncores=${OPTARG}
                        #echo "Using $ncores"
                        ;;
		q)
			qual=${OPTARG}
			#echo "Requiring a minimum mapping quality of $qual"
                        ;;
                l)
                        log=${OPTARG}
                        #echo "Settings are being outputted to $log"
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done
logBWA="${log/%.*}BWA.log"

# Writing the Log File
log | tee $log

####################################
###Trimming the reads with leeHom###
####################################
mkdir -p leeHomOut
mkdir -p leeHomOutLogs
echo "Trimming and Merging Reads"

total=$(find $raw/* | wc -l)
count=2

for f in $raw/*R1*; do
	sample=$(basename $f .fastq.gz)
	sample="${sample%_R1_001}"
	R2=$(echo "${sample}_R2_001.fastq.gz")
	dir=$(dirname $f)

	/opt/local/AdapterRemoval/AdapterRemoval --identify-adapters --file1 $f --file2 $dir/$R2 > tmp.out
	ada1=$(grep "adapter1" tmp.out | sed -e "s/.* //g" -e "s/ .*$//g")
	ada2=$(grep "adapter2" tmp.out | sed -e "s/.* //g" -e "s/ .*$//g")

	leeHomMulti --ancientdna -f $ada1 -s $ada2 --log leeHomOutLogs/$sample.log -t $ncores -fq1 $f -fq2 $dir/$R2 -fqo leeHomOut/$sample

	echo "$sample completed ($count/$total)"
	count=$(echo "$count + 2" | bc)
done

mkdir -p leeHomOut/Failed
mv leeHomOut/*fail* leeHomOut/Failed/

###################
###Pooling Lanes###
###################
mkdir -p PooledLanes
echo "Pooling Lanes Together"

### Getting the unique sample names -- Removing any information about lanes
samples=$(find leeHomOut/ -maxdepth 1 -type f | sed -e "s/leeHomOut\///g" -e "s/.fq.gz//g" -e "s/_L00[1-2].*//g" -e "s/_r[1-2]//g" | sort | uniq)
echo $samples
##

for sample in $samples; do
	if [ $(ls leeHomOut/${sample}* | wc -l) -eq 6 ]; then
		echo "$sample has two lanes"
		cat leeHomOut/${sample}_L001.fq.gz leeHomOut/${sample}_L002.fq.gz > PooledLanes/${sample}.fastq.gz
		cat leeHomOut/${sample}_L001_r1.fq.gz leeHomOut/${sample}_L002_r1.fq.gz > PooledLanes/${sample}_r1.fastq.gz
		cat leeHomOut/${sample}_L001_r2.fq.gz leeHomOut/${sample}_L002_r2.fq.gz > PooledLanes/${sample}_r2.fastq.gz
	else
		echo "$sample has one lane"
		cat leeHomOut/${sample}_L00*.fq.gz > PooledLanes/${sample}.fastq.gz
		cat leeHomOut/${sample}_L00*_r1.fq.gz > PooledLanes/${sample}_r1.fastq.gz
		cat leeHomOut/${sample}_L00*_r2.fq.gz > PooledLanes/${sample}_r2.fastq.gz
		#cp -f leeHomOut/${sample}_L00*.fq.gz PooledLanes/${sample}.fastq.gz # Weird bug.... This line didn't work
		#cp -f leeHomOut/${sample}_L00*_r1.fq.gz PooledLanes/${sample}_r1.fastq.gz
		#cp -f leeHomOut/${sample}_L00*_r2.fq.gz PooledLanes/${sample}_r2.fastq.gz
	fi

done

#####################################
### Mapping against the Reference ###
#####################################

#~/Scripts/V2Folder/QC/BWAalnPairedMapping.sh -r $ref -i PooledLanes -o BWAMapping -k $len -q $qual -l $logBWA -n $ncores 

# Human Mapping
#echo "Filtering against $(basename $ref)"
#mkdir -p UnmappedReads 
#mkdir -p MappedReads
#mkdir -p DeduplicatedMappings
##refFilt="~/Scratch/Covid/BatLib/BatGenomes.fna"
#total=$(find PooledLanes/*r1* | wc -l)
#count=1
#
## Human Merged
#for file in PooledLanes/*r1*; do
#	# Getting set up for the two mappings needed here
#	sample=$(basename $file .fastq.gz)
#	sample="${sample%_r1}"
#	echo "$sample"
#
#	# Performing the mapping for the Merged sequences
#	bwa aln -o 2 -n 0.01 -l 16500 $ref PooledLanes/$sample.fastq.gz -t $ncores 2>/dev/null|\
#		bwa samse $ref /dev/stdin PooledLanes/$sample.fastq.gz 2>/dev/null| samtools sort - |\
#		samtools view -h -m $len -U tmp.sam > Merged.sam
#	samtools fastq tmp.sam > tmp.fastq
#	~/Applications/bbmap/reformat.sh in=tmp.fastq out=Mtmp.fastq minlength=$len overwrite=true 
#
#	# Getting the Forward Reads
#	bwa aln -o 2 -n 0.01 -l 16500 $ref $file -t $ncores 2>/dev/null |\  # First index
#		bwa samse $ref /dev/stdin $file 2>/dev/null | samtools sort - |\
#		samtools view -h -m $len -U tmp.sam > Forward.sam
#	samtools fastq tmp.sam > tmp.fastq
#	~/Applications/bbmap/reformat.sh in=tmp.fastq out=Ftmp.fastq minlength=$len overwrite=true
#
#	# Making the BAM files
#	samtools merge -f -O BAM MappedReads/$sample.bam Merged.sam Forward.sam
#	cat Mtmp.fasta Ftmp.fasta > UnmappedReads/${sample}.fasta
#
#	# QoL progress
#	printf "$sample has been filtered ($count/$total)\n"
#	count=$(echo "$count + 1" | bc)
#	rm -rf Merged.sam Forward.sam *tmp*
#done
#
## Deduplicating Reads
#parallel -j $ncores --bar '/usr/local/biohazard/bin/bam-rmdup -c -o DeduplicatedMappings/{/} {} 2> /dev/null' ::: MappedReads/*bam # Removes Duplicates
#parallel -j $ncores --bar 'samtools index {} > {}.bai' ::: DeduplicatedMappings/*bam # Makes the indices
