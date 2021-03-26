#! /usr/bin/env sh
# These are the files and variables that will be needed
usage() { printf 'Ancient QC Script V0.5
	Note that if planning on assembling contigs, use <FILE> instead
        -i\tThe folder that contains the Raw sequencing files
        -r\tThe BWA index.  Will be mapping reads against it
	-k\tMinimum Fragment Length (Default: 30)
	-m\tBWA MEM mode	
	-n\tNumber of CPU Threads to be used
	-q\tMinimum QUality (Default: 30)
	-l\tLog File Name (Default: $date)
	-s\tString Deduplication (Default: Mapping)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Mapping settings for $(date):
	Log File:\t${log}
	Input folder:\t${raw}
	Reference:\t${ref}
	CPU Threads:\t${ncores}
	BWA Mode:\t${bwa}
	Deduplication=\t${dedupMethod}
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
bwa="aln"
dedupMethod="Mapping"

while getopts "i:r:k:n:q:l:hms" arg; do
        case $arg in
                i)
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
		m)
			bwa="mem"
			;;
		s)
			dedupMethod="String"
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
mkdir -p trimmed_out
mkdir -p merge_logs
mkdir -p failed_out
echo "Trimming and Merging Reads"
njobs=$(echo "scale=0;var1=$ncores/16;var1"|bc) # Will round down!!!

#total=$(find $raw/* | wc -l)
#count=2
#
#for f in $raw/*R1*; do
#        sample=$(basename $f .fastq.gz)
#        sample="${sample%_R1_001}"
#        R2=$(echo "${sample}_R2_001.fastq.gz")
#        dir=$(dirname $f)
#
#        fastp -i $f -I $dir/$R2 --merge \
#        --merged_out trimmed_out/${sample}_merged.fastq.gz \
#	--out1 trimmed_out/${sample}_r1.fastq.gz \
#	--out2 trimmed_out/${sample}_r2.fastq.gz \
#        --adapter_fasta /usr/local/trimmomatic/adapters/TruSeq3-PE-2.fa --correction \
#        --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
#        --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
#        --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
#        --overlap_diff_limit 100 --n_base_limit 0\
#        --overlap_len_require 15 --length_required $len\
#        --html merge_logs/${sample}.html \
#        --json merge_logs/${sample}.json -R $sample --thread 16 -R $sample \
#        --unpaired1 trimmed_out/${sample}_u1.fastq.gz \
#        --unpaired2 trimmed_out/${sample}_u2.fastq.gz\
#        --failed_out failed_out/${sample}_failed.fastq.gz;
#
#        echo "$sample completed ($count/$total)"
#        count=$(echo "$count + 2" | bc)
#done
#
######################
######Pooling Lanes###
######################

mkdir -p PooledLanes
echo "Pooling Lanes"

### Getting the unique sample names -- Removing any information about lanes
samples=$(find trimmed_out/ -maxdepth 1 -type f | sed -e "s/trimmed_out\///g" -e "s/.fq.gz//g" -e "s/_L00[1-2].*//g" -e "s/_r[1-2]//g" | sort | uniq)
##

for sample in $samples; do
                echo "$sample has two lanes"
                cat trimmed_out/${sample}_L00*_merged.fastq.gz > PooledLanes/${sample}.fastq.gz
                cat trimmed_out/${sample}_L00*_r1.fastq.gz > PooledLanes/${sample}_r1.fastq.gz
                cat trimmed_out/${sample}_L00*_r2.fastq.gz > PooledLanes/${sample}_r2.fastq.gz
done

####################################
### Running String Deduplication ###
####################################
if [ $dedupMethod = "String" ]; then
	echo "Deduplication"
	mkdir -p PrinseqLog

	mv PooledLanes PooledLanesDups # Want to keep the non deduplicated section in case there's stuff of interest there
	mkdir -p PooledLanes
	samples=($(ls -1 PooledLanesDups | sed -e "s/_r.*//g" -e "s/\.fastq.*//g" | sort | uniq)) # Getting the list of samples

	for sample in "${samples[@]}"; # Running the loop to String Deduplicate
	do
		echo "$sample"
		merged="PooledLanesDups/$sample.fastq.gz"
		paired="PooledLanesDups/${sample}_r1.fastq.gz"
		
		if [ -f "$merged" ]; then
			zcat $merged > tmp.fastq
			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq tmp.fastq -out_good MergedDedup -derep 14 -log PrinseqLog/$sample.log
			gzip -c MergedDedup.fastq > PooledLanes/$sample.fastq.gz
		fi

		if [ -f "$paired" ]; then
			zcat $paired > Paired_r1.fastq
			zcat ${paired/r1/r2} > Paired_r2.fastq

			perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq Paired_r1.fastq -fastq2 Paired_r2.fastq -out_good PairedDedup -derep 14 -log PrinseqLog/${sample}Paired.log
			gzip -c PairedDedup_1.fastq > PooledLanes/${sample}_r1.fastq.gz
			gzip -c PairedDedup_2.fastq > PooledLanes/${sample}_r2.fastq.gz
		fi
	done
	rm -f Paired*fastq Merged*fastq tmp*fastq # Some Quick Cleanup of the data
fi

#########################################
#### String Deduplication is in Order ###
#########################################
#mkdir -p StringDeduplication
#
#for file in PooledLanes/*r1*; do
#	sample=$(basename $file .fastq.gz)
#	sample="${sample%_r1}"
#
#	zcat $file ${file/_r1/} > tmp.fastq
#	perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq tmp.fastq -out_good StringDeduplication/$sample.fastq -derep 14 -min_len $len
#done
#
#################################
### Filtering out the Bat DNA ###
#################################
if [ $bwa = "aln" ];
then
	~/Scripts/V2Folder/QC/BWAalnPairedMapping.sh -r $ref -i PooledLanes -o BWAAlnMapping -k $len -q $qual -l $logBWA -n $ncores
else
	~/Scripts/V2Folder/QC/BWAmemPairedMapping.sh -r $ref -i PooledLanes -o BWAMemMapping -k $len -q $qual -l $logBWA -n $ncores
fi

## Human Mapping
#echo "Filtering against $(basename $ref)"
#mkdir -p UnmappedReads 
#mkdir -p MappedReads
#mkdir -p DeduplicatedMappings
#total=$(find PooledLanes/*r1* | wc -l)
#count=1
#
## Human Merged
#for file in StringDeduplication/*fastq; do
#	# Getting set up for the two mappings needed here
#	sample=$(basename $file .fastq)
#	echo "$sample"
#
#	# Performing the mapping for the Merged sequences
#	bwa aln -o 2 -n 0.01 -l 16500 $ref $file -t $ncores 2>/dev/null|\
#		bwa samse $ref /dev/stdin $file 2>/dev/null| samtools sort - |\
#		samtools view -h -m $len -q 30 -U tmp.sam > Merged.sam
#	samtools fastq tmp.sam > tmp.fastq
#	~/Applications/bbmap/reformat.sh in=tmp.fastq out=Mtmp.fastq minlength=$len overwrite=true 
#
#	# Getting the Forward Reads
#	#bwa aln -o 2 -n 0.01 -l 16500 $ref $file -t $ncores 2>/dev/null | bwa samse $ref /dev/stdin $file 2>/dev/null | samtools sort - |\
#	#	samtools view -h -m $len -q 30 -U tmp.sam > Forward.sam
#	#samtools fastq tmp.sam > tmp.fastq
#	#~/Applications/bbmap/reformat.sh in=tmp.fastq out=Ftmp.fastq minlength=$len overwrite=true
#
#	# Making the BAM files
##	samtools merge -f -O BAM MappedReads/$sample.bam Merged.sam Forward.sam
##	cat Mtmp.fastq Ftmp.fastq > UnmappedReads/${sample}.fastq
#
#	samtools view -b -h Merged.sam > MappedReads/$sample.bam
#	mv Mtmp.fastq UnmappedReads/${sample}.fastq
#
#	# QoL progress
#	printf "$sample has been filtered ($count/$total)\n"
#	count=$(echo "$count + 1" | bc)
#	rm -rf Merged.sam Forward.sam *tmp*
#done
#
## Deduplicating Reads
##parallel -j $ncores --bar '/usr/local/biohazard/bin/bam-rmdup -c -o DeduplicatedMappings/{/} {} 2> /dev/null' ::: MappedReads/*bam # Removes Duplicates
##parallel -j $ncores --bar 'samtools index {} > {}.bai' ::: DeduplicatedMappings/*bam # Makes the indices
