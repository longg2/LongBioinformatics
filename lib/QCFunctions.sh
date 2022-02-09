#! /usr/bin/env bash

########################
### Ancient Trimming ###
########################
AncientTrimming(){ # This uses a combination leeHom and AdapterRemoval

	if [ "$r2" == "NA" ]; then #Must be single ended
		leeHomMulti --ancientdna --log ${out}leeHomLogs/${sample}.log -t $ncores -fq1 $r1 -fqo ${out}Trimmed/${sample} 2> ${out}Logs/${sample}.log # Don't like this as it is hardcoded....  Will it improve my issue though?
		#leeHomMulti --ancientdna -f /opt/local/trimmomatic/adapters/TruSeq3-SE.fa --log ${out}leeHomLogs/${sample}.log -t $ncores -fq1 $r1 -fqo ${out}Trimmed/${sample} 2> /dev/null

	else # It's paired

		# First step is to identify the Adapters
		/opt/local/AdapterRemoval/AdapterRemoval --identify-adapters --file1 $r1 --file2 $r2 > tmp.out 2> /dev/null
		ada1=$(grep "adapter1" tmp.out | sed -e "s/.* //g" -e "s/ .*$//g")
		ada2=$(grep "adapter2" tmp.out | sed -e "s/.* //g" -e "s/ .*$//g")
	
		# Next we want to perform the trimming
		leeHomMulti --ancientdna -f $ada1 -s $ada2 --log ${out}leeHomLogs/${sample}.log -t $ncores -fq1 $r1 -fq2 $r2 -fqo ${out}Trimmed/${sample} 2> ${out}Logs/${sample}.log
	fi

}
AncientTrimmingAssembly(){ # This uses a combination leeHom and AdapterRemoval
	local r1=$1
	local r2=$2
	local sample=$3
	local out=$4

	if [ "$r2" == "NA" ]; then # If not a paired sample...
		echo "$sample is a SE sample"
        	fastp -i $r1 \
		--out1 ${out}Trimmed/${sample}_r1.fastq.gz \
        	--adapter_fasta /usr/local-centos6/trimmomatic/adapters/TruSeq3-PE-2.fa --correction \
        	--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
        	--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
        	--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
        	--n_base_limit 0\
        	--length_required $len\
        	--html ${out}FastpLogs/${sample}.html \
        	--json ${out}FastpLogs/${sample}.json -R $sample --thread 16 -R $sample \
        	--failed_out ${out}FailedQC/${sample}_failed.fastq.gz;
	else
		echo "$sample is a PE sample"

        	fastp -i $r1 -I $r2 \
		--out1 ${out}Trimmed/${sample}_r1.fastq.gz \
		--out2 ${out}Trimmed/${sample}_r2.fastq.gz \
        	--detect_adapter_for_pe --correction \
        	--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
        	--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
        	--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
        	--n_base_limit 0 --length_required $len\
        	--html ${out}FastpLogs/${sample}.html \
        	--json ${out}FastpLogs/${sample}.json -R $sample --thread 16 -R $sample \
        	--unpaired1 ${out}Trimmed/${sample}_u1.fastq.gz \
        	--unpaired2 ${out}Trimmed/${sample}_u2.fastq.gz\
        	--failed_out ${out}FailedQC/${sample}_failed.fastq.gz;
	fi
}
FastpWrapperAncient(){ # Convenient Wrapper for parallelization
	FileIdentificationInFunction $1 $folder
	FileExtraction $folder
	AncientTrimmingAssembly
}

#######################
### Modern Trimming ###
#######################
Trimming(){ # Performing the trimming
	# Trimming
	local r1=$1
	local r2=$2
	local sample=$3
	local out=$4
	
	# Need to control for differences between info2020 and the rest
	if [ "$HOSTNAME" == "info2020" ]; then
		local localFolder="local-centos6"
	else
		local localFolder="local"

	fi

	if [ "$r2" == "NA" ]; then # If not a paired sample...
		#echo "$sample is a SE sample"
        	fastp -i $r1 \
		--out1 ${out}Trimmed/${sample}_r1.fastq.gz \
        	--adapter_fasta /usr/${localFolder}/trimmomatic/adapters/TruSeq3-PE-2.fa --correction \
        	--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
        	--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
        	--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
        	--n_base_limit 0\
        	--length_required $len\
        	--html ${out}FastpLogs/${sample}.html \
        	--json ${out}FastpLogs/${sample}.json -R $sample --thread 16 -R $sample \
        	--failed_out ${out}FailedQC/${sample}_failed.fastq.gz;
	else
		#echo "$sample is a PE sample"
        	fastp -i $r1 -I $r2 --merge \
        	--merged_out ${out}Trimmed/${sample}_merged.fastq.gz \
		--out1 ${out}Trimmed/${sample}_r1.fastq.gz \
		--out2 ${out}Trimmed/${sample}_r2.fastq.gz \
        	--adapter_fasta /usr/${localFolder}/trimmomatic/adapters/TruSeq3-PE-2.fa --correction \
        	--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
        	--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
        	--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
        	--overlap_diff_limit 100 --n_base_limit 0\
        	--overlap_len_require 15 --length_required $len\
        	--html ${out}FastpLogs/${sample}.html \
        	--json ${out}FastpLogs/${sample}.json -R $sample --thread 16 -R $sample \
        	--unpaired1 ${out}Trimmed/${sample}_u1.fastq.gz \
        	--unpaired2 ${out}Trimmed/${sample}_u2.fastq.gz\
        	--failed_out ${out}FailedQC/${sample}_failed.fastq.gz;
	fi
}
FastpWrapper(){ # Convenient Wrapper for parallelization
	FileIdentificationInFunction $1 $folder
	FileExtraction $folder
	Trimming
}

# String Deduplication Commands.  Parallelization may or may not work
StringDeduplication(){ 
	local sample=$1
	local folderFunction=$2
	local len=$3
	FileIdentificationInFunction $sample $folderFunction
	FileExtractionInFunction $folderFunction

	# Now for the actual deduplication
	if [ "$merged" != "NA" ]; then # If I found a merged file
		if [ "$Dedup" == "TRUE" ]; then
			prinseq -fastq $merged -out_bad null -out_good stdout -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20 | gzip > ${out}StringDedup/${sample}_Merged.fastq.gz
		else
			prinseq -fastq $merged -out_bad null -out_good stdout -min_len $len -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20 | gzip > ${out}StringDedup/${sample}_Merged.fastq.gz
		fi
	fi

	if [[ $r1 != "NA" && $r2 == "NA" ]]; then # If dealing with a single end library
		echo "Single"
		if [ "$Dedup" == "TRUE" ]; then
			prinseq -fastq $r1 -out_bad null -out_good stdout -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 | gzip > ${out}StringDedup/${sample}_r1.fastq.gz
		else
			prinseq -fastq $r1 -out_bad null -out_good stdout -min_len $len -log ${out}prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 | gzip > ${out}StringDedup/${sample}_r1.fastq.gz
		fi
	fi 

	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
		if [ "$Dedup" == "TRUE" ]; then
			prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20
		else
			prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20

		fi

	# Because of how prinseq is coded, I'll need to compress separately	
		gzip -c TMP/${sample}_1.fastq > ${out}StringDedup/${sample}_r1.fastq.gz
		gzip -c TMP/${sample}_2.fastq > ${out}StringDedup/${sample}_r2.fastq.gz

	fi
}
StringDeduplicationParallel(){
	local sample=$1
	local folderFunction=$2
	local len=$3
	FileIdentificationInFunction $sample $folderFunction
	FileExtractionInFunction $folderFunction

	# Now for the actual deduplication
	if [ "$merged" != "NA" ]; then # If I found a merged file
		if [ "$Dedup" == "TRUE" ]; then
			prinseq -fastq $merged -out_bad null -out_good TMP/${sample} -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20 
		else
			prinseq -fastq $merged -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20
		fi
	fi

	if [[ $r1 != "NA" && $r2 == "NA" ]]; then # If dealing with a single end library
		if [ "$Dedup" == "TRUE" ]; then
			prinseq -fastq $r1 -out_bad null -out_good TMP/${sample} -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 
		else
			prinseq -fastq $r1 -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 
		fi
	fi 

	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
		if [ "$Dedup" == "TRUE" ]; then
			prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -derep 14 -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20
		else
			prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20

		fi

	# Because of how prinseq is coded, I'll need to compress separately	
		gzip -c TMP/${sample}.fastq > ${out}StringDedup/${sample}.fastq.gz
		gzip -c TMP/${sample}_1.fastq > ${out}StringDedup/${sample}_r1.fastq.gz
		gzip -c TMP/${sample}_2.fastq > ${out}StringDedup/${sample}_r2.fastq.gz

	fi
}

AssemblyPreparation(){ 
	mkdir -p ${out}prinseqLog
	mkdir -p ${out}LengthFiltered
	# Now for the actual deduplication
	if [ "$merged" != "NA" ]; then # If I found a merged file
		prinseq -fastq $merged -out_bad null -out_good stdout -min_len $len -trim_left $chomp -trim_right $chomp -log ${out}prinseqLog/${sample}Merged.log -lc_method dust -lc_threshold 20 | gzip > ${out}LengthFiltered/${sample}_Merged.fastq.gz
	fi

	if [[ $r1 != "NA" && $r2 == "NA" ]]; then # If dealing with a single end library
		prinseq -fastq $r1 -out_bad null -out_good stdout -min_len $len -log ${out}/prinseqLog/${sample}Single.log -lc_method dust -lc_threshold 20 | gzip > ${out}LengthFiltered/${sample}_r1.fastq.gz
	fi 

	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
		prinseq -fastq $r1 -fastq2 $r2 -out_bad null -out_good TMP/${sample} -min_len $len -log ${out}prinseqLog/${sample}Paired.log -lc_method dust -lc_threshold 20

	fi

	# Because of how prinseq is coded, I'll need to compress separately	
		gzip -c TMP/${sample}_1.fastq > ${out}LengthFiltered/${sample}_r1.fastq.gz
		gzip -c TMP/${sample}_2.fastq > ${out}LengthFiltered/${sample}_r2.fastq.gz
}
AdapterRemovalHeavy(){
	local sample=$1
	local folderFunction=$2
	local out=$3

	# Make the required folders
	mkdir -p ${out}AdaptersFiltered
	mkdir -p ${out}AdaptersSingle

	FileIdentificationInFunction $sample $folderFunction
	FileExtractionInFunction $folderFunction
	if [ "$merged" != "NA" ]; then # If I found a merged file
		zcat $merged | paste - - - - | agrep -v -4 GGAAGAGCGTCGTGTAGGGAAAGAG |\
			agrep -v -4 ATCGGAAGAGCACACGTCTGAACTCCAG |\
			agrep -v -4 ATCGGAAGAGCGTCGTGTAGGGAAAGAG |\
			agrep -v -4 ATCGGAAGAGCACACGTCTGAACTCCAG | tr "\t" "\n" | gzip -c > ${out}AdaptersFiltered/$sample.fq.gz
	fi

	if [[ $r1 != "NA" && $r2 == "NA" ]]; then # If dealing with a single end library
		zcat $r1 | paste - - - - | agrep -v -4 GGAAGAGCGTCGTGTAGGGAAAGAG |\
			agrep -v -4 ATCGGAAGAGCACACGTCTGAACTCCAG |\
			agrep -v -4 ATCGGAAGAGCGTCGTGTAGGGAAAGAG |\
			agrep -v -4 ATCGGAAGAGCACACGTCTGAACTCCAG | tr "\t" "\n" | gzip -c > ${out}AdaptersFiltered/${sample}_r1.fq.gz
	fi 

	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
		zcat $r1 | paste - - - - | agrep -v -4 GGAAGAGCGTCGTGTAGGGAAAGAG |\
			agrep -v -4 ATCGGAAGAGCACACGTCTGAACTCCAG |\
			agrep -v -4 ATCGGAAGAGCGTCGTGTAGGGAAAGAG |\
			agrep -v -4 ATCGGAAGAGCACACGTCTGAACTCCAG | tr "\t" "\n" | gzip -c > TMP/${sample}_r1.fq.gz

		zcat $r2 | paste - - - - | agrep -v -4 GGAAGAGCGTCGTGTAGGGAAAGAG |\
			agrep -v -4 ATCGGAAGAGCACACGTCTGAACTCCAG |\
			agrep -v -4 ATCGGAAGAGCGTCGTGTAGGGAAAGAG |\
			agrep -v -4 ATCGGAAGAGCACACGTCTGAACTCCAG | tr "\t" "\n" | gzip -c > TMP/${sample}_r2.fq.gz

		# In this case, we need to remove all the singletons....
		/home/sam/Applications/bbmap/repair.sh in1=TMP/${sample}_r1.fq.gz in2=TMP/${sample}_r2.fq.gz out1=${out}AdaptersFiltered/${sample}_r1.fq.gz out2=${out}AdaptersFiltered/${sample}_r2.fq.gz \
			outs=${out}AdaptersSingle/${sample}_single.fq.gz showspeed=f repair=t 2>/dev/null

	fi
}
