#! /usr/bin/env bash

AncientTrimming(){ # This uses a combination leeHom and AdapterRemoval

	if [ "$r2" == "NA" ]; then #Must be single ended
		leeHomMulti --ancientdna --log ${out}leeHomLogs/${sample}.log -t $ncores -fq1 $r1 -fqo ${out}Trimmed/${sample} 2> /dev/null # Don't like this as it is hardcoded....  Will it improve my issue though?
		#leeHomMulti --ancientdna -f /opt/local/trimmomatic/adapters/TruSeq3-SE.fa --log ${out}leeHomLogs/${sample}.log -t $ncores -fq1 $r1 -fqo ${out}Trimmed/${sample} 2> /dev/null

	else # It's paired

		# First step is to identify the Adapters
		/opt/local/AdapterRemoval/AdapterRemoval --identify-adapters --file1 $r1 --file2 $r2 > tmp.out 2> /dev/null
		ada1=$(grep "adapter1" tmp.out | sed -e "s/.* //g" -e "s/ .*$//g")
		ada2=$(grep "adapter2" tmp.out | sed -e "s/.* //g" -e "s/ .*$//g")
	
		# Next we want to perform the trimming
		leeHomMulti --ancientdna -f $ada1 -s $ada2 --log ${out}leeHomLogs/${sample}.log -t $ncores -fq1 $r1 -fq2 $r2 -fqo ${out}Trimmed/${sample} 2> /dev/null
	fi

}
ProgressBar() { # From github.com/fearside/ProgressBar
	# Process data
		let _progress=(${1}*100/${2}*100)/100
		let _done=(${_progress}*4)/10
		let _left=40-$_done
	# Build progressbar string lengths
	_done=$(printf "%${_done}s")
	_left=$(printf "%${_left}s")
	printf "\rProgress : [${_done// />}${_left// /-}] ${_progress}%%"

}	

Trimming(){ # Performing the trimming
	# Trimming
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
        	fastp -i $r1 -I $r2 --merge \
        	--merged_out ${out}Trimmed/${sample}_merged.fastq.gz \
		--out1 ${out}Trimmed/${sample}_r1.fastq.gz \
		--out2 ${out}Trimmed/${sample}_r2.fastq.gz \
        	--adapter_fasta /usr/local-centos6/trimmomatic/adapters/TruSeq3-PE-2.fa --correction \
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

