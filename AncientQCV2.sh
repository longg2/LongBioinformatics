#! /usr/bin/env bash

DeduplicateArray(){ # Deduplicating the sample names
	local array=("$@") # This feels weird, but, it'll allow you pass an array to it
	local sampleNames=$(for name in ${array[@]}; do
		tmp="$(echo $name |sed -e 's/_r1.*//I' -e 's/_r2.*//I' -e 's/_merged.*//I' -e 's/\.fa.*//I')";
		echo $tmp;
       	done) # Getting only the sample names
	samples=( $(echo ${sampleNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') )

}

FileIdentification(){ # Extract Files from an array using results from another array
	local sample=$1 # Basename of the file
	local arrayFiles=$files

	# Finding the indices which have the same samplename
	sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))
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
	#printf '%s\n' "${sampleFiles[@]}" #> .hiddenlist.list

	# Identifying the files
	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r1"; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r1')
	        r1="$folder/$fileName"
	fi

	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r2"; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r2')
	        r2="$folder/$fileName"
	fi

	# Two cases for the merged.  Want to control for shenanigans
	if printf '%s\n' "${sampleFiles[@]}" | grep -q -i "$sample\.f.*"; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -v -q "_\.f*") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -i "$sample\.f.*")
	        merged="$folder/$fileName"
	fi

	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged"; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'merged')
	        merged="$folder/$fileName"
	fi

}

AncientTrimming(){ # This uses a combination leeHom and AdapterRemoval

	# First step is to identify the Adapters
	AdapterResults=$(/opt/local/AdapterRemoval/AdapterRemoval --identify-adapters --file $r1 --file2 $r2)
	ada1=$(echo $AdapterResults | grep "adapter1" | sed -e "s/.* //g" -e "s/ .*$//g" )
	ada2=$(echo $AdapterResults | grep "adapter2" | sed -e "s/.* //g" -e "s/ .*$//g" )

	# Next we want to perform the trimming
	leeHomMulti --ancientdna -f $ada1 -s $ada2 --log ${out}leeHomLogs/${sample}.log -t $ncores -fq1 $r1 -fq2 $r2 -fqo ${out}Trimmed/${sample}

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

usage() { printf 'Ancient QC Script V1
	Very simple.  Want this to focus primarily the workup
	prior to analysis.
        -i\tThe folder that contains the Raw sequencing files
        -n\tNumber of CPU Threads to be used
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Ancient Mapping settings for $(date):
	Log File:\t${log}
	Input folder:\t${raw}
	CPU Threads:\t${ncores}\n"; exit 0;
}

##################################
#Default Values
out="BWAMappingScript"
ncores=8
log="$(date +'%Y%m%d').log"

while getopts "i:n:h" arg; do
        case $arg in
                i)
                        raw=${OPTARG}
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                n)
                        ncores=${OPTARG}
                        #echo "Using $ncores"
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
###################################################################
# Writing the Log File
log | tee $log

mkdir -p ${out}Trimmed
mkdir -p ${out}PooledLanes
###################################################################

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

echo "Trimming the ${#samples[@]} found in ${folder}.  This may take a while"
total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in ${samples[@]}; do
	FileIdentification $sample
	FileExtraction
	AncientTrimming
	
	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done

###################
###Pooling Lanes###
###################

pooledNames=$(for name in ${samples[@]}; do # This was lifted from DeduplicateArray().
	tmp="$(echo $name |sed -e 's/_L00.*//I')";
		echo $tmp;
       	done) # Getting only the sample names
pooledNames=( $(echo ${pooledNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') ) # This was lifted from Deduplicate Array

echo "Pooling the lanes together"
parallel --bar -j $ncores "cat ${out}Trimmed/{}*merged.fastq.gz > ${out}PooledLanes/{}.fastq.gz;
	cat ${out}Trimmed/{}*r1.fastq.gz > ${out}PooledLanes/{}_r1.fastq.gz;
	cat ${out}Trimmed/{}*r2.fastq.gz > ${out}PooledLanes/{}_r2.fastq.gz;" ::: "${pooledNames[@]}" # "${samples[@]}"
