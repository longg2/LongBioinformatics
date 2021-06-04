#! /usr/bin/env bash
# These are the functions that actually do the work.  Mean to make parallelization easier to do
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
	printf '%s\n' "${sampleFiles[@]}" > .hiddenlist.list

	# Identifying the files
	if grep -P -i -q "_r1" .hiddenlist.list; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i '_r1')
	        r1="$folder/$fileName"
	fi

	if grep -P -i -q "_r2" .hiddenlist.list; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i '_r2')
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
	samples=( $(echo ${sampleNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') )

}

SPAdesFunction(){
	spades --isolate -1 $r1 -2 $r2 --merged $merged -t $ncores -o ${out}/$sample
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
# These are the files and variables that will be needed
usage() { printf 'SPAdes Modern Assembly Script V0.5
	Absolute minimal script to get SPAdes running
        -i\tThe folder that contains the Raw sequencing files
	-o\tThe output prefix (Default: ModernIsolate)
	-n\tNumber of CPU Threads to be used (Default: 16)
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "SPAdes Isolate Settings (Modern) for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

out="ModernIsolate"
ncores=16
log="$(date +'%Y%m%d').log"

while getopts "i:n:o:l:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                o)
                        out=${OPTARG}
			export $out
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                n)
                        declare -i ncores=${OPTARG}
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

if [ -z ${folder} ]; then # Testing if I have an input folder.  That's all I need to get running
	printf "The input was not set"
	usage
	exit 1
fi

# Writing the Log File
log | tee $log

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs the variable samples

mkdir -p ${out}
mkdir -p ${out}/SPAdesLogs
#########################
###Making the Assembly###
#########################

total=${#samples[@]}
count=0
ProgressBar $count $total
for sample in ${samples[@]}; do
	FileIdentification $sample
	FileExtraction
	SPAdesFunction > ${out}/SPAdesLogs/$sample.log

	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
done
