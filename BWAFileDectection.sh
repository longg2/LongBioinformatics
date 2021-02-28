# Functions that I will be using here

FileExtraction(){ # Extract Files from an array using results from another array
	local sample=$1 # Basename of the file
	local arrayFiles=$files

	# Finding the indices which have the same samplename
	sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))
#	local numFiles=$(echo "$sampleFiles" | wc -l )
#	#echo "$numFiles"
#
#	if [ $numFiles = 1 ]; then
#		echo "$sample only has one file.  Assuming all Merged"
#		sampleFiles=($( echo $sampleFiles | tr '\012' ' ' ))
#		#FastaorFastq "$folder/$sampleFiles"	
#		## What file is it?
#		#if [ $? == "0" ]; then
#		#	echo "It's a fastq file.  Proceed to mapping."
#		#elif [ $? == "1" ]; then
#		#	echo "It's a fasta file.  Won't be mapping it."
#		#else
#		#	echo "Unknown file.  Halting."
#		#	exit 1
#		#fi
#
#		# Need to test if it's fasta or fastq
#	elif [ $numFiles = 2 ]; then
#		echo "$sample has two files.  Assuming only Paired"
#		sampleFiles=($( echo $sampleFiles | tr '\012' ' ' ))
#
#		#FastaorFastq "$folder/${sampleFiles[0]}"	
#		## What file is it?
#		#if [ $? == "0" ]; then
#		#	echo "It's a fastq file.  Proceed to mapping."
#		#elif [ $? == "1" ]; then
#		#	echo "It's a fasta file.  Won't be mapping it."
#		#else
#		#	echo "Unknown file.  Halting."
#		#	exit 1
#		#fi
#	elif [ $numFiles = 3 ]; then
#		echo "$sample has three files.  Assuming both paired and merged reads"
#		sampleFiles=($( echo $sampleFiles | tr '\012' ' ' ))
#
#		#echo "$folder/${sampleFiles[0]}"
#		#FastaorFastq "$folder/${sampleFiles[0]}"
#		## What file is it?
#		#if [ $? == "0" ]; then
#		#	echo "It's a fastq file.  Proceed to mapping."
#		#elif [ $? == "1" ]; then
#		#	echo "It's a fasta file.  Won't be mapping it."
#		#else
#		#	echo "Unknown file.  Halting."
#		#	exit 1
#		#fi
#	else
#		echo "$sample has more than three files.  Will be skipped."
#	fi

}

DeduplicateArray(){ # Deduplicating the sample names
	local array=("$@") # This feels weird, but, it'll allow you pass an array to it
	#echo "$array"
	local sampleNames=$(for name in ${array[@]}; do tmp=$(echo ${name/_*}); echo ${tmp/%.f*}; done) # Getting only the sample names
	samples=( $(echo ${sampleNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') )

}

FastaorFastq(){ # Figuring out if fasta or fastq
	local tmp=$1
	if file $tmp | grep -q "compressed"; then 
		local firstChar=$(zcat $1 | head -c 1)
	else 
		local firstChar=$(cat $1 | head -c 1)
	fi

	if [ "$firstChar" == "@" ];then
		return 0
	elif [ "$firstChar" == ">" ];then
		return 1
	else
		return 2
	fi
}

# This here will be the main section of the script.  Would be replaced by what I'm currently using.
declare -r folder=$1 # Don't want any shenanigans going on here
declare -r files=$(find TestingArea/* -type f -printf "%f\n") # Making an array of files
#echo "$files"

DeduplicateArray "${files[@]}" # Deduplicating the array.  Outputs to the variable sample

for sample in ${samples[@]}; do # Iterating over an array of Samples
	FileExtraction $sample # Extracting the file names.  Will be saved as $sampleFiles
	unset merged
	unset r1
	unset r2
	
	printf '%s\n' "${sampleFiles[@]}" > .hiddenlist.list
	echo $sample
	# Need to get the files into their correct variables
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
	printf "MERGED:$merged\nR1:$r1\nR2:$r2\n-------\n"

	#sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep -P "$sample" | tr '\012' ' '))
done

rm .hiddenlist.list
