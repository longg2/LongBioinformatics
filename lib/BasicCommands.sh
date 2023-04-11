#! /usr/bin/env bash
DeduplicateArray(){ # Deduplicating the sample names
	local array=("$@") # This feels weird, but, it'll allow you pass an array to it
	local sampleNames=$(for name in ${array[@]}; do
	tmp="$(echo $name |sed -E -e 's/_(r1|r2|merged).*//I' -e 's/\.fa.*//I' -e 's/\.fna//I' -e 's/\.fq.*//I' -e 's/_(1|2)\.f.*//I' -e 's/_(1|2)$//I')";
		echo $tmp;
       	done) # Getting only the sample names
	samples=( $(echo ${sampleNames[@]} | tr  ' ' '\n' | uniq | tr '\n' ' ') )

}
FileIdentificationInFunction(){ # Since I can't export Arrays, this here is a workaround using find
	local sample=$1 # Basename of the file
	local location=$2

	# Finding the indices which have the same samplename
	sampleFiles=( $(find -L $location -regex "${sample}\..*|${sample}_.*" -type f -exec basename {} \;) )
	#sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "$sample" | tr '\012' ' '))
}
FileIdentification(){ # Extract Files from an array using results from another array
	local sample=$1 # Basename of the file
	local arrayFiles=$files

	# Finding the indices which have the same samplename
	sampleFiles=($(printf '%s\n' "${arrayFiles[@]}" | grep "${sample}\.|${sample}_" | tr '\012' ' '))
}
FileExtractionInFunction(){ # Assign files to their variables.  Assumes that $sampleFiles and $sample exists
	# Unsetting variables in case they're already defined from a previous run 
	# unset merged
	# unset r1
	# unset r2

	local folderExtraction=$1
	# This here is a fix that's needed if -v doesn't
	merged="NA"
	r1="NA"
	r2="NA"

	# Creating a hidden text file of file names
	#printf '%s\n' "${sampleFiles[@]}" #> .hiddenlist.list

	# Identifying the files
	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r1\.f*|_1\.f*|_r1_0.*|_1$"; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r1\.f*|_1\.f*|_r1_0.*|_1$')
	        r1="$folderExtraction/$fileName"
	fi

	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r2\.f*|_2\.f*|_r2_0.*|2$"; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r2\.f*|_2\.f*|_r2_0.*|2$')
	        r2="$folderExtraction/$fileName"
	fi

	# Two cases for the merged.  Want to control for shenanigans
	if printf '%s\n' "${sampleFiles[@]}" | grep -q -i "$sample\.f.*"; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -v -q "_\.f*") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -i "$sample\.f.*")
	        merged="$folderExtraction/$fileName"
	fi

	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged"; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'merged')
	        merged="$folderExtraction/$fileName"
	fi

}
FileExtraction(){ # Assign files to their variables.  Assumes that $sampleFiles and $sample exists
	# Unsetting variables in case they're already defined from a previous run 
#	unset merged
#	unset r1
#	unset r2

	# This here is a fix that's needed if -v doesn't
	merged="NA"
	r1="NA"
	r2="NA"

	# Creating a hidden text file of file names
	#printf '%s\n' "${sampleFiles[@]}" #> .hiddenlist.list

	# Identifying the files
	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r1\.f*|_1\.f*|_r1_0.*|_1$"; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r1\.f*|_1\.f*|_r1_0.*|_1$')
	        r1="$folder/$fileName"
		unset $filename
	fi

	if printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "r2\.f*|_2\.f*|_r2_0.*|2$"; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'r2\.f*|_2\.f*|_r2_0.*|2$')
	        r2="$folder/$fileName"
		unset $filename
	fi

	# Two cases for the merged.  Want to control for shenanigans
	if printf '%s\n' "${sampleFiles[@]}" | grep -P -q -i "$sample\.f.*"; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -v -q "_\.f*") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -i "$sample\.f.*")
	        merged="$folder/$fileName"
		unset $filename
	elif printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged"; then
	#if [ $(printf '%s\n' "${sampleFiles[@]}" | grep -P -i -q "merged") ]; then
		fileName=$(printf '%s\n' "${sampleFiles[@]}" | grep -P -i 'merged')
	        merged="$folder/$fileName"
		unset $filename
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
GzipDetection(){ # Need to ID Gzipped files and decompress if needed
	local file=$1
	local folder=$2
	local name=${file/%.gz*}

	if file $folder/$file | grep -q "compressed"; then
		gunzip -c $folder/$file > IntGzip/$name
	else
		cp $folder/$file IntGzip/$name
	fi
}
RandomFastaSelection(){
	local fasta=$1
	local subsample=$2
	awk 'BEGIN{RS=">";FS="\n"}NR>1{seq="";for (i=2;i<=NF;i++) seq=seq""$i; print ">"$1"\n"seq}' $fasta | awk '{if ((NR%2) == 0) print prev"\t"$0; prev=$0}' | shuf | head -n $subsample | sed 's/\t/\n/g'
}
