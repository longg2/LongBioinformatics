#! /usr/bin/env sh
# Load the variable
export declare out=$1
ncores=$2

LCA() {	# This will perform the LCA analysis based on the Blast Results
	
	# Variables and folders
	local in=$1 # The sample name
	mkdir -p ${out}LCA

	# The actual work
	echo "Getting the counts"

	mapfile -t counts < <(cut -f 1,13 $in | sort --compress-program gzip | uniq -c | sed -e "s/^ *//g" -e "s/ /\t/g" | sort -k3 )

	mapfile -t blast < <(cut -f 3 <(printf "%s\n" "${counts[@]}") | uniq | taxonkit lineage | taxonkit reformat -t -R "NA" | cut -f 1,4)
	#local blast=$(cut -f 3 <(echo ${counts[@]}) | uniq | taxonkit lineage | taxonkit reformat -t -R "NA" | cut -f 1,4)
	join -1 3 -2 1 <(printf "%s\n" "${counts[@]}") <(printf "%s\n" "${blast[@]}") | cut -f 2- -d " " | sed -e "s/ /\t/g" -e "s/;/\t/g" | sort -k 2 > ${out}LCA/$(basename $in)
}

export -f LCA

echo "Determining the LCA for each hit"
# Because taxonkit defaults to 4 cores, we don't want to accidentally swamp the machines
let taxonKitcores=$ncores/8 # 8 Because I'm using two instances of taxonkit per script
parallel -j $taxonKitcores --bar "LCA {}" ::: ${out}/*
