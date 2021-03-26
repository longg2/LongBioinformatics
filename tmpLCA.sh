# Load the variable
export declare -r out=$1
ncores=$2

LCA() {	# This will perform the LCA analysis based on the Blast Results
	
	# Variables and folders
	local in=$1 # The sample name
	mkdir -p ${out}LCA

	# The actual work
	local counts=$(cut -f 1,13 $in |\
		sort --compress-program gzip |\ # Can be memory inefficient at times 
	       	uniq -c |\
	       	sed -e "s/^ *//g" -e "s/ /\t/g" |\ # Preventing errors
		sort -k3)# Getting the Counts

	local blast=$(cut -f 3 $counts | uniq | taxonkit lineage | taxonkit reformat -t -R "NA" | cut -f 1,4)
	join -1 3 -2 1 <($counts) <($blast) | cut -f 2- -d " " | sed -e "s/ /\t/g" -e "s/;/\t/g" | sort -k 2 > ${out}LCA/$in.tab
}

export -f LCA

echo "Determining the LCA for each hit"
# Because taxonkit defaults to 4 cores, we don't want to accidentally swamp the machines
let taxonKitcores=$ncores/8 # 8 Because I'm using two instances of taxonkit per script
parallel -j $taxonKitcores --bar "LCA {}" ::: ${out}/*
