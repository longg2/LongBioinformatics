#! /usr/bin/env bash

# Getting the commands I neeed
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
usage() { printf 'NCBI Protein Query
        -i\tThe file with the genes to be searched
	-o\tThe Output Prefix
        -h\tShow this help message and exit\n' 1>&2; exit 1; }


##################################
#Default Values
out="${Org}Proteins.list"
while getopts "i:o:l:h" arg; do
        case $arg in
                i)
			readarray -t accession < ${OPTARG}
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                o)
                        out="${OPTARG}"
                        #echo "Settings are being outputted to $log"
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Check if gene supplied
if [ -z ${accession+x} ]
then
	printf "List of proteins not provided\n"
	exit 1
fi
# Removing previous gene list with same name to prevent confusion
if [ -f "$out" ]; then
	rm $out
fi

# Want to get the headers setup
#printf "ProteinID\tTaxon\n" > $out

# This is the main workhorse here
((timeEstimate=${#accession[@]}/120))
echo "Pulling out the descriptions of ${#accession[@]} proteins. This will take at least $timeEstimate minute(s)."
total=${#accession[@]}
count=0
ProgressBar $count $total
for acc in ${accession[@]}; do
	esearch -db nuccore -query "${acc}" |\
		esummary |\
		xtract -pattern DocumentSummary -element Caption,TaxId >> $out

	((count=count+1))
	#count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
	sleep 0.5
done
