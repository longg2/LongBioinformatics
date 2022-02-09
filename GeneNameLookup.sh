#! /usr/bin/env bash

# Getting the commands I neeed
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
usage() { printf 'NCBI Gene Query
        -i\tThe file with the genes to be searched
	-n\tName of the organism
	-o\tThe Output Prefix
        -h\tShow this help message and exit\n' 1>&2; exit 1; }


##################################
#Default Values
out="${Org}Genes.list"
while getopts "i:o:n:l:h" arg; do
        case $arg in
                i)
			readarray -t genes < ${OPTARG}
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                n)
                        Org=${OPTARG}
                        #echo "Using $ncores"
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
if [ -z ${genes+x} ]
then
	printf "List of genes not provided\n"
	exit 1
fi
# Removing previous gene list with same name to prevent confusion
if [ -f "$out" ]; then
	rm $out
fi

# This is the main workhorse here
((timeEstimate=${#genes[@]}/120))
echo "Pulling out the descriptions of ${#genes[@]} genes from ${Org}."
echo "This will take at least $timeEstimate minute(s)."
total=${#genes[@]}
count=0
ProgressBar $count $total
for gene in ${genes[@]}; do
	esearch -db gene -query "(${gene}) AND ${Org}[Organism]" |\
		efetch -format tabular | grep -v "discontinued" >> $out

	((count=count+1))
	#count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
	sleep 0.5
done

