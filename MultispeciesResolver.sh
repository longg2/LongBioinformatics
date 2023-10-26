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
			blastFile="${OPTARG}"
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                o)
                        out="${OPTARG}"
                        #echo "Settings are being outputted to $log"
                        ;;
		g)
			genus=${OPTARG}
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Check if gene supplied
if [ -z ${blastFile+x} ]
then
	printf "List of proteins not provided\n"
	exit 1
fi
# Removing previous gene list with same name to prevent confusion
if [ -f "$out" ]; then
	rm $out
fi

# Let's first filter the blast file so that only MULTISPECIES are included from the genus we're interested
grep "MULTISPECIES" $blastFile | grep -v "#" | grep "$genus" > MultispeciesBlastOnly.tab

# Next let's get an array to work with which will have the pertinent match information
mapfile -t MultiHits < <(sed -e "s/ref|//" -e "s/|//" MultispeciesBlastOnly.tab | cut -f 1,6,7)

# Want to get the headers setup
printf "#ProteinID\tName\tOrganism\n" > $out

# This is the main workhorse here
((timeEstimate=${#prots[@]}/75))
echo "Pulling out the descriptions of ${#prots[@]} proteins. This will take at least $timeEstimate minute(s)."
total=${#prots[@]}
count=0
ProgressBar $count $total
for hypo in ${MultiHits[@]}; do
	hypoID="$(cut -f 2 $hypo)"

	esearch -db protein -query "${hypo}" |\
		esummary |\
		xtract -pattern DocumentSummary -element Caption,Title,Organism >> $out

	((count=count+1))
	#count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
	sleep 0.75
done
