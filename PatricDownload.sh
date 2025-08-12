#! /usr/bin/env bash
# These are the files and variables that will be needed.
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)
export script_full_path

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.

usage() { printf 'Patric Genome Download
        -i\tThe Metadata
	-o\tThe Output Folder Suffix
        -h\tShow this help message and exit\n' 1>&2; exit 1; }
log() {	printf "MapDamage Settings for $(date):
	Log File:\t${log}
	Input folder:\t${in}
	Output Prefix:\t${out}
	-------------------------------------\n"; exit 0;
}

ncores=2
log="$(date +'%Y%m%d').log"
while getopts "i:l:o:r:n:h" arg; do
        case $arg in
                i)
                        in=${OPTARG}
                        ;;
                o)
                        Outfolder=${OPTARG}
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Getting the folders setup
mkdir -p $Outfolder

# Getting the genomes extracted from the Metadata
#local sampleArray=( $(cut -f 1 -d "," $in | tail -n +2 ))
cut -f 1 -d "," $in | tail -n +2 | tr -d '"'  > tmp.txt # We're extracting the ID codes for each genome here
readarray -t sampleArray < tmp.txt # Putting it in an array

# The loop
cd $Outfolder

total=${#sampleArray[@]}
count=0

ProgressBar $count $total
for i in ${sampleArray[@]}; do
	wget -nv "ftp://ftp.bvbrc.org/genomes/$i/$i.fna"
       	wget -nv "ftp://ftp.bvbrc.org/genomes/$i/$i.PATRIC.faa"
       	wget -nv "ftp://ftp.bvbrc.org/genomes/$i/$i.PATRIC.ffn"
       	wget -nv "ftp://ftp.bvbrc.org/genomes/$i/$i.PATRIC.gff"
	sleep 6

	printf "\n"
	count=$(echo "$count + 1" | bc)
	ProgressBar $count $total
	printf "\n"
done

# Going back to the original folder
printf "\n"
cd -
rm -f tmp
