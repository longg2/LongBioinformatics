#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)
export script_full_path

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.

usage() { printf 'Poppunk cluster assigning and updating wrapper v0.1
	Due to memory constraints, you cannot reliably run more than 100
	genomes through Poppunk at a time. This is trying to solve this issue.
	-i\tThe folder the fasta files (REQUIRED)
	-o\tThe final Poppunk model
	-r\tThe original poppunk model
	-g\tMaximum number of genomes per run (Default: 100)
	-n\tNumber of CPU Threads to be used (Default: 8)
	-s\tCharacter to split on (Default: _)
	-l\tLog File Name (Default: $date)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

log() {	printf "BWA mem settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Reference:\t${ref}
	CPU Threads:\t${ncores}
	-------------------------------------
	Max Genomes:\t${gen}
	Sample Delimiter:\t${splitter}
	-------------------------------------\n"; exit 0;
}

######################
### Default values ###
######################

declare -i gen=100
out="PoppunkOut$(date +'%Y%m%d')"
declare -i ncores=8
log="$(date +'%Y%m%d').log"
dedup="FALSE"
splitter="_"

##############
### The UI ###
##############

while getopts "i:s:o:g:r:l:n:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			#declare -r files=$(find -L $folder/* -type f -printf "%f\n") # Making an array of files
                        #echo "The raw sequencing files are located in $in"
                        ;;
                o)
                        #out=${OPTARG}
			declare -r out=$(echo ${OPTARG} | sed "s/\/$//")
                        #echo "$out is the output folder"
                        ;;
		g)
			declare -i gen=${OPTARG}
			#echo "Requiring a minimum mapping quality of $qual"
                        ;;
		s)
			declare -r splitter=${OPTARG}
			#echo "Requiring a minimum mapping quality of $qual"
                        ;;
                r)
			declare -r ref=$(echo ${OPTARG} | sed "s/\/$//")
                        #echo "We're mapping the QCed reads against $ref"
                        ;;
                l)
                        log=${OPTARG}
                        #echo "Settings are being outputted to $log"
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        #echo "Using $ncores CPU Threads"
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

# Testing if files options are missing
if [ -z ${folder+x} ] || [ -z ${ref+x} ]; then
	echo "You are missing either the Input folder or the reference"
	exit 1
fi

#################
### The Setup ###
#################

log | tee $log # The inital log file

# We need to convert the folder of files into a manifest file for poppunk in the form of Sample\t
# R1/Genome\tR2

find -L $folder* |\
       	awk -v s=$splitter '{ind = split($0,file,"/"); split(file[ind],sample,s);print sample[1],$0}' |\
	$script_full_path/lib/PoppunkManifestCreation.awk > AllSamples.list

# We're dealing with two potential situations. One where we have less than the maximum sample count
# and another where we're above it. The former is easy enough (just run poppunk), the other is the reason
# why I've written this whole script

if [ $(cat AllSamples.list | wc -l) -ge $gen ]; then
	printf "$folder has more than 100 samples. Splitting the manifest file into 100 line chunks\n" | tee -a $log

	split -l $gen AllSamples.list SplitSamples 

	# Making it wasy to find what I'm looking for
	splitSamp=$(find ./SplitSamples* -type f -printf "%f\n")
	# Setting up the loop
	total=$(find ./SplitSamples* -type f | wc -l)
	count=0

	cp -r -f $ref old # Preparing for the loop
	for f in old/*; do mv -f $f ${f/$ref/old}; done
	
	ProgressBar $count $total

	for manifest in ${splitSamp[@]}; do
		# Now, because of the way the authors wrote poppunk, we need to start renaming things
		if [ $count -ge 1 ]; then
			rm -rf old
			cp -r newOut old
	
			for f in old/*; do mv -f $f ${f/newOut/old}; done
		fi

		poppunk_assign --query $manifest --update-db --db old --out newOut --write-references --overwrite --threads $ncores 2>> $log

		if [ $? -gt 0 ]; then
			printf "\nPoppunk encoutered an error at $manifest. Please go through the logs to figure out where and rerun from the beginning.\n"
			exit 1
		fi

		count=$(echo "$count + 1" | bc)
		ProgressBar $count $total
	done

	# Now for the final move!

	mv -f newOut $out

	for f in $out/*; do mv -f $f ${f/newOut/$out}; done
	
else
	printf "$folder has more than 100 samples. Splitting the manifest file into 100 line chunks\n" | tee -a $log
	poppunk_assign --query Allsamples.list --update-db --db $ref --out $out --write-references --overwrite --threads $ncores 2>> $log
fi

# Some clean up now
rm -rf SplitSamples* old*
