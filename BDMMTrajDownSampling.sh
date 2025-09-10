#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)
source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.

BDMMTrajDownSampling(){
	local trajfile=$1
	local out=$2
	local burnin=$3
	local downsample=$4
	
	# Performing the downsampling
	awk -v burnin="$burnin" -v downsample="$downsample" 'BEGIN{OFS = "\t"; FS = "\t"} NR==1{print $0; next} NR>1{if($1 >= burnin && $1 % downsample == 0){print $0}}' $trajfile |\
	       	gzip > ${out}.traj.gz
}

usage() { printf "BDMM Trajectory Downsampling
	-i\tThe trajectory file. Can be gzipped (REQUIRED)
	-o\tOutput prefix
	-n\tLength of chain
	-b\tBurn-in fraction (Default: $burninFrac)
	-s\tDownsample number (Default: $samples)
	-l\tLog File Name (Default: $log)
        -h\tShow this help message and exit\n" 1>&2; exit 1; }

log() {	printf "BDMM downsampling settings for $(log):
	Log File:\t${log}
	Input file:\t${file}
	Output prefix:\t${out}
	-------------------------------------
	Length of Run:\t${length}
	Burn-in:\t${burninFrac}
	DownsampleSize:\t${samples}
	-------------------------------------\n"; exit 0;
}

######################
### Default values ###
######################
out="DownsampledBDMMTraj"
log="$(date +'%Y%m%d').log"
declare samples=10000
declare burninFrac=0.1

##############
### The UI ###
##############

while getopts "i:o:l:b:s:n:h" arg; do
        case $arg in
                i)
                        declare -r file=${OPTARG}
                        ;;
                o)
                        out=${OPTARG}
                        #echo "$out is the output folder"
                        ;;
                n)
                        declare length=${OPTARG}
                        ;;
		b)
			declare burninFrac=${OPTARG}
			;;
                s)
                        declare samples=${OPTARG}
                        ;;
                l)
                        log=${OPTARG}
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

#################
### The Setup ###
#################

# Testing if files/options are missing
if [ -z ${file+x} ] || [ -z ${length+x} ]; then
	echo "You are missing either the trajectory or the number of chains"
	exit 1
fi

#if [ ${burninFrac} -ge 1 ]; then
#	echo "The burn-in fraction is >= 1! Please choose a smaller fraction"
#	exit 1
#fi

#log | tee $log # The inital log file

##################
### The Script ###
##################

# Getting the burnin length and downsampling
burnin=$(echo "scale=0;$length * $burninFrac" | bc)
downsample=$(echo "scale=0;$length / $samples" | bc)

# Next, we need to test if the file is already compressed
if file $file | grep -q "compressed"; then
	BDMMTrajDownSampling <(zcat $file) $out $burnin $downsample
else
	BDMMTrajDownSampling $file $out $burnin $downsample
fi

