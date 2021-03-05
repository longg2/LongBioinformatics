#! /usr/bin/env sh
# To be done:
# 	• Remove reliance of seqkit
# 	• Integration (with option for ignoring) the LCA script for blast
######################################
### Functions that I'll be calling ###
######################################

GzipDetection(){ # Need to ID Gzipped files and decompress if needed
	local file=$1
	local folderV2=$2
	local name=${file/%.gz*}

	if file $folderV2/$file | grep -q "compressed"; then
		gunzip -c $folderV2/$file > IntGzip/$name
	else
		cp $folderV2/$file IntGzip/$name
	fi
}
FastaorFastq(){ # Figuring out if fasta or fastq and converting to fasta.
	# All local as I'm avoiding accidental overwrites
	local tmp=$1
	local outFolder=$2
	local firstChar=$(head -c 1 $tmp)
	local name=$(basename $tmp | sed -e "s/\.fastq//" -e "s/\.fq//")
	
	if [ "$firstChar" == ">" ];then
		mv $tmp $outFolder/$(basename $tmp)
	elif [ "$firstChar" == "@" ];then # Need to convert the file to fasta
		sed -n '1~4s/^@/>/p;2~4p' $tmp > $outFolder/$name.fasta
	else
		echo "Skipping $tmp as it's not fasta/fastq" >> FastaFastq.log
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
blastCMD() { # The Meat and Potatoes of the script
	local in=$1
	local sample=$(basename $in .fasta)
	mkdir -p ${out}
		
	echo "Running $blast for $sample"
	if [ "$blast" == "blastn" ]; then
		blastn -db $db -query $in -outfmt "6 std staxid" -evalue $eval -num_threads $ncores -perc_identity $Pident -task blastn > ${out}/$sample.tab 2>> BlastNWarnings.log  &
		pid=$! # Getting the PID of the blast run
		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well

		local nblines=$(wc -l $in | cut -f 1 -d " ")
		echo "Waiting...."
		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
			if [ -s ${out}/$sample.tab ]; then
				local curquery=$(tail -1 ${out}/$sample.tab | cut -f 1)
				local curline=$(fgrep -n $curquery $in |  cut -f 1 -d ':')
				ProgressBar $curline $nblines
				sleep 10
			else
				sleep 60 # Want to give it time to think
			fi
		done
		
		printf "\n$in is Complete\n"


	elif [ "$blast" == "blastp" ]; then
		blastp -db $db -query $in -outfmt "6 std staxid" -evalue $eval -num_threads $ncores > ${out}/$sample.tab 2>> BlastPWarnings.log &
		pid=$! # Getting the PID of the blast run
		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well

		local nblines=$(wc -l $in | cut -f 1 -d " ")
		echo "Waiting...."
		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
			if [ -s ${out}/$sample.tab ]; then
				local curquery=$(tail -1 ${out}/$sample.tab | cut -f 1)
				local curline=$(fgrep -n $curquery $in |  cut -f 1 -d ':')
				ProgressBar $curline $nblines
				sleep 10
			else
				sleep 60 # Want to give it time to think
			fi
		done
		
		printf "\n$in is Complete\n"
	elif [ "$blast" == "blastx" ]; then
		blastx -db $db -query_gencode 11 -query $in -outfmt "6 std staxid" -evalue $eval -num_threads $ncores  > ${out}/${sample}.tab 2>> BlastXWarnings.log &

		pid=$! # Getting the PID of the blast run
		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well

		local nblines=$(wc -l $in | cut -f 1 -d " ")
		echo "Waiting...."
		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
			if [ -s ${out}/$sample.tab ]; then
				local curquery=$(tail -1 ${out}/$sample.tab | cut -f 1)
				local curline=$(fgrep -n $curquery $in |  cut -f 1 -d ':')
				ProgressBar $curline $nblines
				sleep 10
			else
				sleep 60 # Want to give it time to think
			fi
		done
		
		printf "\n$in is Complete\n"
	else
		echo "Please choose either blastn or blastp"
		return 1
	fi

	return 0
}
usage() { printf "BlastN/P Wrapper Script V0.9
	Outputs tab deliminated BlastN/P report file in the form of std staxid.  Taxa counts
	for each step are also outputted.  Need seqkit for subsampling.	String Deduplication occurs.
	TBD: Remove Seqkit Completely, integration of LCA script with BlastN
        -i\tInput Folder
	-o\tOutput Folder (Default = BlastOut)
	-b\tblastn, blastp, or blastx? (Default = blastn)
	-d\tBlast database (Default = /1/scratch/blastdb/nt)
	-e\tMaximum Evalue (Default = 1e-5)
	-k\tRead Length (Default = 30)
	-l\tLog File
	-p\tPercent Identity (Default = 90)
	-s\tReads to Subsample (Default = 0)
	-n\tNumber of CPU Threads to be used (Default = 8)
        -h\tShow this help message and exit\n" 1>&2; exit 1; }

log() {	printf "Blast settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Blast Type:\t${blast}
	Blast Database:\t${db}
	Computer:\t${HOSTNAME}
	-------------------------------------
	E-value:\t${eval}
	Percent Identity:\t${Pident}
	Subsampled Reads:\t${subsample}
	Min Read Length:\t${len}	
	CPU Threads:\t${ncores}\n"; exit 0;
}

export -f FastaorFastq # important as I'm running this in parallel
export -f GzipDetection # important as I'm running this in parallel
#export -f ProgressBar

# Preset Variables
#declare -r folder=$1
#declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files

blast="blastn"
db="/1/scratch/blastdb/nt"
declare -i Pident=90
declare -i len=30
declare -i subsample=0
declare -i ncores=8
eval="1e-5"
log="$(date +'%Y%m%d').log"
out="BlastOut"

while getopts "i:e:d:n:o:b:p:l:s:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
			declare -r files=$(find $folder/* -type f -printf "%f\n") # Making an array of files
                        ;;
		e)
			eval=${OPTARG}
			;;
                d)
                        db=${OPTARG}
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        ;;
		o)
			out=${OPTARG}
			;;
		p)
			declare -i Pident=${OPTARG}
			;;
		b)
			blast=${OPTARG}
			;;
		s)
			declare -i subsample=${OPTARG}
			;;
		l)
			log=${OPTARG}
			;;
		k)
                        declare -i len=${OPTARG}
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

##################
### The Workup ###
##################
log | tee $log

# We need to determine if the file is gzipped
echo "Decompressing the files"
mkdir -p IntGzip
parallel -j $ncores --bar "GzipDetection {} $folder" ::: "${files[@]}"

# Now to detect if the file in IntGzip is fastq, fasta, or other and conver the file
mkdir -p ${out}FastaOnly
echo "Converting FastQ to Fasta"
parallel -j $ncores --bar "FastaorFastq {} ${out}FastaOnly" ::: IntGzip/*
rm -rf IntGzip

# Now to string deduplicate the files as I'd like to speed up the blast runs
mkdir -p ${out}StringDedup/
mkdir -p ${out}prinseqLog
echo "String Deduplciation with prinseq"
parallel -j $ncores --bar "perl /home/sam/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fasta {} -out_good ${out}StringDedup/{/.} -out_bad null -min_len $len -derep 14 -log ${out}prinseqLog/{/.}.log 2> /dev/null" ::: ${out}FastaOnly/*

# The final step is to test if StringDeduplication will occur
if [[ $subsample != 0 ]]; then
	echo "Subsampling the files to ~$subsample reads"
	mkdir -p ${out}Subsample
	parallel -j $ncores --bar "seqkit sample {} -n $subsample > ${out}Subsample/{/} 2> /dev/null" ::: ${out}StringDedup/*
fi

#####################################
### Now to run the Blast commands ###
#####################################

rm -rf ${out} # Delete the folder if it already exists
if [[ $subsample != 0 ]]; then
	for final in ${out}Subsample/*; do
		blastCMD $final	
	done
else
	for final in ${out}StringDedup/*; do
		blastCMD $final	
	done
fi
