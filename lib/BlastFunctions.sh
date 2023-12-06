#! /usr/bin/env bash

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

# I'm not the one who figured out how to guestimate progress... Which biostars forum?
blastCMD() { # The Meat and Potatoes of the script
	local in=$1
	local sample=$(basename $in .fasta)
	#mkdir -p ${out}

	echo "Running $blast for $sample"
	if [ "$blast" == "blastn" ]; then
		blastn -db $db -query $in -outfmt "6 std staxid" -evalue $Eval -num_threads $ncores -perc_identity $Pident -task blastn > ${out}BlastResults/$sample.tab 2>> BlastNWarnings.log  &
		pid=$! # Getting the PID of the blast run
		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well

		local nblines=$(wc -l $in | cut -f 1 -d " ")
		echo "Waiting...."
		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
			if [ -s ${out}BlastResults/$sample.tab ]; then
				local curquery=$(tail -1 ${out}BlastResults/$sample.tab | cut -f 1)
				local curline=$(fgrep -n $curquery $in |  cut -f 1 -d ':')
				ProgressBar $curline $nblines
				sleep 10
			else
				sleep 10 # Want to give it time to think
			fi
		done
		
		printf "\n$in is Complete\n"


	elif [ "$blast" == "blastp" ]; then
		blastp -db $db -query $in -outfmt "6 std staxid" -evalue $Eval -num_threads $ncores > ${out}BlastResults/$sample.tab 2>> BlastPWarnings.log &
		pid=$! # Getting the PID of the blast run
		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well

		local nblines=$(wc -l $in | cut -f 1 -d " ")
		echo "Waiting...."
		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
			if [ -s ${out}BlastResults/$sample.tab ]; then
				local curquery=$(tail -1 ${out}BlastResults/$sample.tab | cut -f 1)
				local curline=$(fgrep -n $curquery $in |  cut -f 1 -d ':')
				ProgressBar $curline $nblines
				sleep 10
			else
				sleep 10 # Want to give it time to think
			fi
		done
		
		printf "\n$in is Complete\n"
	elif [ "$blast" == "blastx" ]; then
		blastx -db $db -query_gencode 11 -query $in -outfmt "6 std staxid" -evalue $Eval -num_threads $ncores  > ${out}BlastResults/${sample}.tab 2>> BlastXWarnings.log &

		pid=$! # Getting the PID of the blast run
		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well

		local nblines=$(wc -l $in | cut -f 1 -d " ")
		echo "Waiting...."
		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
			if [ -s ${out}BlastResults/$sample.tab ]; then
				local curquery=$(tail -1 ${out}BlastResults/$sample.tab | cut -f 1)
				local curline=$(fgrep -n $curquery $in |  cut -f 1 -d ':')
				ProgressBar $curline $nblines
				sleep 10
			else
				sleep 10 # Want to give it time to think
			fi
		done
		
		printf "\n$in is Complete\n"
	else
		echo "Please choose either blastn or blastp"
		return 1
	fi

	return 0
}

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

diamondCMD() { # The Meat and Potatoes of the script
	local in=$1
	local sample=$(basename $in .fasta)
	#mkdir -p ${out}

	echo "Running Diamond for $sample"
	if [ "$blast" == "blastp" ]; then
		diamond blastp --db $db --query $in --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids --evalue $Eval --threads $ncores --out ${out}BlastResults/$sample.tab 2>> BlastPWarnings.log &
		pid=$! # Getting the PID of the blast run
		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well

		local nblines=$(wc -l $in | cut -f 1 -d " ")
		echo "Waiting...."
		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
			if [ -s ${out}BlastResults/$sample.tab ]; then
				local curquery=$(tail -1 ${out}BlastResults/$sample.tab | cut -f 1)
				local curline=$(fgrep -n $curquery $in |  cut -f 1 -d ':')
				ProgressBar $curline $nblines
				sleep 10
			else
				sleep 10 # Want to give it time to think
			fi
		done
		
		printf "\n$in is Complete\n"
	elif [ "$blast" == "blastx" ]; then
		diamond blastx --db $db --query-gencode 11 --query $in --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids -b 8 -c 1 --evalue $Eval --threads $ncores --out ${out}BlastResults/${sample}.tab 2>> BlastXWarnings.log &

		pid=$! # Getting the PID of the blast run
		trap "kill $pid 2> /dev/null" EXIT # If the script is killed, kill the blast run as well

		local nblines=$(wc -l $in | cut -f 1 -d " ")
		echo "Waiting...."
		while kill -0 $pid 2> /dev/null; do # What I want the script to do while I'm waiting
			if [ -s ${out}BlastResults/$sample.tab ]; then
				local curquery=$(tail -1 ${out}BlastResults/$sample.tab | cut -f 1)
				local curline=$(fgrep -n $curquery $in |  cut -f 1 -d ':')
				ProgressBar $curline $nblines
				sleep 10
			else
				sleep 10 # Want to give it time to think
			fi
		done
		
		printf "\n$in is Complete\n"
	else
		echo "Please choose either blastp or blastx"
		return 1
	fi

	return 0
}
