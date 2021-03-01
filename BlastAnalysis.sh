#! /usr/bin/env sh
######################################
### Functions that I'll be calling ###
######################################

GzipDetection(){ # Need to ID Gzipped files and decompress if needed

}

FastaorFastq(){ # Figuring out if fasta or fastq
	local tmp=$1

	if [ "$firstChar" == "@" ];then
		exit 0
	elif [ "$firstChar" == ">" ];then
		exit 1
	else
		exit 2
	fi
}

blastCMD() {
	if [ "$blast" == "blastn" ]; then
		blastn -db $db -query tmp.fasta -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores -perc_identity $Pident
	elif [ "$blast" == "blastp" ]; then
		blastp -db $db -query tmp.fasta -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores
	elif [ "$blast" == "blastx" ]; then
		blastx -db $db -query_gencode 11 -query tmp.fasta -out ${out}/$sample.tab -outfmt "6 std staxid" -evalue $eval -num_threads $ncores 
	else
		echo "Please choose either blastn or blastp"
		return 1
	fi

	return 0
}

