
SPAdesFunction(){
	spades -1 $r1 -2 $r2 --merged $merged -t $ncores -o ${out}/$sample
}

SPAdesAncientFunctionV2(){
	local sample=$1
	local folderFunction=$2

	FileIdentificationInFunction $sample $folderFunction
	printf "${sampleFiles[0]}\n"
	FileExtractionInFunction $folderFunction

	printf "Forward:\t$r1\nReverse:\t$r2\nMerged:\t$merged\n"

	/home/sam/bin/spades --meta -1 $r1 -2 $r2 --merged $merged -t $ncores -o ${out}/$sample -k 15,21,29 # The kmers are what let it run with the smaller read lengths

}

SPAdesAncientFunction(){
	local sample=$1
	local folderFunction=$2

	FileIdentificationInFunction $sample $folderFunction
	FileExtractionInFunction $folderFunction

	if [ "$merged" != "NA" ]; then
		/home/sam/bin/spades --careful -1 $r1 -2 $r2 --merged $merged -t $ncores -o ${out}/$sample --phred-offset 33 # Might be only due to E. coli samples....
	else
		/home/sam/bin/spades --meta -1 $r1 -2 $r2 -t $ncores -o ${out}/$sample --phred-offset 33 # Might be only due to E. coli samples....

	fi
}
