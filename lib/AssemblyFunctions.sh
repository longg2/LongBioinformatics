
SPAdesFunction(){
	spades --isolate -1 $r1 -2 $r2 --merged $merged -t $ncores -o ${out}/$sample
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
