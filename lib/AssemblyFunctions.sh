
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

# Using Shovill
shovillFunction(){
	local sample=$1
#	local folder=$2
#	local out=$3

	# This is being done to facilitate parallelization
	FileIdentificationInFunction $sample $folder
	FileExtractionInFunction $folder

	shovill --R1 $r1 --R2 $r2 --cpus 8 --outdir ${out}/$sample 2>&1 > $out/Logs/${sample}.log

}

#Successive Pilon Improvements
pilonFunction(){
	local sample=$1
	local iteration=$2
#	local folder=$2
#	local out=$3

	# Now for other local functions to make life easier
	local WorkFolder="${out}/$sample/Refining/latest"

	# This is being done to facilitate parallelization
	FileIdentificationInFunction $sample $folder
	FileExtractionInFunction $folder

	# We need to make sure that we're doing the work in the right folder!
	mkdir -p ${out}/$sample/Refining/$iteration
	ln -fns $iteration $WorkFolder

	# Getting the genome from the right spot
	if [[ $iteration -eq 1 ]]; then
		ln -fs ../../contigs.fa $WorkFolder/InputGenome.fasta
	else
		ln -fs ../$((iteration-1))/pilon.fasta $WorkFolder/InputGenome.fasta
	fi

	# Now we need to make the bwa index
	mkdir -p $WorkFolder/Reference
	bwa index -p $WorkFolder/Reference/InputGenome $WorkFolder/InputGenome.fasta
	
	# Doing the mapping and deduplicating in one fell swoop
	bwa mem $WorkFolder/Reference/InputGenome $r1 $r2 -t 8 |\
		samtools view -b -h | samtools sort -n |\
	        samtools fixmate -m - - | samtools sort - |\
 	        samtools markdup -r -S - $WorkFolder/MappedReads.bam

	samtools index $WorkFolder/MappedReads.bam

	# Running Pilon
	pilon --genome $WorkFolder/InputGenome.fasta --bam $WorkFolder/MappedReads.bam --outdir $WorkFolder/ --output pilon
}
