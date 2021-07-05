#! /usr/bin/env bash

KrakenAnalysis(){ # This here will do the heavy lifting for the script

	local sample=$1 # Need to get the sample name
	local folderFunction=$2 # Here we'll get the StringDedup Folder
	local ncores=$3 # How many cores will be used?
	FileIdentificationInFunction $sample $folderFunction
	FileExtractionInFunction $folderFunction

	# Now to actually run Kraken
	if [ "$merged" != "NA" ]; then # If I found a merged file
		kraken2 --db $db --threads $ncores --report ${out}Reports/${sample}Merged.tab \
			--classified-out ${out}Class/${sample}Merged.fastq --unclassified-out ${out}Unclass/${sample}Merged.fastq \
		       --use-names $merged --output ${out}Output/${sample}Merged.out --confidence 0.1
	fi

	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If dealing with a single end library
		kraken2 --db $db --threads $ncores --report ${out}Reports/${sample}Single.tab \
			--classified-out ${out}Class/${sample}Single.fastq --unclassified-out ${out}Unclass/${sample}Single.fastq \
		       --use-names $r1 --output ${out}Output/${sample}Single.out --confidence 0.1
	fi 
	
	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
		kraken2 --db $db --threads $ncores --report ${out}Reports/${sample}Paired.tab \
			--classified-out ${out}Class/${sample}Paired.fastq --unclassified-out ${out}Unclass/${sample}Paired.fastq \
		       --use-names $r1 --output ${out}Output/${sample}Paired.out --confidence 0.1
	fi
		
}
TaxaExtraction(){ # If requested, we need to get the reads which mapped against a particular taxa
	local sample=$1 # Need to get the sample name
	local folderFunction=$2 # Here we'll get the StringDedup Folder
	FileIdentificationInFunction $sample $folderFunction
	FileExtractionInFunction $folderFunction

	if [ "$merged" != "NA" ]; then # If I found a merged file
		extract_kraken_reads.py -k ${out}Output/${sample}Merged.out -r ${out}Reports/${sample}Merged.tab \
			-o ${out}TaxaInterest/$sample.fastq -t $taxa2 --include-children --fastq-output \
			-s $merged
	fi

	if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If dealing with a single end library
		extract_kraken_reads.py -k ${out}Output/${sample}Single.out -r ${out}Reports/${sample}Single.tab \
			-o ${out}TaxaInterest/${sample}_single.fastq -t $taxa2 --include-children --fastq-output \
			-s $r1
	fi 

	if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If paired
			extract_kraken_reads.py -k ${out}Output/${sample}Paired.out -r ${out}Reports/${sample}Paired.tab \
				-o ${out}TaxaInterest/${sample}_r1.fastq -o2 ${out}TaxaInterest/${sample}_r2.fastq \
			       	-t $taxa2 --include-children --fastq-output \
				-s $r1 -s2 $r2
	fi

}
