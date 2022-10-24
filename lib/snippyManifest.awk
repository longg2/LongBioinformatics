#!/usr/bin/awk -f

# Take the results from find and convert them to a tab list
BEGIN{OFS = "\t"}
{
	name=$1
	gsub(/.*\/|\.(fna|fasta|fastq|bam)/,"",name)
	gsub(/\.(gz)/,"",name) # Just in case it's compressed...
	#sub(/\.*$/,"",name)
	print name,$1
}
