#!/usr/bin/awk -f

# Take the results from find and convert them to a tab list
BEGIN{OFS = "\t"}
{
	name=$1
	gsub(/.*\/|\.(fasta|fastq)/,"",name)
	#sub(/\.*$/,"",name)
	print name,$1
}
