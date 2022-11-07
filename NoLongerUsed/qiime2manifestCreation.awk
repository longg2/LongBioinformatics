#!/usr/bin/awk -f

# HOW TO USE:
# find ~+/FOLDER/*R1* -type f | THIS FILE > MANIFEST.tab

# Take the results from find and convert them to a tab list
BEGIN{
	OFS = "\t"
	print "sample-id","forward-absolute-filepath","reverse-absolute-filepath"

}
{
	name=$1
	read1=$1
	gsub(/.*\/|_R.*/,"",name) # Getting the sample name
	read2=gensub(/R1/,"R2","g",read1)


	#sub(/\.*$/,"",name)
	print name,read1,read2
}
