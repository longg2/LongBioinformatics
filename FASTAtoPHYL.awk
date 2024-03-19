#! /usr/bin/awk -f

# This file converts a fasta alignment into a relaxed phylip file.
# Order of the file will be flipped at the end as it is assumed that order doesn't matter!

## Functions being used
function revSort(file){
	cmd="tac tmp.tmp > "file
	system(cmd)
	}

function basename(file, a, n) {
	n = split(file, a, "/")
	return a[n]
}


BEGIN{
	RS=">"
	OFS="      "
}
BEGINFILE{

	outputFile = FILENAME
	gsub(/.fa.*/,".phy",outputFile)
}

{
	# This does the printing
	len = 0
	printf $1	> "tmp.tmp"
	printf OFS > "tmp.tmp"
	for (i=2; i<=NF;i++){
		len+=length($i)
		printf $i > "tmp.tmp"
		}
	printf "\n" > "tmp.tmp"
}

ENDFILE{
	# Some final sorting
	printf FNR-1" "len"\n" > "tmp.tmp"
	revSort(outputFile) # Sequence order shouldn't matter, so we're just going to flip the file

	printf "\n"basename(FILENAME)"\n"
	print "Number of Sequences = ",FNR - 1
	print "Length of Sequences = ",len
}
END{
	system("rm -f tmp.tmp")
	}
