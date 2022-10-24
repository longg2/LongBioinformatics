#!/usr/bin/awk -f
function basename(file, a, n) {
	    n = split(file, a, "/")
	        return a[n]
	}

BEGIN{
	RS=">"
	gene=ARGV[1] # The first argument is the gene I'm interested in
	ARGV[1]="" # Need to remove it so that AWK starts with the "Second File"
	}

$0 ~ gene{
	name=basename(FILENAME)
	split(name,a,".f")
	print ">"a[1]; # Need to get the filename here
	for (i=2; i<=NF; i++) print $i # This here will print all but the gene name
	}
