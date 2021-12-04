#!/usr/bin/awk -f

/^@/{# The header
	print $0
	next	
}

{
	split($13,Edit,/:/) # Splitting the edit distance by the colons
	if (Edit[3] <= 2){
		print $0
		}
	}


