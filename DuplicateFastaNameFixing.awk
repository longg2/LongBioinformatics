#!/usr/bin/awk -f
BEGIN{
	OFS = "\t"
	}
{
	if(FNR == 1){ # If the first file there's no point in complicating things
		a[$1] = 1
		print $1"-V"a[$1],$2
	}else if(a[$1] >= 1){ # If we find a match
		a[$1] += 1 # Need to increment the number of matches here
		print $1"-V"a[$1],$2
	}else{ #If no match
		a[$1] = 1
		print $1"-V"a[$1],$2
	}
}

