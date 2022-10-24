#!/usr/bin/awk -f

/^@/{# The header
	print $0
	next	
}

{
	if ($5 == 0){
		print $0
		}
	}
