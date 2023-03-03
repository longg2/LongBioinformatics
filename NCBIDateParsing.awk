#!/usr/bin/awk -f
function DateParsing(first, second){
	# First, need to test if one of the fields are empty....
	if(first == ""){
		print $1,second
		next
	}else if(second == ""){
		print $1,first
		next
	}

	split(first,field1,"-") # This is the biosample Date
	split(second,field2,"-") # This is the assembly Date

	if(field1[1] < field2[1]){ # Testing the years
		print $1,first
	}else if(field1[1] > field2[1]){
		print $1,second
	}else if(field1[2] < field2[2]){ # Testing the months
		print $1,first
	}else if(field1[2] > field2[2]){
		print $1,second
	}else if(field1[3] < field2[3]){ # Testing the days
		print $1,first
	}else if(field1[3] > field2[3]){
		print $1,second
	}else{ # Otherwise, they must be the same date!
		print $1,first
	}

	}

BEGIN{
	OFS="\t"
	}

NR==FNR{ # This is for the file which contains the collections field
		genomes[$1]	
		if($5 ~ /^[0-9]/){ # If there's an actual date in the collection field, we want that!
			gsub(/T.*/,"",$5) # Not interested in the time zone or time or anything like that
			print $1,$5	
		}else{ # Otherwise, we need to start figuring out which date is the earliest!
			gsub(/T.*/,"",$2) # Not interested in the time zone or time or anything like that
			gsub(/T.*/,"",$3) # Not interested in the time zone or time or anything like that
			DateParsing($2,$3)
		}
	}
NR != FNR{ # Now we're catching everything that didn't have any entry in the collection date field
		if(!($1 in genomes)){
			gsub(/T.*/,"",$2) # Not interested in the time zone or time or anything like that
			gsub(/T.*/,"",$3) # Not interested in the time zone or time or anything like that
			DateParsing($2,$3)
		}
	}

