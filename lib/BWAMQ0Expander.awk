#!/usr/bin/awk -f
# This is an awk script to expand the results from multiple hits
BEGIN{
	FS="\t"
	OFS="\t"
	}

$0 ~ /XA:Z:/{ # We're only interested in those which have the XA flag
	sub("XA:Z:","",$20) # Removing the XA:Z: from the field of interest
	sub("NM:i:","",$12) # Removing the XA:Z: from the field of interest
	split($20,XA, ";") # Let's split the XA variable in to an array
	lowestNM=$12

	for(i in XA){ # First, let's make sure we only get the reads with the loest NM
		split(XA[i],XA_i ,",") # A new sub array!

		if(XA_i[4] < lowestNM){
			lowestNM=XA_i[4]
		}

	}

	for(i in XA){ # Now let's print those new arrays!
		split(XA[i],XA_i ,",") # A new sub array!
		
		# Testing if the NM is greater than our current lowest
		if(XA_i[4] > lowestNM){
			continue
		}		

		# Checking if its a reverse strand mapping
		if(XA_i[2] > 0 & $2 == 0){
			qual=$2 + 16
			sub("-","",XA_i[2])
			qual=$2
		}else{
			qual=$2 + 16
			sub("-","",XA_i[2])
		}

		# Printing the new record
		printf $1"-v"i"\t"qual"\t"XA_i[1]"\t"XA_i[2]"\t"XA_i[3]
		printf "%s-v%s\t%s\t%s\t%s\t%s",$1,i,qual,XA_i[1],XA_i[2],XA_i[3]
		for(x=7;x<=12;x++) printf "%s\t", $x
		printf "NM:i:%s\t",XA_i[4]
		for(x=14;x<=19;x++) printf "%s\t", $x
		printf "\n"


	}

	if($12 <= lowestNM) print $0
}
