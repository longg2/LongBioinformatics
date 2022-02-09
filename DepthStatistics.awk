#!/usr/bin/awk -f
function stdev(mean, len ,values){
	for(i = 1; i <= len; i++){
		ss+=(values[i] - mean)^2
		}
	return sqrt(ss/len)
	}

function basename(file, a, n) {
	    n = split(file, a, "/")
	        return a[n]
	}

BEGIN{
	OFS="\t"
	#print "File", "Fragment", "Mean", "Std. Dev", "CV","Percent Coverage"
	}

{
	if (NR == 1){ # Setting up if its the first line of the file
		chrom=$1
		len=1
		totDepth=$3
		valarray[$2] = $3
		if ($3 > 0){
			PCov = 1
		}else{
			PCov = 0
		}
	}else if($1 == chrom){ # If still in the same fragment
		len+=1
		totDepth+=$3
		valarray[$2] = $3
		if ($3 > 0){ PCov+=1 }
	}else{ # If not, we need to report the mean and PCov
		if (totDepth/len > 0){ # Don't want to be dividing by zero
			print basename(FILENAME),chrom, totDepth/len, stdev(totDepth/len, len, valarray), stdev(totDepth/len, len,valarray) /(totDepth/len) , PCov/len
		}else{
			print basename(FILENAME), chrom, 0,0,NA,0
		}
		
		delete valarray 
		chrom=$1
		len=1
		totDepth=$3
		valarray[$2] = $3
		if ($3 > 0){
			PCov = 1
		}else{
			PCov = 0
		}

	}
}
END{
	if (FNR == 0){ # If the file is empty we want out
		exit
	}else if (totDepth/len > 0){ # Don't want to be dividing by zero
		print basename(FILENAME),chrom, totDepth/len, stdev(totDepth/len, len, valarray), stdev(totDepth/len, len,valarray) /(totDepth/len) , PCov/len
	}else{
		print basename(FILENAME), chrom, 0,0,NA,0
	}

}

