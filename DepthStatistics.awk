#! /usr/bin/awk -f

#function stdev(mean, len ,values){
#	for(i = 1; i <= len; i++){
#		ss+=(values[i] - mean)^2
#		}
#	return sqrt(ss/len)
#	}
function stdev(sum, sumsq, len){
	sdev=sqrt((sumsq / len) - (sum / len)^2)
	return(sdev)
	}
function mean(sum, len){
	Mean=sum/len
	return(Mean)
	}
function CI(stdev, len){
	Error=1.96 * stdev/sqrt(len) # Looking for a 95% CI
	return(Error)
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
		sumsq=$3^2
		if ($3 > 0){
			PCov = 1
		}else{
			PCov = 0
		}
	}else if($1 == chrom){ # If still in the same fragment
		len+=1
		totDepth+=$3
		sumsq+=$3^2
		if ($3 > 0){ PCov+=1 }
	}else{ # If not, we need to report the mean and PCov
		Mean=mean(totDepth,len)
		if (Mean > 0){ # Don't want to be dividing by zero
			STDEV=stdev(totDepth,sumsq,len)
			CV=STDEV/Mean
			pcov=PCov/len
			Error=CI(STDEV,len)
			print basename(FILENAME),chrom, Mean, STDEV, Error, CV, pcov
		}else{
			print basename(FILENAME), chrom, 0,0,0,"NA",0
		}
		
		chrom=$1
		len=1
		totDepth=$3
		sumsq=$3^2
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
	}else if (mean(totDepth,len) > 0){ # Don't want to be dividing by zero
		Mean=mean(totDepth,len)
		STDEV=stdev(totDepth,sumsq,len)
		CV=STDEV/Mean
		pcov=PCov/len
		Error=CI(STDEV,len)
		print basename(FILENAME),chrom, Mean, STDEV, Error, CV, pcov
	}else{
		print basename(FILENAME), chrom, 0,0,0,"NA",0
	}

}

