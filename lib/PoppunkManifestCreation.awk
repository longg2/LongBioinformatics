#! /usr/bin/awk -f

BEGIN{
	OFS="\t"
}

{ # Want to set this up so that it works for everything
	if(FNR==1){ # If the first run
		sample=$1
		r1=$2
	}else if(sample == $1){ # Dealing with a paired end
		r2=$2
	}else{
		if(r2 == ""){
			print sample,r1
		}else{
			print sample,r1,r2
		}
		sample = $1
		r1 = $2
	}
}

END{ # Need to make sure the final samples squeaks through
	if(r2 == ""){
		print sample,r1
	}else{
		print sample,r1,r2
	}

}
