#! /usr/bin/awk -f
#
{
	if(NR % 4 == 1 || NR % 4 == 3){
		print $0
	}else{
		#totalLength = length($0) - 3*2
		print substr($0,0, length($0) - 12)
	}
}
