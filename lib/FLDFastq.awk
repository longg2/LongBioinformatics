#!/usr/bin/awk -f
{ if(NR%4 == 2){
	{print length($0)}
	}
}
