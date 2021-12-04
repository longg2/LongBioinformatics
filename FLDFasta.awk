#!/usr/bin/awk -f
/^>/{# Found the header
	if(seqlen){
			print seqlen
		}
		seqlen = 0
		next
}
# Seq Data
{
	seqlen += length($0)
	}
# If we have left over data at the end
	END{if(seqlen){print seqlen}}
