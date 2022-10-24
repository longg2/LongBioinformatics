#!/usr/bin/awk -f

function bits2str(bits, data, mask) # Frow the gnu manual does what I'm looking for
{
	    if (bits == 0)
		            return "0"
	        mask = 1
		    for (; bits != 0; bits = rshift(bits, 1))
			            data = (and(bits, mask) ? "1" : "0") data
		    
		        while ((length(data) % 12) != 0)
				        data = "0" data
			    return data
}

# The actual command.  Replaces samflag with a binary string
BEGIN{OFS = "\t"}
{$2=bits2str($2); print}
