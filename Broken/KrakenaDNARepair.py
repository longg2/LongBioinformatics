#from Bio.Seq import Seq

# The testing variables
fastqSeq='ATCGGTCTATAATCGACTAGCATCGATGCCGGCCGTATACGATGCAAGC' # The test string
kmerLen=25
window=5

# The actual work here
beginString=fastqSeq[0:kmerLen] # This gets the first kmer
endString=fastqSeq[(len(fastqSeq)-kmerLen):len(fastqSeq)] # This gets the last kmer

# Now we need to get all possible substitutions
allTs=[pos for pos, char in enumerate(beginString) if char == "T"]
potentialDeamination=list(filter(lambda index: index <= 4, allTs)) # The filter

# Now to print all the potential combinations of the beginning string

