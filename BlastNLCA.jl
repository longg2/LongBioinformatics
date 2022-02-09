# Testing if julia would be a better idea to do the blastn analysis

# This is where I would load up the packages I need

using CSV
using DataFrames

# This is where I get the variables from the command line

print("Reading the Data.\n")
df = CSV.read("TMP.file", DataFrame, delim = "\t")
ColNames = ["Count", "Sequence", "superkingdom", "phylum","class","order","family","genus","species"]
rename!(df, ColNames)

print("Splitting the data into a dictionary.\n")
dfGroup = groupby(df, :Sequence)
k = keys(dfGroup)

print("Performing the Weighted LCA Search: This can take a while \n")
