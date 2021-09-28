#! /usr/bin/Rscript

suppressPackageStartupMessages(library(methods))
# Loading the Libraries and installing if need be
listPackages <- c("dplyr", "pbapply", "tibble", "parallel")
newPackages <- listPackages[!(listPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
#suppressPackageStartupMessages(library(tibble))
#suppressPackageStartupMessages(library(taxonomizr))
suppressPackageStartupMessages(library(pbapply))

# Need to get those variables inserted
arguments <- commandArgs(trailingOnly = T)
blastTable <- arguments[1]
ncores <- arguments[2]
outName <- arguments[3]
composition <- as.numeric(arguments[4])

# Reading the blast output
cat("Reading the Data.\n")
blastTable <- as_tibble(read.delim(blastTable, header = F, stringsAsFactors = F, na.strings = ""))
colnames(blastTable) <- c("Count", "Sequence", "superkingdom", "phylum","class","order","family","genus","species") 
#blastTable <- blastTable[!is.na(blastTable$genus),] # Less stringent than other options.  Also speeding things up!
invisible(gc())
#taxa <- read.table(taxa,head = F, stringsAsFactors = F)[,1]

cat("Splitting the dataframe into a list.  This should be quick.\n")
blastTable <- split(blastTable, f = blastTable$Sequence)

# This is where it splits into multiple cores
op <- pboptions(type = "timer")
cat("Performing a Weighted LCA Search: This can take a while \n")
WeightedLCATable <- pblapply(blastTable, cl = ncores, function(tmp){
#WeightedLCATable <- do.call(bind_rows,pblapply(blastTable, cl = ncores, function(tmp){
	if(nrow(tmp) == 1){ # If only one
		return(tmp %>% select(-c("Count")))
	}else if(nrow(tmp) == 0){ # If NOTHING
		return(NA)
	}else{
		# Making a weighted LCA table
		thresh <- ceiling(sum(tmp$Count) * composition)
		speciesTest <- as.data.frame(tmp %>% group_by(species) %>% summarize(Count = sum(Count), .groups = "drop_last"))
		if(any(speciesTest$Count >= thresh)){ # If species A OK
			index <- speciesTest %>% filter(Count >= thresh) %>% pull(species)
			tmpDat <- tmp %>% filter(species == index) %>% select(-c("Count")) %>% distinct()
			return(tmpDat)
		}else{ # We dig deeper
			genusTest <- as.data.frame(tmp %>% group_by(genus) %>% summarize(Count = sum(Count), .groups = "drop_last"))
			if(any(genusTest$Count >= thresh)){ # If genus A OK
				index <- genusTest %>% filter(Count >= thresh) %>% pull(genus)
				tmpDat <- tmp %>% filter(genus == index) %>% select(-c("Count", "species")) %>% distinct()
				tmpDat$species <- NA
				return(tmpDat)
			}else{ # Family?
				familyTest <- as.data.frame(tmp %>% group_by(family) %>% summarize(Count = sum(Count), .groups = "drop_last"))
				if(sum(familyTest$Count >= thresh)){
					index <- familyTest %>% filter(Count >= thresh) %>% pull(`family`)
					tmpDat <- tmp %>% filter(`family` == index) %>% select(-c("Count", "species","genus")) %>% distinct()
					tmpDat$species <- NA
					tmpDat$genus <- NA
					return(tmpDat)
				}else{
					tmp[1,3:9] <- NA
					tmp <- as_tibble(tmp[1,-1])
					return(tmp)
				}
			}
		}
	}
})
#}))
cat("The Table is complete, now writing the results!\n")
WeightedLCATable <- WeightedLCATable[!(is.na(WeightedLCATable))]
final.df <- do.call(bind_rows,WeightedLCATable)

write.table(final.df, file = outName, row.names = F, col.names = T, quote = F, sep = "\t")