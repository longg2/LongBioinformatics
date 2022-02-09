#! /usr/bin/Rscript

suppressPackageStartupMessages(library(methods))
# Loading the Libraries and installing if need be
listPackages <- c("dplyr", "pbapply", "tibble", "parallel")
newPackages <- listPackages[!(listPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
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
columnNames <- c("Count", "Sequence", "superkingdom", "phylum","class","order","family","genus","species") 
colnames(blastTable) <- columnNames 
#blastTable <- blastTable[!is.na(blastTable$genus),] # Less stringent than other options.  Also speeding things up!
#taxa <- read.table(taxa,head = F, stringsAsFactors = F)[,1]

cat("Splitting the dataframe into a list.  This should be quick.\n")
blastTable <- split(blastTable, f = blastTable$Sequence)
invisible(gc())

# This is where it splits into multiple cores
op <- pboptions(type = "timer")
cat("Performing a Weighted LCA Search: This can take a while \n")
WeightedLCATable <- pblapply(blastTable, cl = ncores, function(tmp){
	if(nrow(tmp) == 1){ # If only one
		return(tmp %>% select(-c("Count")))
	}else if(nrow(tmp) == 0){ # If NOTHING
		return(NA)
	}else{
		# Making the while loop
		found <- F
		colIndex <- 9
		thresh <- ceiling(sum(tmp$Count) * composition)
		while(found == F){
			summedCount <- as.data.frame(tmp %>% group_by_at(columnNames[colIndex]) %>%
				      summarize(Count =sum(Count), .groups = "drop_last"))
			if(any(summedCount$Count >= thresh)){
				found <- T
			}else if(colIndex != 3){
				colIndex <- colIndex - 1
			}else{ # If we didn't find anything
				colIndex <- -1
				found <- T
			}
		}
		 # When the loop is finished
		if(colIndex == 9){ # If it was the species that we found it in
			index <- summedCount %>% filter(Count >= thresh) %>% pull(species)
			tmpDat <- tmp %>% filter(species == index) %>% select(-c("Count")) %>% distinct()
			return(tmpDat)
		}else if(colIndex == -1){ # If we found nothing matching our specifications
			tmp[1,3:9] <- NA
			tmp <- as_tibble(tmp[1,-1])
			return(tmp)
		}else{# If we found at least something
			index <- summedCount %>% filter(Count >= thresh) %>% pull(get(columnNames[colIndex]))
			tmpDat <- tmp %>% filter(get(columnNames[colIndex]) == index)
			tmpDat[1,columnNames[(colIndex + 1):9]] <- NA
			return(tmpDat[1,-1])
		}

	}
})

cat("The Table is complete, now writing the results!\n")
WeightedLCATable <- WeightedLCATable[!(is.na(WeightedLCATable))]
final.df <- do.call(bind_rows,WeightedLCATable)

write.table(final.df, file = outName, row.names = F, col.names = T, quote = F, sep = "\t")
