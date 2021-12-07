library(seqinr) # Reads Fasta Data
library(dplyr)
library(zoo)
library(ggplot2)
library(ggpubr)
#library(minpack.lm) # This is for the Levenberg-Marquardt NLS Algorithm used by Brown et al. UNNEEDED! lm() gets the same answer here

# Making the functions to do the hardwork
CorrectingGC <- function(depthData){ # From the iREP Paper
	GCModel <- lm(Mean ~ GCContent, data = depthData)
	resids <- abs(summary(GCModel)$residuals) # Trying to find the top 1% largest residuals. Don't care about location
	residFilt <- resids < quantile(resids,0.99)
	GCModelFilt <- lm(Mean ~ GCContent, data = depthData[residFilt,])
	adjr2 <- summary(GCModel)$adj.r.squared
	beforeGC <- depthData[residFilt,] %>% ggplot(aes(x = GCContent, y = Mean)) +
		geom_point() + geom_density_2d(colour = "grey") + theme_bw() + geom_smooth(method = "lm") +
	       	ggtitle(bquote("Before GC Correction"~R[adj]^2 == .(round(adjr2, 3))), subtitle = "Top 1% Residuals Filtered") + 
	       	ylab("Mean Coverage") + xlab("GC Content")
	if(adjr2 > 0){ # We make the correction
		covAverage <- mean(depthData$Mean)
		correction <- covAverage - predict(GCModel, data = depthData$GCContent)
		depthData <- depthData %>% mutate(GCCorrected = Mean + correction)
	}else{
		depthData <- depthData %>% mutate(GCCorrected = Mean)
	}
	NewGCModel <- lm(GCCorrected ~ GCContent, data = depthData[residFilt,])
	adjr2 <- summary(NewGCModel)$adj.r.squared
	afterGC <- depthData[residFilt,] %>% ggplot(aes(x = GCContent, y = GCCorrected)) +
		geom_point() + geom_density_2d(colour = "grey") + theme_bw() + geom_smooth(method = "lm") +
	       	ggtitle(bquote("After GC Correction"~R[adj]^2 == .(round(adjr2, 3))), subtitle = "Top 1% Residuals Filtered") +
	       	ylab("Mean Coverage") + xlab("GC Content")
	GCPlots <- ggarrange(beforeGC, afterGC, ncol = 2, align = "hv")
	return(list("Plot" = GCPlots, UpdatedDepths = depthData))

}

# First step is Loading the data
depths <- read.table("ContigsMapping/KaeroEcoliDepths.tab", header = F, col.names = c("Scaffold","Position","Coverage")) %>% as_tibble()
fasta <- read.fasta("ReferenceContigAssembly/contigs.fasta")
depths <- depths %>% split(f = depths$Scaffold) %>% 
                lapply(function(x){x %>% mutate(Base = as.character(fasta[[x$Scaffold[1]]]))}) %>% bind_rows()

PerCoverage <- depths %>% summarize(sum(Coverage > 0)/length(Coverage)) %>% pull() # Needs to be >= 98% to continue
FragmentsperMbp <- depths %>% pull(Scaffold) %>% unique() %>% length() / nrow(depths) * 1e6 # Per Brown et al 2016 needs to be <= 175
# Need to filter the smallest Fragments...
filteredcontigs <- depths %>% count(Scaffold, name = "Length") %>% arrange(Length) %>% filter(Length >=250) %>% pull(Scaffold)

# The authors remove the first and last 100bp of the scaffolds
depthsWithoutEnd <- depths %>% filter(Scaffold %in% filteredcontigs) %>% group_by(Scaffold) %>% filter(`Position` > 100, `Position` < (length(`Position`) - 100)) %>%
	ungroup() %>% mutate(NewPosition = 1:length(`Position`)) %>% mutate(Base = grepl("g|c", Base, ignore.case = T) %>% as.numeric())

# Now to get the rolling mean in here. Want to also do this to the GC Content. We're only interested in every 100bp
RollingMean <- depthsWithoutEnd %>% summarize(Mean = rollmean(Coverage, 5e3), GCContent = rollsum(Base, 5e3)/5e3)
RollingMean <- RollingMean %>% mutate(Position = 1:nrow(RollingMean)) %>% filter(Position %% 100 == 0)
extremes <- RollingMean %>% summarize(abs(Mean/median(Mean)) >= 8) # Need to filter out the extremes (8 fold difference)
RollingMean <- RollingMean[!extremes,]

histPlot <- RollingMean %>% ggplot(aes(x = Mean)) + geom_histogram(bins = 100, colour = "white") + theme_bw() +
	ylab("5KB Windows") + xlab("Mean Read Depth") + ggtitle("Mean Read Depth Histogram")

CorrectedGC <- CorrectingGC(RollingMean) # The GC Correction
CorrectedGC$Plot # Make the plot here

SortedRolling <- CorrectedGC$UpdatedDepths %>% arrange(GCCorrected) %>% mutate(SortedPosition = 1:nrow(RollingMean) * 100)

# Getting the 5% exlcusion of the sorted data
exclude <- quantile(SortedRolling$GCCorrected, c(0.05,0.95))
tmp3 <- SortedRolling %>% filter(GCCorrected >= exclude[1], GCCorrected <= exclude[2])

# Let's fit the nls model of the new data here
linearModel <- lm(log2(GCCorrected) ~ SortedPosition, data = tmp3)
summary(linearModel)
# Need to get the residuals
#lmFUN <- function(xx,parS){return(xx*parS$m + parS$b)}
#residFun <- function(p,observed,xx){observed - lmFUN(xx = xx, parS = p)}
#parStart <- list(b=1, m=1)
#nlFit <- nls.lm(par = parStart, fn = residFun, xx = as.numeric(tmp3$SortedPosition), observed = log2(tmp3$GCCorrected), control = nls.lm.control(nprint=1))
#
#lines(tmp3$SortedPosition, lmFUN(tmp3$SortedPosition,as.list(coef(nlFit))), col = 2, lwd = 2)


# Let's make the linear model here
#test <- predict(linearModel)
adjr2 <- summary(linearModel)$adj.r.squared
#summary(linearModel)
iRep = 2^(coef(linearModel)[2]*nrow(depths))

# Plotting the data
iRepPlot <- ggplot(data = SortedRolling, aes(x = Position, y = log2(GCCorrected))) + geom_line(colour = "grey60", alpha = 0.5) + 
	theme_bw() + geom_line(aes(x = SortedPosition)) + geom_hline(yintercept = log2(exclude), col = "red") +
	geom_line(data = tmp3, aes(x = SortedPosition, y = log2(GCCorrected)), colour = "green") +
	geom_smooth(data = tmp3 %>% filter(GCCorrected > exclude[1], GCCorrected < exclude[2]), aes(x = SortedPosition, y = log2(GCCorrected)), method = "lm", colour = "purple") +
	ylab(bquote(log[2](~"Mean Coverage"~))) +
	ggtitle(label = "Sorted Depth Plot",
		subtitle = bquote("iRep ="~.(round(iRep,2))~R[adj] ^2== .(round(adjr2, 2)) ~"Mean Coverage ="~.(round(mean(SortedRolling$GCCorrected), 2))~"Percent Coverage ="~.(round(PerCoverage * 100,2))~"Fragments/Mbp =" ~ .(FragmentsperMbp)))

pdf("~/Test.pdf", width = 12, height = 9)
print(histPlot)
print(CorrectedGC$Plot)
print(iRepPlot)
dev.off()
