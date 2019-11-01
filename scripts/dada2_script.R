#!/usr/bin/env Rscript
#######	load libraries #######
library(dada2); packageVersion("dada2")

#######	Set working directory and build file lists #######
path <- getwd()
fnFs <- sort(list.files(path, pattern="_R1.fastq"))	# note the quoted pattern -- if you have
fnRs <- sort(list.files(path, pattern="_R2.fastq"))	# different naming conventions, modify
									# these quoted sections accordingly
sample.names <- sapply(strsplit(fnFs, "_R1.fastq"), `[`, 1)
sample.names
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
fnFs
fnRs

#######	Visualize quality scores #######
p<-1; q<-4		# generates quality plots for samples 1-4. modify to view others.
for (a in p:q)
	{
	pdf(paste("QualPlot",a,".pdf", sep=""))
	print(plotQualityProfile(c(fnFs[a],fnRs[a])))
	dev.off()
	}

#######	Filtering	#######
# dada2 offers extensive filtering options but we've already done all of this
# using general-purpose sequence processing tools.
# here we run the minimum filtering for dada2, removing reads with Ns and truncating
# the beginning of each read to remove primer sequences

filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_R1.filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R2.filt.fastq"))
filtFs
filtRs

filt.out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs,
	trimLeft=c(20,21), compress=F)
# Note that if your primers have a different lengths, you will
# need to change the above (19,27) to the apropriate lengths.
filt.out

#######	Learn error rates	#######
# Error rates are learned by alternating between sample inference (?dada) and
# loess estimation of error rates (?loessErrfun) until convergence
errF <- learnErrors(filtFs)
errR <- learnErrors(filtRs)

#######	Evaluate error rate estimates	#######
# visualize estimated error rates
# error rates should decline with increasing qual score
# red line is based on definition of quality score alone
# black line is estimated error rate after convergence
# dots are observed error rate for each quality score

pdf("ErrorEstimatesF.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf("ErrorEstimatesR.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

#######	Dereplicate reads	#######
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#######	Identify amplicon sequence variants (ASVs) #######
# this setting is recommended by software authors for ITS data. it affects the
# behavior of Needleman-Wunsch alignments for regions with high indel rates, like ITS.
# For 16S data, you should instead use the default (by commenting out the following line)
setDadaOpt(BAND_SIZE=32)

# This is the core dada2 computation, analyzing sequence variation in
# the dereplicated amplicon sequencing reads and inferring composition of
# each sample
dadaFs <- dada(derepFs, err=errF)
dadaRs <- dada(derepRs, err=errR)
dadaFs
dadaRs

#######	Merge paired reads #######
# To further cull spurious sequence variants, we merge the forward and reverse reads
# Paired reads that do not match at the specified level of agreement are removed
# if you have no overlap remaining after trimming, you may have to ignore overlap
# and instead concatenate the sequences (minOverlap=0)
dadaM <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,
	minOverlap=30, maxMismatch=3)

#######	Construct ASV table	#######
# a higher-resolution version of the OTU table produced by clustering approaches
seqtab <- makeSequenceTable(dadaM)
dim(seqtab) 	# this is the nuber of samples and sequence variants

# Inspect and trim distribution of sequence lengths
table(nchar(getSequences(seqtab)))

pdf("ASVlengths.pdf")
plot(table(nchar(getSequences(seqtab))))
dev.off()

# compare this distribution of amplicons with the expected sze (300 for our ITS libraries).
# Some length variation is probably biologically real. You may wish to remove extreme outliers.
# To illustrate how to do so, here we remove anything below 272 bp or above 301 bp.
# I've used 293 and 302 bp range based on the ASV distribution

trim.seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(293,302)]
table(nchar(getSequences(trim.seqtab)))
dim(seqtab)
dim(trim.seqtab) # comparing these tells you how many sequence variants were removed.

#######	Remove chimeras	#######
# The core dada method removes substitution and indel errors, but chimeras may remain.
# Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier
# than it is when dealing with clusters: all sequences which can be exactly reconstructed as
# a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(trim.seqtab, method="consensus", verbose=TRUE)
dim(seqtab.nochim)

# fraction of variants remaining after chimera removal
dim(seqtab.nochim)[2]/dim(trim.seqtab)[2]

# fraction of reads remaining after chimera removal
sum(seqtab.nochim)/sum(trim.seqtab)

# Chimeras may reflect a large fraction of sequence vaariants but are exected to reflect
# a very small fraction of total reads.

# reformat the trimmed ASV table and convert to proportions
renamed<-t(seqtab.nochim)
rownames(renamed)<-c(1:nrow(renamed))
head(renamed)
seqtab.prop<-t(t(renamed)/colSums(renamed))
head(seqtab.prop)
colSums(seqtab.prop)

# make a simple barplot describing these communities
pdf("ASVsBySample.pdf")
b<-barplot(head(seqtab.prop, n=20), col=rainbow(n=20), axisnames=F)
text(colnames(seqtab.prop), x=b, y=-0.1, srt=90, xpd=T)
dev.off()
head(seqtab.prop)

####### Track numbers of reads passing each step of the pipeline #######
# Since we started with HQ reads that had already been extensively filtered,
# we don't expect to lose many reads at any of these steps. If you find that
# you're losing a large fraction of your reads at any step in this pipeline,
# you should investigate that step.
getN <- function(x) sum(getUniques(x))
track <- cbind(filt.out, sapply(dadaFs, getN), sapply(dadaM, getN),
	rowSums(trim.seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track,file="ReadFilterStats.csv",row.names=TRUE,quote=FALSE)

#######	Write outputs for downstream analysis	#######
# There are two main outputs from this pipeline:
# 1. a FASTA file containing the amplicon sequence variants described in the ASV table, and
# 2. The ASV table, tab-delimited text describing the abundance of each ASV in each sample.

uniquesToFasta(seqtab.nochim, "ASVs.fasta", mode = "w")

ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
seqtab.nochim.renamed<-seqtab.nochim
colnames(seqtab.nochim.renamed)<-ids
head(seqtab.nochim.renamed)
write.table(t(seqtab.nochim.renamed),file="ASVtable.txt",quote=F, sep="\t",
	col.names=NA, row.names=TRUE)

####### Save outputs files for R #######
# This way we can come back to the analysis stage at a later point
# or filter ASVs to different lengths if desired
saveRDS(seqtab.nochim, file="seqtab_nochim.rds")
saveRDS(seqtab, file="seqtab.rds")

# If you need to read in previously saved datafiles uncomment and run the following
# lines in R.	You can now run	any futher analyses without re-calculating ASVs.
#seqtab.nochim = readRDS("seqtab_nochim.rds")
#seqtab = readRDS("seqtab.rds")
