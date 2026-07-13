####R SESSION SETUP####

# This section was adapted from the DADA2 tutorial:
# https://benjjneb.github.io/dada2/tutorial.html

# Load packages
library(phyloseq)
library(devtools)
library(ggplot2)
library(RColorBrewer)
library(permute)
library(lattice)
library(vegan)
library(reshape2)
library(agricolae)
library(ggfortify)
library(ape)
library(phangorn)
library(pheatmap)
library(qdapTools)
library(gridExtra)
library(plyr)
library(cowplot)
library(utils)
library(Rcpp)
library(dada2)

###PRIMARY SEQUENCE ANALYSIS###

##PREPARE THE FILES###
# Set the working directory where the input files are located and output files will be saved
setwd("Working/directory")

# Define the directory containing the samples to be analyzed
path_bact <- "Samples/directory" 

# Define variables containing the forward and reverse read filenames
fnFs_bact2 <- sort(list.files(path_bact, 
                              pattern="_1.fastq.gz", full.names = TRUE))
fnRs_bact2 <- sort(list.files(path_bact, 
                              pattern="_2.fastq.gz", full.names = TRUE))
sample.names_bact <- sapply(strsplit(basename(fnFs_bact2), "_"), `[`, 1)

# Get the number of forward FASTQ files
num_files <- length(fnFs_bact2) # Check that the expected samples are being analyzed
print(num_files)

# Inspect sequence quality
# Modify the number of samples to visualize if needed
a_bact <- plotQualityProfile(fnFs_bact2[1:5])
b_bact <- plotQualityProfile(fnRs_bact2[1:5])
gridExtra::grid.arrange(a_bact, b_bact) # Display forward and reverse quality profiles

# Filter sequences
# Define the output filenames for filtered reads
filtFs_bact2 <- file.path(path_bact, "filtered", 
                          paste0(sample.names_bact, "_F_filt.fastq.gz"))
filtRs_bact2 <- file.path(path_bact, "filtered", 
                          paste0(sample.names_bact, "_R_filt.fastq.gz"))

# Filter sequences according to the quality profiles
# If too few reads pass the filter, slightly increase maxEE,
# especially for reverse reads (e.g., 2.5)
# truncLen can also be reduced, but avoid excessive trimming
# since sufficient overlap between paired reads is required
out_bact2 <- filterAndTrim(fnFs_bact2, filtFs_bact2, fnRs_bact2, 
                           filtRs_bact2, truncLen= c(250,240),
                           maxN=0, maxEE=c(2,2),truncQ=2, rm.phix=TRUE, trimLeft = 0,
                           compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out_bact2)

# The following code can be used if filtering is not required
# (Do not use if the previous filtering step has been performed)
#out_bact_F <- file.copy(from = fnFs_bact, to = filtFs_bact)
#out_bact_R <- file.copy(from = fnRs_bact, to = filtRs_bact)

# DADA2 step 1: Error rate learning
# DADA2 builds an error model by estimating sequencing error rates from the data
errF_bact2 <- learnErrors(filtFs_bact2, multithread=TRUE)
errR_bact2 <- learnErrors(filtRs_bact2, multithread=TRUE)

# Plot the estimated error rates to verify that they follow the expected trend
# Error rates should decrease as quality scores increase
plotErrors(errF_bact2, nominalQ=TRUE)

# DADA2 step 2: Dereplication
# Collapse identical sequencing reads into unique sequences
derepFs_bact2 <- derepFastq(filtFs_bact2, verbose=TRUE)
derepRs_bact2 <- derepFastq(filtRs_bact2, verbose=TRUE)

# DADA2 step 3: ASV inference
# Infer ASVs using the unique sequences and the learned error model
dadaFs_bact2 <- dada(derepFs_bact2, err=errF_bact2, multithread=FALSE)
dadaRs_bact2 <- dada(derepRs_bact2, err=errR_bact2, multithread=FALSE)
dadaFs_bact2[[1]]
dadaRs_bact2[[1]]

# Merge paired-end reads
mergers_bact <- mergePairs(dadaFs_bact2, derepFs_bact2, 
                            dadaRs_bact2, derepRs_bact2, verbose=TRUE)
head(mergers_bact2[[1]])

# Construct the ASV sequence table
seqtab_bact <- makeSequenceTable(mergers_bact)
??makeSequenceTable
dim(seqtab_bact)
table(nchar(getSequences(seqtab_bact)))

# Remove chimeric sequences
seqtab.nochim_bact <- removeBimeraDenovo(seqtab_bact, 
                                         method="consensus", multithread=FALSE, verbose=TRUE)
??removeBimeraDenovo
dim(seqtab.nochim_bact)

# Check the proportion of non-chimeric reads (this value should be high)
sum(seqtab.nochim_bact)/sum(seqtab_bact)

# Track read counts throughout the pipeline
# A high proportion of reads should successfully merge
# Poor merging may indicate that filtering was either too stringent or too permissive
# The goal is to balance the number of retained reads with sequence quality
getN <- function(x) sum(getUniques(x))
track_bact <- cbind(out_bact2, 
                    sapply(dadaFs_bact2, getN), sapply(dadaRs_bact2, getN), 
                    sapply(mergers_bact2, getN), rowSums(seqtab.nochim_bact))
colnames(track_bact) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_bact) <- sample.names_bact
head(track_bact)

# Export the read tracking table
write.table(track_bact, "track_bact.txt",
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)