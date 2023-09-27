## Original codes are from https://benjjneb.github.io/dada2/tutorial.html
## The following code was modified by Youn-Jeong Choi (2020)
#load dada2 package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")
library(dada2); packageVersion("dada2")

#Change to the directory containing files.
path1 <- "submission/WT_KO_I/R1" 
path2 <- "submission/WT_KO_I/R2"
list.files(path1)
list.files(path2)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path1, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path2, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path1, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path2, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# filter and trim sequences and On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

#Learn error rate
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Merge R1 &R2 reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=30, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# make sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", minFoldParentOverAbundance = 8, multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "WT_KO_I_data/track.csv")

# Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v138_train_set.fa.gz", multithread=TRUE,tryRC=TRUE, minBoot=80)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.fa.gz", allowMultiple=TRUE, verbose=TRUE, tryRC=TRUE)
head(unname(taxa))
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
table(is.na(taxa.print[,"Genus"]))
write.csv(taxa.print, "WT_KO_I_data/taxa.print.csv")

# Adding ASV to sequences
rownames(taxa.print) <- paste0("ASV_",1:nrow(taxa.print))

identical(rownames(taxa.print), rownames(asv_tax))
identical(rownames(taxa.print), rownames(asv_tab))

write.table(taxa.print, "WT_KO_I_data/ASVs_taxonomy.txt", sep="\t", col.names=NA, row.names=T, quote=F)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "WT_KO_I_data/ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "WT_KO_I_data/ASVs_counts.txt", sep="\t", col.names=NA, row.names=T, quote=F)

## transfer to wynton for qiime. The end of dada2