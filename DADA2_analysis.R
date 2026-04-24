# Analysis of microbial communities in meat processing facilities
# DADA2 analysis is based on the https://benjjneb.github.io/dada2/tutorial.html

# R version 4.5.0

# SEQUENCE ANALYSIS

# Install packages
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install(c("zlibbioc", "Biostrings", "dada2"))
install.packages("ggplot2")
install.packages("deldir")

# Load DADA2 packages
library(zlibbioc)
library(deldir)
library(dada2)
library(ggplot2)

# Set working directory
setwd('/storage/home/jzk303/work/LMU/meatmicro_reads')

# Path to input files
path <- '/storage/home/jzk303/work/LMU/meatmicro_reads' 
list.files(path)

# Forward and reverse fastq filenames
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Check read quality for forward and reverse reads
qc_F<-plotQualityProfile(fnFs[1:5])
ggsave("qc_F.png", plot=qc_F, device='png')
qc_F

qc_R<-plotQualityProfile(fnRs[1:5])
ggsave("qc_R.png", plot=qc_R, device='png')

# Save filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter sequences
out <- filterAndTrim(
  fwd         = fnFs, 
  filt        = filtFs,
  rev         = fnRs, 
  filt.rev    = filtRs,
  maxN        = 0,
  maxEE       = c(2, 6),
  truncQ      = 2,
  rm.phix     = TRUE,
  compress    = TRUE,
  multithread = TRUE
)

head(out)

# Ensure sample names line up
basename(fnFs)
basename(fnRs)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Visualize estimated error rates
error_F<-plotErrors(errF, nominalQ=TRUE) 
error_R<-plotErrors(errR, nominalQ=TRUE)
ggsave("error_F.png", plot=error_F, device='png')
ggsave("error_R.png", plot=error_R, device='png')

# ASV inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Check the returned DADA object
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Check the mergers data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Check the distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove sequences outside of the expected length
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 402:431]

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "/storage/work/jzk303/LMU/meatmicro_reads/track.csv", row.names = TRUE)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/storage/home/jzk303/work/LMU/meatmicro_reads/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# Load libraries
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

# Create a phyloseq object
ps_all <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   tax_table(taxa))

# Replace sequence names with ASV labels
dna_all <- Biostrings::DNAStringSet(taxa_names(ps_all))
names(dna_all) <- taxa_names(ps_all)
ps_all <- merge_phyloseq(ps_all, dna_all)
taxa_names(ps_all) <- paste0("ASV", seq(ntaxa(ps_all)))
ps_all

# Save data
asvs_all<-as.data.frame(otu_table(ps_all))
Taxon_all<-as.data.frame(tax_table(ps_all))
DNA<-as.data.frame(refseq(ps_all))

write.csv(asvs_all, file = "ASV_all.csv")
write.csv(Taxon_all, file = "Taxon_all.csv")
write.csv(DNA, file="DNA_all.csv")
