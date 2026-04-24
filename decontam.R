# DECONTAMINATION

# Import data
asvs_16s <- read.csv('ASV_all.csv', header = TRUE, row.names = 1)
taxon_16s <- read.csv('Taxon_all.csv', header = TRUE, row.names = 1)
metadata_16s <- read.csv('metadata_16s_clean.csv', header = TRUE, row.names = 1)

# Install packages
install.packages("remotes")
install.packages("devtools")
devtools::install_github("benjjneb/decontam")
remotes::install_github("joey711/phyloseq")
install.packages("dplyr")

# Load packages
library(decontam)
library(phyloseq)
library(remotes)
library(dplyr)

# Convert ASV and taxonomy tables to matrices
asvs_16s <- as.matrix(asvs_16s)
taxon_16s <- as.matrix(taxon_16s)

# Transpose ASV table (ASVs are columns, samples are rows)
asvs_16s <- t(asvs_16s)

# Confirm that all ASVs in asvs_16s are present in taxon_16s
stopifnot(all(colnames(asvs_16s) %in% rownames(taxon_16s)))

# Confirm that all sample names in metadata_16s are present in asvs_16s
stopifnot(all(rownames(metadata_16s) %in% rownames(asvs_16s)))

# Create a phyloseq object
ps_16s <- phyloseq(otu_table(asvs_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata_16s))

# Remove Chloroplast and Mitochondria reads
physeq_16s <- ps_16s %>% subset_taxa(Order != "Chloroplast" | is.na(Order))
physeq_16s <- physeq_16s %>% subset_taxa(Family != "Mitochondria" | is.na(Family))

# Confirm negative controls are present
stopifnot(any(sample_names(physeq_16s) %in% c("neg1", "neg2", "neg3", "neg5")))

# Determine which samples are negative controls
sample_data(physeq_16s)$is.neg <- sample_names(physeq_16s) %in% c("neg1", "neg2", "neg3", "neg5")

# Run decontam using the prevalence method
contamdf.prev <- isContaminant(physeq_16s, method = "prevalence", neg = "is.neg")

# View summary of contaminants
table(contamdf.prev$contaminant)

# List ASVs identified as contaminants
contaminant_asvs <- rownames(contamdf.prev)[contamdf.prev$contaminant == TRUE]
print(contaminant_asvs)