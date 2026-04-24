# STATISTICAL ANALYSIS

# Install CRAN packages
install.packages(c("ape", "vegan", "cowplot", "tidyr", "dplyr", "compositions", "zCompositions", "viridis", "readxl", "psych", "svglite", "pairwiseAdonis"))

# Install additional packages
install.packages("remotes")
install.packages("tidyverse")
install.packages("zCompositions")
install.packages("vegan")
install.packages("mnormt")
install.packages("BiocManager")
BiocManager::install("phyloseq")

# Install pairwiseAdonis
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Load packages
library(tidyverse)
library(ape)
library(vegan)
library(ggplot2)
library(phyloseq)
library(cowplot)
library(tidyr)
library(dplyr)
library(compositions)
library(zCompositions)
library(viridis)
library(readxl)
library(pairwiseAdonis)
library(mnormt)
library(psych)
library(svglite)
library(phyloseq)

set.seed(336)

# Set working directory
setwd("/storage/work/jzk303/LMU/meatmicro_reads")

# Import data
asvs_16s<-read.csv('ASV_all.csv', header = TRUE, row.names = 1)
taxon_16s<-read.csv('Taxon_all.csv', header = TRUE, row.names = 1)
metadata_16s<-read.csv('metadata.csv', header=TRUE, row.names = 1)

# Ensure taxon_16s is a data frame
taxon_16s <- as.data.frame(taxon_16s)

# Remove ASVs identified as contaminants by decontam (decontam.R)
# Updated vector of ASV names to remove
remove_asvs <- c("ASV734", "ASV763", "ASV2413", "ASV2443", "ASV2896", "ASV3148", "ASV3686", "ASV5157", 
                 "ASV5349", "ASV5634", "ASV5947",  "ASV7350", "ASV7492", "ASV8535", "ASV8561",
                 "ASV8696", "ASV9580", "ASV10030", "ASV16266")

# Vector of sample names to remove
remove_samples <- c("neg1", "neg2", "neg3", "neg5")

# Remove rows (samples) and columns (ASVs)
asvs_16s <- asvs_16s[!(rownames(asvs_16s) %in% remove_samples), 
                     !(colnames(asvs_16s) %in% remove_asvs)]

# Confirm that all sample names in metadata_16s are present in asvs_16s
all(rownames(metadata_16s) %in% rownames(asvs_16s))

# Confirm that the negative controls were removed
any(rownames(asvs_16s) %in% c("neg1", "neg2", "neg3", "neg5"))
any(rownames(metadata_16s) %in% c("neg1", "neg2", "neg3", "neg5"))

# Clean up taxonomy tables
# Add '_unclassified' marker to NAs in the taxonomy table
taxon_16s$Phylum<-ifelse(is.na(taxon_16s$Phylum), paste(taxon_16s$Kingdom, "unclassified", sep = '_'), taxon_16s$Phylum)
taxon_16s$Class<-ifelse(is.na(taxon_16s$Class), paste(taxon_16s$Phylum, "unclassified", sep = '_'), taxon_16s$Class)
taxon_16s$Order<-ifelse(is.na(taxon_16s$Order), paste(taxon_16s$Class, "unclassified", sep = '_'), taxon_16s$Order)
taxon_16s$Family<-ifelse(is.na(taxon_16s$Family), paste(taxon_16s$Order, "unclassified", sep = '_'), taxon_16s$Family)
taxon_16s$Genus<-ifelse(is.na(taxon_16s$Genus), paste(taxon_16s$Family, "unclassified", sep = '_'), taxon_16s$Genus)

# Remove redundant _unclassified
taxon_16s$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Class)
taxon_16s$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Order)
taxon_16s$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Order)
taxon_16s$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Genus)

# Convert ASV and taxonomy tables to matrices
asvs_16s<-as.matrix(asvs_16s)
taxon_16s<-as.matrix(taxon_16s)

# Create a phyloseq object
ps_16s<-phyloseq(otu_table(asvs_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata_16s))

# Remove Chloroplast and Mitochondria reads from ASV table
physeq_16s <- ps_16s %>%  subset_taxa( Order!="Chloroplast" | is.na(Order) )
physeq_16s <- physeq_16s %>% subset_taxa( Family!= "Mitochondria" | is.na(Family))

# Extract the ASV table from phyloseq object
asv_16s<-as.data.frame(t(otu_table(physeq_16s)))
tail(rowSums(asv_16s))

# Extract the taxonomy table from phyloseq object
taxon.16s<-as.matrix(tax_table(physeq_16s))

# Identify ASVs present in both taxon_16s and asv_16s objects
keep_ids <- intersect(rownames(taxon_16s), rownames(asv_16s))

# Subset taxon_16s to those ASVs
taxon_16s_keep <- taxon_16s[keep_ids, , drop = FALSE]

# Sanity checks
length(keep_ids)                                
sum(rownames(taxon_16s) %in% rownames(asv_16s))
dim(taxon_16s_keep)

# Save ASV, taxon, and metadata tables
write.csv(asv_16s, file="ASV_16s_clean.csv")
write.csv(taxon_16s_keep, file="Taxon_16s_clean.csv")
write.csv(metadata_16s, file="metadata_16s_clean.csv")

# DIFFERENTIAL ABUNDANCE ANALYSIS
asv_16s<-read.csv('ASV_16s_clean.csv', header = TRUE, row.names = 1)
taxon_16s<-read.csv('Taxon_16s_clean.csv', header = TRUE, row.names = 1)
metadata_16s<-read.csv('metadata_16s_clean.csv', header=TRUE, row.names = 1)

# Install MaAsLin 3 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biobakery/maaslin3")

# Load library
library(maaslin3)

# Convert three columns in the metadata data frame from plain character strings into ordered factors
# Variables are in columns, samples are in rows
metadata_16s$Listeria.monocytogenes <-
  factor(metadata_16s$Listeria.monocytogenes, levels = c('negative', 'positive'))

metadata_16s$Listeria.spp <-
  factor(metadata_16s$Listeria.spp, levels = c('negative', 'positive'))

metadata_16s$Sampling.site.detail.drain <-
  factor(metadata_16s$Sampling.site.detail.drain, levels = c('Drain', 'Other'))


# Run MaAsLin 3 
set.seed(1)

fit <- maaslin3fit <- maaslin3(
  input_data = asv_16s,
  input_metadata = metadata_16s,
  output = "Maaslin3_Lm",
  formula = "~ Listeria.monocytogenes",
  normalization = "TSS",
  transform = "LOG",
  standardize = TRUE,
  max_significance = 0.1,
  cores = 7
)



fit1 <- maaslin3fit <- maaslin3(
  input_data = asv_16s,
  input_metadata = metadata_16s,
  output = "Maaslin3_Lspp",
  formula = "~ Listeria.spp",
  normalization = "TSS",
  transform = "LOG",
  standardize = TRUE,
  max_significance = 0.1,
  cores = 7
)

# Analysis on a subset of drain samples
fit <- maaslin3fit <- maaslin3(
  input_data = asv_16s,
  input_metadata = metadata_16s,
  output = "Maaslin3_Lm_Drain",
  formula = "~ Listeria.monocytogenes",
  normalization = "TSS",
  transform = "LOG",
  standardize = TRUE,
  max_significance = 0.1,
  cores = 7
)


fit1 <- maaslin3fit <- maaslin3(
  input_data = asv_16s,
  input_metadata = metadata_16s,
  output = "Maaslin3_Lspp_Drain",
  formula = "~ Listeria.spp",
  normalization = "TSS",
  transform = "LOG",
  standardize = TRUE,
  max_significance = 0.1,
  cores = 7
)



# Transpose to have samples in rows, ASVs in columns 
head(asv_16s) 
asv_16s<-t(asv_16s)

# Replace zeros with pseudo-counts
library(zCompositions)
asv.n0_16s<-t(cmultRepl(asv_16s, label=0, method="CZM", z.warning = 0.9, z.delete = 0, output="p-counts")) 

# Replace any negative values with 1e-6
asv_n0_16s <- ifelse(asv.n0_16s < 0, 1e-6, asv.n0_16s)

# Check the format of the new object (samples in columns, ASVs in rows)
head(asv_n0_16s) 

# CLR transformation
asv.n0.clr_16s<-apply(asv_n0_16s, 2, function(x){log(x)-mean(log(x))})
asv.n0.clr_16s<-t(asv.n0.clr_16s)

# Save CLR-transformed table
write.csv(asv.n0.clr_16s, file = "asv.n0.clr_16s.csv")


# IDENTIFY CORE ASVs PER FACILITY
# Transpose to have sample names in columns
asv_16s_prop <- t(asv_16s_prop)

common_samples <- intersect(colnames(asv_16s_prop), rownames(metadata_16s))
asv_prop_sub <- asv_16s_prop[, common_samples]
metadata_sub <- metadata_16s[common_samples, , drop = FALSE]

# Subset taxonomy
taxon_16s_sub <- taxon_16s[rownames(asv_prop_sub), , drop = FALSE]

# Unique facilities 
facilities <- unique(metadata_sub$Facility)

# Store results
core_taxa_by_facility <- list()

for (fac in facilities) {
  # Samples from this facility
  fac_samples <- rownames(metadata_sub)[metadata_sub$Facility == fac]
  
  # Subset abundance matrix
  fac_abund <- asv_prop_sub[, fac_samples, drop = FALSE]
  
  # Identify core ASVs (≥ 0.001 in ≥ 50% of samples)
  core_asvs <- rownames(fac_abund)[rowSums(fac_abund >= 0.0001) >= 0.80 * length(fac_samples)]
  
  # Calculate mean and standard deviation
  if (length(core_asvs) > 0) {
    mean_abund <- rowMeans(fac_abund[core_asvs, , drop = FALSE])
    sd_abund <- apply(fac_abund[core_asvs, , drop = FALSE], 1, sd)
    
    # Combine with Genus
    result <- data.frame(
      ASV = core_asvs,
      Genus = taxon_16s_sub[core_asvs, "Genus", drop = TRUE],
      Mean_Abundance = mean_abund,
      SD_Abundance = sd_abund,
      row.names = NULL
    )
  } else {
    result <- data.frame(
      ASV = character(),
      Genus = character(),
      Mean_Abundance = numeric(),
      SD_Abundance = numeric()
    )
  }
  
  # Save to list and file
  core_taxa_by_facility[[as.character(fac)]] <- result
  write.csv(result,
            file = paste0("core_taxa_facility", fac, ".csv"),
            row.names = FALSE)
}

# Return the object
core_taxa_by_facility

# IDENTIFY TOP 10 TAXA PER FACILITY
# The 10 ASVs that ever reached the highest single‑sample relative abundance within that facility (ranking by max across samples)
common_samples <- intersect(colnames(asv_16s_prop), rownames(metadata_16s))
asv_prop_sub <- asv_16s_prop[, common_samples]
metadata_sub <- metadata_16s[common_samples, , drop = FALSE]

# Subset taxonomy
taxon_16s_sub <- taxon_16s[rownames(asv_prop_sub), , drop = FALSE]

# Unique facilities
facilities <- unique(metadata_sub$Facility)

# Store results
top_taxa_by_facility <- list()

for (fac in facilities) {
  # Get samples for this facility
  fac_samples <- rownames(metadata_sub)[metadata_sub$Facility == fac]
  
  # Subset abundance matrix
  fac_abund <- asv_prop_sub[, fac_samples, drop = FALSE]
  
  # Max abundance per ASV across facility samples
  max_abund <- apply(fac_abund, 1, max)
  
  # Get top 10 ASVs by max abundance
  top_asvs <- names(sort(max_abund, decreasing = TRUE))[1:min(10, length(max_abund))]
  
  # Compute summary stats
  mean_abund <- rowMeans(fac_abund[top_asvs, , drop = FALSE])
  sd_abund <- apply(fac_abund[top_asvs, , drop = FALSE], 1, sd)
  
  # Combine results
  result <- data.frame(
    ASV = top_asvs,
    Genus = taxon_16s_sub[top_asvs, "Genus", drop = TRUE],
    Max_Abundance = max_abund[top_asvs],
    Mean_Abundance = mean_abund,
    SD_Abundance = sd_abund,
    row.names = NULL
  )
  
  # Store and save
  top_taxa_by_facility[[as.character(fac)]] <- result
  write.csv(result,
            file = paste0("top10_taxa_facility.csv", fac, ".csv"),
            row.names = FALSE)
}

# Output object with all results
top_taxa_by_facility

# Collapse ASV relative abundances into genus relative abundances
# Find shared ASVs between abundance and taxonomy tables
shared_asvs <- intersect(rownames(asv_16s_prop), rownames(taxon_16s))

# Subset abundance and taxonomy
asv_prop <- asv_16s_prop[shared_asvs, , drop = FALSE]
taxa_sub <- taxon_16s[shared_asvs, , drop = FALSE]

# Extract and clean genus names
genus_vec <- taxa_sub[, "Genus"]
genus_vec[is.na(genus_vec) | genus_vec == ""] <- "Unclassified"

# Convert abundance matrix to data frame and add genus column
asv_prop_df <- as.data.frame(asv_prop)
asv_prop_df$Genus <- genus_vec

# Sum ASV abundances by genus across all samples
library(dplyr)
genus_prop <- asv_prop_df %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  as.data.frame()

# Set row names to genus and remove genus column
rownames(genus_prop) <- genus_prop$Genus
genus_prop$Genus <- NULL

# Save
write.csv(genus_prop, file = "genus_prop.csv")


# IDENTIFY CORE GENERA BY FACILITY
# Ensure sample names match between genus-level matrix and metadata
common_samples <- intersect(colnames(genus_prop), rownames(metadata_16s))
genus_prop_sub <- genus_prop[, common_samples, drop = FALSE]
metadata_sub <- metadata_16s[common_samples, , drop = FALSE]

# Unique facilities
facilities <- unique(metadata_sub$Facility)

# Store results
core_genera_by_facility <- list()

for (fac in facilities) {
  # Samples from this facility
  fac_samples <- rownames(metadata_sub)[metadata_sub$Facility == fac]
  
  # Subset genus abundance matrix
  fac_abund <- genus_prop_sub[, fac_samples, drop = FALSE]
  
  # Identify core genera
  core_genera <- rownames(fac_abund)[rowSums(fac_abund >= 0.0001) >= 0.90 * length(fac_samples)]
  
  if (length(core_genera) > 0) {
    mean_abund <- rowMeans(fac_abund[core_genera, , drop = FALSE])
    sd_abund <- apply(fac_abund[core_genera, , drop = FALSE], 1, sd)
    
    result <- data.frame(
      Genus = core_genera,
      Mean_Abundance = mean_abund,
      SD_Abundance = sd_abund,
      row.names = NULL
    )
  } else {
    result <- data.frame(
      Genus = character(),
      Mean_Abundance = numeric(),
      SD_Abundance = numeric()
    )
  }
  
  # Store and write
  core_genera_by_facility[[as.character(fac)]] <- result
  write.csv(result,
            file = paste0("core_genera_0.0001_90_facility_", fac, ".csv"),
            row.names = FALSE)
}

# Return the object
core_genera_by_facility



# IDENTIFY CORE GENERA ACROSS ALL SAMPLES
# Ensure sample names match between genus-level matrix and metadata
common_samples <- intersect(colnames(genus_prop), rownames(metadata_16s))
genus_prop_sub <- genus_prop[, common_samples, drop = FALSE]

# Identify core genera
min_prevalence <- 0.9 * ncol(genus_prop_sub)

core_genera <- rownames(genus_prop_sub)[rowSums(genus_prop_sub >= 0.0001) >= min_prevalence]

# Compute mean and SD across all samples
if (length(core_genera) > 0) {
  mean_abund <- rowMeans(genus_prop_sub[core_genera, , drop = FALSE])
  sd_abund <- apply(genus_prop_sub[core_genera, , drop = FALSE], 1, sd)
  
  core_genus_summary <- data.frame(
    Genus = core_genera,
    Mean_Abundance = mean_abund,
    SD_Abundance = sd_abund,
    row.names = NULL
  )
} else {
  core_genus_summary <- data.frame(
    Genus = character(),
    Mean_Abundance = numeric(),
    SD_Abundance = numeric()
  )
}

# Save
write.csv(core_genus_summary, file = "core_genera_all_samples.csv", row.names = FALSE)

# View
core_genus_summary


# PCA
# Perform PCA on the CLR-transformed ASV count matrix
pc.clr_16s<-prcomp(asv.n0.clr_16s)

# Calculate total variance in the data
library(compositions)
mvar.clr_16s<-mvar(asv.n0.clr_16s)

# Extract the first two principal components
library(dplyr)
row_16s<-rownames(asv.n0.clr_16s)
pc_out_16s<-as.data.frame(pc.clr_16s$x[,1:2])

# Load metadata and align it with PCA output using sample names
metadata_16s<-read.csv('metadata_16s_clean.csv', header=TRUE, row.names = 1)

pc_out_16s <- pc_out_16s[match(rownames(metadata_16s), rownames(pc_out_16s)), ]

# Ensure that rows match between PCA and metadata
stopifnot(all(rownames(pc_out_16s) == rownames(metadata_16s)))

# Merge PCA results with metadata for plotting and interpretation
pc_out_meta_16s <- cbind(pc_out_16s, metadata_16s)

# Convert relevant metadata columns to factors for plotting
pc_out_meta_16s$Facility<-as.factor(pc_out_meta_16s$Facility)
pc_out_meta_16s$Species<-as.factor(pc_out_meta_16s$Species)
pc_out_meta_16s$Food.contact.surface<-as.factor(pc_out_meta_16s$Food.contact.surface)
pc_out_meta_16s$Area<-as.factor(pc_out_meta_16s$Area)
pc_out_meta_16s$Listeria.monocytogenes<-as.factor(pc_out_meta_16s$Listeria.monocytogenes)
pc_out_meta_16s$Listeria.spp<-as.factor(pc_out_meta_16s$Listeria.spp)

# PCA plots with first 2 PCs
library(ggplot2)

# Facility
PCA_16s_1 <- ggplot(pc_out_meta_16s, aes(x=PC1,y=PC2, color=Facility))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,linewidth=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  #ggtitle("Bacteria", subtitle = "PCA by Facility and Food contact surface")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16s_1
ggsave("PCA_Facility.png", plot =PCA_16s_1, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Facility.svg", plot =PCA_16s_1, device="svg", width=6, height=5, units="in",dpi=600)

# L. monocytogenes
PCA_16s_2 <- ggplot(pc_out_meta_16s, aes(x=PC1,y=PC2, color=Listeria.monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,linewidth=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  #ggtitle("Bacteria", subtitle = "PCA by Facility and Listeria monocytogenes")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16s_2
ggsave("PCA_Lm.png", plot =PCA_16s_2, device="png", width=10, height=5, units="in",dpi=600)
ggsave("PCA_Lm.svg", plot =PCA_16s_2, device="svg", width=10, height=5, units="in",dpi=600)

# Listeria spp.
PCA_16s_3 <- ggplot(pc_out_meta_16s, aes(x=PC1,y=PC2, color=Listeria.spp))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,linewidth=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  #ggtitle("Bacteria", subtitle = "PCA by Facility and Listeria monocytogenes")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16s_3
ggsave("PCA_Lspp.png", plot =PCA_16s_3, device="png", width=10, height=5, units="in",dpi=600)
ggsave("PCA_Lspp.svg", plot =PCA_16s_3, device="svg", width=10, height=5, units="in",dpi=600)


# PERMANOVA by ASV
#Calculate Aitchinson distance
dist_16s<-dist(asv.n0.clr_16s, method='euclidean')
dist_mat <- as.matrix(dist_16s)
write.csv(dist_mat, file = "dist_mat.csv", row.names = FALSE)

# Convert Facility and Species to factors
metadata_16s$Facility <- as.factor(metadata_16s$Facility)
metadata_16s$Species <- as.factor(metadata_16s$Species)
metadata_16s$Listeria.monocytogenes <- as.factor(metadata_16s$Listeria.monocytogenes)
metadata_16s$Listeria.spp <- as.factor(metadata_16s$Listeria.spp)
metadata_16s$Sampling.site.detail <- as.factor(metadata_16s$Sampling.site.detail)

# Check the number of levels
levels(metadata_16s$Facility)
levels(metadata_16s$Species)
levels(metadata_16s$Listeria.monocytogenes)
levels(metadata_16s$Listeria.spp)
levels(metadata_16s$Sampling.site.detail)

# Reorder metadata to match dist_16s labels
metadata_16s <- metadata_16s[labels(dist_16s), ]

# Sanity check
all(rownames(metadata_16s) == labels(dist_16s))

# Install packages
install.packages("devtools")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Load packages 
library(devtools)
library(pairwiseAdonis)

# Test the main effect of facility + interaction with (animal) species
permanova_16s <- pairwise.adonis2(
  dist_16s ~ Facility * Species,
  data = metadata_16s,
  perm = 999,
  p.adjust.m = "bonferroni"
)

permanova_16s

# Test the main effect of Area + interaction with facility
permanova_16s1 <- pairwise.adonis2(
  dist_16s ~ Area * Facility,
  data = metadata_16s,
  perm = 999,
  p.adjust.m = "bonferroni"
)

permanova_16s1

# Test the main effect of Listeria.monocytogenes + interaction with Facility and Area
permanova_16s3 <- pairwise.adonis2(
  dist_16s ~ Listeria.monocytogenes * Facility * Area,
  data = metadata_16s,
  perm = 999,
  p.adjust.m = "bonferroni"
)

permanova_16s3

# Test the main effect of Listeria.spp + interaction with Facility and Area
permanova_16s4 <- pairwise.adonis2(
  dist_16s ~ Listeria.spp * Facility * Area,
  data = metadata_16s,
  perm = 999,
  p.adjust.m = "bonferroni"
)

permanova_16s4

# Test the main effect of Sampling site detail
permanova_16s5 <- pairwise.adonis2(
  dist_16s ~ Sampling.site.detail,
  data = metadata_16s,
  perm = 999,
  p.adjust.m = "bonferroni"
)

permanova_16s5

# Test the main effect of Listeria.monocytogenes + interaction with Sampling site detail
permanova_16s6 <- pairwise.adonis2(
  dist_16s ~ Listeria.monocytogenes * Sampling.site.detail,
  data = metadata_16s,
  perm = 999,
  p.adjust.m = "bonferroni"
)

permanova_16s6