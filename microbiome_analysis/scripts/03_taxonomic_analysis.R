###TERTIARY SEQUENCE ANALYSIS###

# Load the files required to create the phyloseq object
# (OTU table, taxonomy table, and metadata)
otu_bact <- read.table("ASV_bact.txt", header = TRUE, row.names = 1)
tax_bact <- as.matrix(read.table("taxa_bact.txt", header = TRUE, row.names = 1, sep = "\t"))
meta <- import_qiime_sample_data("meta.txt")
meta <- data.frame(meta)

# Run the following code if you want to create a new "Species" column
# combining the genus and species names
tax_bact <- data.frame(tax_bact)
tax_bact$Species_Last_Word <- sapply(strsplit(as.character(tax_bact$Species), " "), tail, 1)
tax_bact$Species2 <- paste(tax_bact$Genus, tax_bact$Species_Last_Word)
tax_bact <- as.matrix(tax_bact)
tax_bact <- tax_bact[, -c(7,8)]
colnames(tax_bact) <- gsub("Species2", "Species", colnames(tax_bact))
tax_bact <- as.matrix(tax_bact)

# Build a phylogenetic tree if required
# This step is computationally intensive and may take a long time
library(DECIPHER)

seqs <- getSequences(seqtab.nochim_bact)
names(seqs) <- seqs

alignment <- AlignSeqs(DNAStringSet(seqs),
                       anchor = NA,
                       verbose = FALSE)

library(phangorn)

phanAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phanAlign)
treeNJ <- NJ(dm)

fit <- pml(treeNJ, data = phanAlign)
fitGTR <- update(fit, k = 4, inv = 0.2)

fitGTR <- optim.pml(fitGTR,
                    model = "GTR",
                    optInv = TRUE,
                    optGamma = TRUE,
                    rearrangement = "stochastic",
                    control = pml.control(trace = 0))

## MODIFY FILES TO ENSURE MATCHING SAMPLE NAMES
## BEFORE CREATING THE PHYLOSEQ OBJECT

# Convert the OTU table into a data frame for sorting
otu_df <- as.data.frame(otu_bact)

# Sort the OTU table by sample names
otu_df <- otu_df[order(rownames(otu_df)), ]

# Sort the metadata using the same sample order
meta_sorted <- meta[order(rownames(meta)), ]

# Convert the metadata back to a sample_data object
meta_sorted <- sample_data(meta_sorted)

# Verify that sample names are correctly aligned
all(rownames(otu_df) == rownames(meta_sorted))

# Remove the ".gz" suffix from sample names
rownames(otu_df) <- gsub("\\.gz$", "", rownames(otu_df))

# Convert the sorted OTU table back into a phyloseq otu_table
otu_bact_sorted <- otu_table(as.matrix(otu_df),
                             taxa_are_rows = FALSE)

# Create the phyloseq object
phy_bact <- phyloseq(
  otu_table(otu_bact_sorted, taxa_are_rows = FALSE),
  tax_table(tax_bact),
  sample_data(meta_sorted)
)

# Use this alternative command if required
# phy_bact <- phyloseq(
#   otu_table(otu_bact_sorted, taxa_are_rows = FALSE),
#   tax_table(tax_bact),
#   meta_sorted
# )

# Replace missing taxonomic assignments with "Unidentified"
phy_bact@tax_table@.Data[
  is.na(phy_bact@tax_table@.Data)
] <- "Unidentified"

# Remove eukaryotic DNA (mitochondria and chloroplast sequences)
phy_bact <- subset_taxa(
  phy_bact,
  Family != "Mitochondria" &
    Order != "Chloroplast"
)

# Remove taxa that could not be identified
phy_bact <- subset_taxa(phy_bact, Species != "Unidentified")
phy_bact <- subset_taxa(phy_bact, Genus != "Unidentified")
phy_bact <- subset_taxa(phy_bact, Family != "Unidentified")
phy_bact <- subset_taxa(phy_bact, Order != "Unidentified")
phy_bact <- subset_taxa(phy_bact, Class != "Unidentified")

# Since one phyloseq object was generated for each project,
# they must be loaded and merged into a single object
phy_bact <- saveRDS(phy_bact, file = "phy_bact1")  # Save each phyloseq object

phy_bact_rr <- readRDS("phy_bact1")                # Load each phyloseq object

phy_bact_combined <- merge_phyloseq(
  phy_bact,
  phy_bact2
)                                                # Merge phyloseq objects

# Convert read counts into relative abundances
phyt_bact <- transform_sample_counts(
  phy_bact,
  function(x) x / sum(x)
)
# phy_bact contains raw read counts,
# whereas phyt_bact contains relative abundances

##FILTER ASVs TO REDUCE NOISE##

install.packages("remotes")
remotes::install_github("Russel88/MicEco")
library("MicEco")  # Required package for ASV filtering

phy_bact <- subset_taxa(phy_bact, Species != "NA NA")

prueba <- ps_prune(phy_bact, min.reads = 10, min.samples = 3)

# Reduce the number of ASVs by retaining only those with
# at least 10 reads and present in at least 3 samples.
# This filtering step helps reduce the false discovery rate (FDR)
# by decreasing the number of statistical tests performed.

# The resulting phyloseq object can be used for downstream
# biodiversity and statistical analyses.
# Before proceeding, verify that no taxa have zero counts
# across all samples to avoid potential errors.

any(taxa_sums(phy_bact) == 0)

phy_bact.2 <- prune_taxa(taxa_sums(phy_bact) > 0, phy_bact)

any(taxa_sums(phy_bact.2) == 0)

####PHYLOSEQ ANALYSIS####

####ALPHA DIVERSITY####

# Calculate alpha diversity indices.
# Specify the diversity metrics of interest in the
# 'measures' argument.
alfa_div <- estimate_richness(
  phy_bact,
  measures = c("Shannon", "Simpson", "Observed")
)

# Save the alpha diversity results
write.table(
  alfa_div,
  "alpha.txt",
  row.names = TRUE,
  sep = "\t"
)

##PLOTS##

{

  # Select the diversity indices to display.
  # If the 'measures' argument is omitted,
  # all eight default indices will be calculated.
  # Update the grouping variable as required.

  p <- plot_richness(
    phy_bact,
    x = "Age_gender",
    color = "Age_gender",
    measures = c("Observed", "Shannon", "Simpson")
  )

  p + geom_point(size = 5, alpha = 0.7)

  p$layers <- p$layers[-1]

  p +
    geom_point(size = 5, alpha = 0) +
    geom_boxplot(
      aes(
        x = Age_gender,
        color = Age_gender,
        fill = Age_gender
      ),
      size = 1
    ) +
    scale_fill_manual(
      values = c(
        "thistle3",
        "steelblue3",
        "steelblue1",
        "thistle1",
        "thistle3",
        "steelblue3",
        "steelblue1",
        "thistle1"
      ),
      name = NULL
    ) +
    scale_color_manual(
      values = c(
        "thistle4",
        "steelblue4",
        "steelblue3",
        "thistle2",
        "thistle4",
        "steelblue4",
        "steelblue3",
        "thistle2"
      ),
      name = NULL
    ) +
    ggtitle("Alpha diversity") +
    labs(
      x = "NA",
      y = "Alpha diversity measures"
    ) +
    scale_x_discrete(
      labels = c(
        "ALZ" = "MCI",
        "Control" = "Control"
      )
    ) +  # Modify x-axis labels if required
    theme_classic() +
    scale_alpha_manual(values = c(1, 0.5)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        size = 8,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 10,
        color = "black"
      ),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),
      legend.key.width = unit(1, "cm"),
      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      )
    )

}

##STATISTICAL ANALYSIS##

# Execute the following only if an unwanted character
# (e.g., "X") appears at the beginning of sample names
# in the alpha diversity table.
row.names(alfa_div) <- gsub("X", "", row.names(alfa_div))

# Merge the alpha diversity table with the metadata
# to enable group comparisons.
sampledf <- data.frame(phy_bact@sam_data)

alfa_div$Sample_ID <- rownames(alfa_div)

# Sort the alpha diversity table by sample ID
alfa_div <- alfa_div[
  order(alfa_div$Sample_ID),
]

# Sort the metadata by sample ID
sampledf <- sampledf[
  order(sampledf$Sample_ID),
]

# Standardize the Sample_ID format in the metadata
sampledf$Sample_ID <- gsub(
  "-",
  ".",
  sampledf$Sample_ID
)

# Merge the alpha diversity table with the metadata
alfa_div2 <- merge(
  alfa_div,
  sampledf,
  by = "Sample_ID"
)

# Perform a Wilcoxon–Mann–Whitney test
# for independent non-parametric samples.
# Modify both the diversity index and the
# grouping variable according to your analysis.
resultado_prueba <- wilcox.test(
  Observed ~ Age_gender,
  data = alfa_div2,
  paired = FALSE
)

resultado_prueba

####BETA DIVERSITY####

# PCoA analysis.
# A phylogenetic tree is required to calculate UniFrac distances.
# If no phylogenetic tree is available, use NMDS with Bray–Curtis distance instead.
ord_PCoA_unifrac = ordinate(phy_bact, "PCoA", "unifrac", weighted = TRUE)

# NMDS analysis.
# This method was used because no phylogenetic tree was generated.
beta_dist_bray <- ordinate(phy_bact, method = "NMDS", distance = "bray")

# Remove empty samples (rows)
phy_bact <- prune_samples(sample_sums(phy_bact) > 0, phy_bact)

# Convert the OTU table into an abundance matrix
abund_table <- otu_table(phy_bact)

# Remove samples containing missing values
phy_bact <- prune_samples(!is.na(sample_sums(phy_bact)), phy_bact)

###PLOT###
# Remember to specify the appropriate ordination method.

# Generate the ordination plot.
# The ordination method (e.g. NMDS, PCoA) and distance metric
# (e.g. Bray, UniFrac, etc.) should be selected as appropriate.
library(ggrepel)
install.packages("ggdist")
library(ggdist)

b <- plot_ordination(phy_bact, beta_dist_bray,
                     shape = "Gender",
                     color = "Gender") +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(1, "cm")
  ) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("thistle3",
                                "steelblue3",
                                "steelblue1",
                                "steelblue4")) +
  theme(
    legend.title = element_blank(),   # Remove the legend title
    plot.title = element_text(
      hjust = 0.5,
      size = 12,
      face = "bold"
    ),
    plot.title.position = "plot"
  ) +
  # Add 95% confidence ellipses around each group
  stat_ellipse(
    aes(group = Gender, color = Gender),
    type = "norm",
    level = 0.95,
    geom = "polygon",
    linetype = 1,
    fill = NA
  )

# Add the plot title
b <- b + ggtitle("Beta diversity")

# Display the plot
print(b)

# Count the number of samples displayed in the ordination plot
ordination_data <- plot_ordination(
  phy_bact,
  beta_dist_bray,
  shape = "Age_gender",
  color = "Age_gender",
  justDF = TRUE
)

num_points <- nrow(ordination_data)

print(paste("Number of samples in the plot:", num_points))

# Generate an ordination plot with sample labels.
# This is useful for identifying potential outliers.
b <- plot_ordination(
  phy_bact,
  ord_PCoA_unifrac,
  shape = "Age_gender",
  color = "Age_gender",
  label = "Sample_ID"
) +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(1, "cm")
  ) +
  geom_point(size = 2) +
  scale_color_manual(values = c("steelblue3",
                                "plum2",
                                "blue")) +
  theme(
    legend.title = element_blank(),   # Remove the legend title
    plot.title = element_text(
      hjust = 0.5,
      size = 12,
      face = "bold"
    ),
    plot.title.position = "plot"
  )

# Add the plot title
b <- b + ggtitle("Beta diversity")

# Display the plot
print(b)

# PERMANOVA is used to determine whether the overall microbial
# community composition differs significantly between groups.

phy_bactotutable <- phy_bact@otu_table
metadf <- data.frame(sample_data(phy_bact))

beta.dist <- vegdist(phy_bactotutable, method = "bray")

set.seed(36)

beta.dist <- adonis2(
  beta.dist ~ Age_range,
  data = metadf
)

# Display the PERMANOVA results,
# including the p-value and R².
beta.dist

# Check for missing values in the distance matrix if necessary
# sum(is.na(beta.dist))

################ TAXON-BY-TAXON ANALYSES ################

### Generate abundance tables at different taxonomic levels ###

{
  # Aggregate taxa at the desired taxonomic rank.
  # Select the taxonomic level using the 'taxrank' argument.
  
  Géneros <- tax_glom(phyt_bact, taxrank = "Genus", NArm = FALSE)
  Familia <- tax_glom(phyt_bact, taxrank = "Family", NArm = FALSE)
  Filo <- tax_glom(phyt_bact, taxrank = "Phylum", NArm = FALSE)
  Especie <- tax_glom(phyt_bact, taxrank = "Species", NArm = FALSE)
  ASVs <- tax_glom(phy_bact, taxrank = "ASV", NArm = FALSE)

  # Remove unidentified taxa
  Géneros <- subset_taxa(Géneros, Genus != "Unidentified")
  Familia <- subset_taxa(Familia, Family != "Unidentified")
  Filo <- subset_taxa(Filo, Phylum != "Unidentified")
  Especie <- subset_taxa(Especie, Species != "Unidentified")
  Especie <- subset_taxa(Especie, Species != "NA NA")
  ASVs <- subset_taxa(ASVs, ASV != "Unidentified")

  # Extract abundance tables and assign taxonomic names
  # as column names.
  # The column index corresponds to the selected taxonomic rank
  # in tax_bact (e.g. column 6 = Genus).

  otu_bact1 <- data.frame(Géneros@otu_table@.Data)
  colnames(otu_bact1) <- Géneros@tax_table[,6]

  otu_bact2 <- data.frame(Familia@otu_table@.Data)
  colnames(otu_bact2) <- Familia@tax_table[,5]

  otu_bact3 <- data.frame(Filo@otu_table@.Data)
  colnames(otu_bact3) <- Filo@tax_table[,2]

  otu_bact4 <- data.frame(Especie@otu_table@.Data)
  colnames(otu_bact4) <- Especie@tax_table[,7]

  otu_bact5 <- data.frame(ASVs@otu_table@.Data)
  colnames(otu_bact5) <- ASVs@tax_table[,8]

  # Merge abundance tables with metadata

  otu_bact1$Sample_ID <- rownames(otu_bact1)
  otu_bact1 <- merge(otu_bact1, meta, by = "Sample_ID")

  # Alternative merging approach
  otu_bact1 <- merge(otu_bact1, meta, by = 0)

  meta <- meta[order(meta$Sample_ID),]

  otu_bact2$Sample_ID <- rownames(otu_bact2)
  # otu_bact2 <- merge(otu_bact2, sampledf, by = "Sample_ID")
  otu_bact2 <- merge(otu_bact2, meta, by = 0)

  otu_bact3$Sample_ID <- rownames(otu_bact3)
  # otu_bact3 <- merge(otu_bact3, sampledf, by = "Sample_ID")
  otu_bact3 <- merge(otu_bact3, meta, by = 0)

  otu_bact4$Sample_ID <- rownames(otu_bact4)
  # otu_bact4 <- merge(otu_bact4, sampledf, by = "Sample_ID")
  otu_bact4 <- merge(otu_bact4, meta, by = 0)

  otu_bact5$Sample_ID <- rownames(otu_bact5)
  # otu_bact5 <- merge(otu_bact5, sampledf, by = "Sample_ID")
  otu_bact5 <- merge(otu_bact5, meta, by = 0)

  # Save abundance tables for each taxonomic level

  write.table(
    otu_bact4,
    "otu_bact4.txt",
    row.names = TRUE,
    sep = "\t"
  )
}

########## Statistical analysis using all taxonomic levels ##########
########## (without generating separate abundance tables beforehand) ##########

# The following methods should be applied to raw read counts,
# not relative abundances, since each method performs its own
# normalization internally.

# IMPORTANT

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiomeMarker")

library(microbiomeMarker)

mm_lefse <- run_lefse(
  phy_bact,
  group = "Age_gender",
  taxa_rank = "Species",
  kw_cutoff = 0.05
)

mm_lefse

plot <- plot_ef_bar(mm_lefse) +
  scale_fill_manual(values = c(
    "thistle3",
    "steelblue3",
    "steelblue1",
    "thistle1"
  )) +
  ggtitle("LEfSe")

plot

write.table(
  mm_lefse@marker_table,
  paste0("lefse.test.", Y, ".txt"),
  sep = "\t"
)

library(dplyr)

# Convert the phyloseq object into a long-format data frame
physeq_df <- psmelt(prueba)

# Remove taxa with missing phylum assignments
physeq_df <- physeq_df %>%
  filter(!is.na(Phylum))

# Identify the 10 most abundant phyla
top_genera <- physeq_df %>%
  group_by(Phylum) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Keep only the top 10 most abundant phyla
physeq_top_genera_df <- physeq_df %>%
  filter(Phylum %in% top_genera) %>%
  group_by(Age_range, Phylum) %>%
  summarise(
    total_abundance = sum(Abundance),
    .groups = "drop"
  )

# Calculate the total abundance within each group
# to obtain relative abundances
physeq_total_by_stage_df <- physeq_top_genera_df %>%
  group_by(Age_range) %>%
  summarise(
    total_abundance_stage = sum(total_abundance),
    .groups = "drop"
  )

# Merge the data frames and calculate percentage abundances
physeq_top_genera_df <- physeq_top_genera_df %>%
  left_join(
    physeq_total_by_stage_df,
    by = "Age_range"
  ) %>%
  mutate(
    percentage_abundance =
      (total_abundance / total_abundance_stage) * 100
  )

# Define a custom color palette
stronger_pastel_colors <- c(
  "#FF9999",
  "#66B2FF",
  "#99FF99",
  "#FFCC99",
  "#FF99FF",
  "#FFDD99",
  "#9999FF",
  "#FF6666",
  "#66FF66",
  "#6699FF"
)

# Generate a stacked bar plot showing the relative abundance
# of the 10 most abundant phyla
ggplot(
  physeq_top_genera_df,
  aes(
    x = Age_range,
    y = percentage_abundance,
    fill = Phylum
  )
) +
  geom_bar(
    stat = "identity",
    position = "stack"
  ) +
  theme_classic() +
  labs(
    x = "Group",
    y = "Relative abundance (%)",
    title = "Relative abundance of the 10 most abundant phyla"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  scale_fill_manual(values = stronger_pastel_colors)




