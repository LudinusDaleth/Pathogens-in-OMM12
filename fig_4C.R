# Import the required libraries
library(tidyverse)
library(sjmisc)
library(ggplot2)
library(forcats)
library(phyloseq)

# Load the phyloseq object containing the data (ps3 = rdp)
load(file = "ps2_rdp_19.Rdata")

# Subset the phyloseq object to include only samples where mouse_microbiome is OMM12
ps2 <- subset_samples(ps2, mouse_microbiome == "OMM12")

# Convert to relative abundance
ps2 = phyloseq::transform_sample_counts(ps2, function(x){x / sum(x)})

# Extract the required components
otu_df <- otu_table(ps2)

# Check if taxa are columns, if so transpose
if (taxa_are_rows(ps2) == TRUE) {
  otu_df <- t(otu_df) 
}

otu_df <- otu_df %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

tax_df <- tax_table(ps2) %>%
  as.data.frame() %>%
  rownames_to_column("OTU")

sam_df <- sample_data(ps2) %>%
  data.frame()

if (!"sample_id" %in% names(sam_df)) {
  sam_df <- sam_df %>%
    rownames_to_column("sample_id")
}

long_df <- otu_df %>%
  gather(OTU, Abundance, -sample_id) %>%
  left_join(tax_df, by = "OTU") %>%
  left_join(sam_df, by = "sample_id")

# Summarize data by sample_id and Genus
genus_metadata <- long_df %>%
  group_by(sample_id, Genus) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")

# Join with sample metadata
genus_metadata_all <- genus_metadata %>%
  left_join(sam_df, by = "sample_id")

# Filter out rows where Genus is NA
genus_metadata_filtered <- genus_metadata %>%
  filter(!is.na(Genus))

# Identify the top Genera based on their total abundance
top_genera <- genus_metadata_filtered %>%
  group_by(Genus) %>%
  dplyr::summarise(Total_Abundance = sum(Abundance), .groups = "drop") %>%
  top_n(11, Total_Abundance) %>%
  pull(Genus)

# Handle NA Genus names by replacing with higher taxa levels
long_df <- long_df %>%
  mutate(Taxonomy = case_when(
    !is.na(Genus) ~ Genus,
    is.na(Genus) & !is.na(Order) ~ paste("Unclassified (", Order, ")", sep = ""),
    is.na(Genus) & is.na(Order) & !is.na(Order) ~ paste("Unclassified (", Class, ")", sep = ""),
    TRUE ~ "Unclassified"
  ))

# Include specific genera like Listeria
specific_genera <- c("")
top_genera <- unique(c(top_genera, specific_genera))

# Add "Others" category for Genera not in the top_genera list
long_df <- long_df %>%
  mutate(Taxonomy = if_else(Taxonomy %in% top_genera, Taxonomy, "Others"))

# Summarize data by sample_id and Taxonomy
final_df <- long_df %>%
  group_by(sample_id, Taxonomy) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")

# Join with sample metadata
final_df_all <- final_df %>%
  left_join(sam_df, by = "sample_id")

# Print the final dataframe
print(final_df_all)

# Create Genus_abdc dataframe
Genus_abdc <- final_df %>% dplyr::rename(Genus = Taxonomy)

# Prepare data for plotting
data <- final_df_all %>%
  select(sample_id, sample_day, mouse_microbiome, Taxonomy, Abundance, group) %>%
  mutate(sample_day = as.character(sample_day)) %>%
  group_by(sample_day, mouse_microbiome, Taxonomy, group) %>%
  dplyr::summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  ungroup() %>%
  mutate(
    sample_day = as.character(sample_day),
    Taxonomy = as.character(Taxonomy),
    Taxonomy = factor(Taxonomy, levels = rev(c(top_genera, "Others"))),
    sample_day = fct_inorder(sample_day),
    group = fct_relevel(group, "PBS", "Salmonella", "Citrobacter", "Campylobacter", "Listeria"),
  )

# rename Taxonomy to Genus
data <- data %>% dplyr::rename(Genus = Taxonomy)


colors = c("#57375D","#FF3FA4","#9400FF","#8DD3C7","#4daf4a",
           "#377eb8","#BEBADA","#FB8072","#FDB462","#BC41A4",
           "#BC80BD","#FCCDE5","#FFED6F","#CD4F39","#CCEBC5",
           "#4F94CD","#E41A1C","#00CD66")

mycolors <- c("#D9D9D9",rev(colors[1:length(unique(data$Genus))-1]))

#######################
# Create a named vector with only the Genus names you want to change
custom_labels <- c("Akkermansiaceae(Verrucomicrobiales)" = "Akkermansiaceae",
)

####################################

# Reorder Genus based on Abundance
data$Genus <- reorder(data$Genus, data$Abundance, FUN=median)  # using median, but you can change it

# Reorder Genus based on Abundance, but ensure "Others" is always on top
data$Genus <- reorder(data$Genus, data$Abundance, FUN=median)
data$Genus <- factor(data$Genus, levels = c("Others", setdiff(levels(data$Genus), "Others")))

# Define colors ensuring that "Others" is gray
mycolors <- c("#D9D9D9",rev(colors[1:length(unique(data$Genus))-1]))

# Plot
p1 <- ggplot(data, aes(x=sample_day, y=Abundance, fill=Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolors, labels = custom_labels) +
  theme(legend.background = element_rect(size = 0.1),
        legend.key.size= unit(5,"mm"),
        legend.text = element_text(size = 8)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        axis.title.x=element_text(vjust=0),  # Adjust this value
        axis.text.x=element_text(size=8, color="black", hjust=0.5),  # Adjust hjust here
        axis.line=element_line(colour ="grey"),
        panel.background = element_rect(fill = "transparent",colour = "grey", size = 1, linetype = 1),
        panel.border=element_rect(fill="transparent",color="grey")) +  # The same theme you've provided
  facet_grid(mouse_microbiome ~ group, scales = "free_x", space = "free_x") +
  ylab("Relative abundance") +
  xlab("Time (Days)")

p1

# Saving the plot
ggsave("plots/combined_genus_avg_omm12.png", width = 8, height = 4, dpi = 500)

