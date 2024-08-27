# Individual box plots for each genus in a genera list
# Using ggsignif package to calculate significance
##############################################
# Subset the phyloseq object to only include one group
######################################
# Calculate the top genera first

library(phyloseq)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(RColorBrewer)


# Load the data
load("ps2_gg.RData")

# Subset the phyloseq object to include only samples where mouse_microbiome is OMM12 and group is "Salmonella" or "PBS"
ps2_subset <- subset_samples(ps2, 
                             (mouse_microbiome == "OMM12" & group == "PBS" & sample_day == 3) | 
                               (mouse_microbiome == "OMM12" & group == "Salmonella" & sample_day == 4))

# Convert to relative abundance
ps2_subset <- phyloseq::transform_sample_counts(ps2_subset, function(x) x / sum(x))

# Agglomerate to genus-level and rename
ps_genus_subset <- tax_glom(ps2_subset, "Genus")
taxa_names(ps_genus_subset) <- tax_table(ps_genus_subset)[, "Genus"]

# Melt the phyloseq object to long format
df_long_subset <- psmelt(ps_genus_subset)

# Filter specifically for day 4 of Salmonella for top genera calculation
df_Salmonella_day4 <- df_long_subset %>%
  filter(sample_day == "4" & group == "Salmonella")

# Find the most abundant genera based on total abundance in day 4 of Salmonella
top_genera_df <- df_Salmonella_day4 %>%
  group_by(Genus) %>%
  dplyr::summarise(Total_Abundance = sum(Abundance), .groups = "drop") %>%
  top_n(9, Total_Abundance) # Change 12 to the number of top genera you want

# Extract the top genera names
top_genera <- top_genera_df$Genus

# Filter the overall melted data for the top genera from day 4 of Salmonella
df_top_genera_subset <- df_long_subset %>%
  filter(Genus %in% top_genera)

# Convert 'sample_day' and 'group' to factors if they're numeric
df_top_genera_subset$sample_day <- as.factor(df_top_genera_subset$sample_day)
df_top_genera_subset$group <- as.factor(df_top_genera_subset$group)

# Create a new variable combining group and day
df_top_genera_subset$group_day <- paste(df_top_genera_subset$group, df_top_genera_subset$sample_day, sep = " (")

# Close parenthesis
df_top_genera_subset$group_day <- paste0(df_top_genera_subset$group_day, ")")

# Ensure that 'PBS' comes before 'Salmonella' by setting factor levels explicitly
df_top_genera_subset$group_day <- factor(df_top_genera_subset$group_day, 
                                         levels = c("PBS (3)", "Salmonella (4)"))

# Create the box plots for each genus, faceted by group
boxplot_genus <- ggplot(df_top_genera_subset, aes(x = group_day, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  labs(title = "Relative Abundance Comparison Between Day 3\nof PBS and Day 4 of Salmonella in OMM12 Mice",
       x = "", y = "Relative Abundance") +
  facet_wrap(~Genus, scales = "free_y") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 8))

# Specify group comparisons for the significance test
group_comparisons <- list(c("PBS (3)", "Salmonella (4)"))

# Add significance annotations for the specified groups
boxplot_genus <- boxplot_genus +
  geom_signif(
    comparisons = group_comparisons,
    map_signif_level = TRUE,
    test = "wilcox.test",
    test.args = list(p.adjust.method = "BH", paired = FALSE),
    step_increase = 0.1,
    textsize = 4.5,
    family = "serif",
    vjust = -0.2
  )

# Print the plot
print(boxplot_genus)

# Saving the plot
ggsave("plots/omm12_pbs_Salmonella_abun_comparison.png", plot = boxplot_genus, width = 8, height = 6, dpi = 500)

###################################################################
###########################################################


