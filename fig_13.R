# Using unfiltered, untransformed data
load(file = "ps0_rdp_19.Rdata")

ps2 <- ps0
###################
# Alpha Diversity #
###################

library(phyloseq)
library(plyr)
library(ggplot2)

# Calculate alpha diversity measures
alpha_div <- estimate_richness(ps2, measures = c("Shannon", "Chao1"))

# Extract sample data and add a sample ID column
sample_data <- data.frame(sample_data(ps2))
sample_data$SampleID <- rownames(sample_data)

# Add a sample ID column to the alpha_div data frame
alpha_div$SampleID <- rownames(alpha_div)

# Merge sample data with alpha diversity measures
merged_data <- join(sample_data, alpha_div, by = "SampleID")

# Convert Day variable to factor for plotting
merged_data$sample_day <- as.factor(merged_data$sample_day)

# Reordering the levels of the 'group' factor to have 'PBS' first
merged_data$group <- factor(merged_data$group, levels = c("PBS", unique(merged_data$group)[unique(merged_data$group) != "PBS"]))

###########################################

# Convert mouse_microbiome to a factor if it isn't already
merged_data$mouse_microbiome <- as.factor(merged_data$mouse_microbiome)

# Plot
shannon_plot <- ggplot(merged_data, aes(x = sample_day, y = Shannon, fill = group)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), stroke = 1, alpha = 0.6) +
  geom_boxplot(aes(fill = group), position = position_dodge(width = 0.8), alpha = 0.6) +
  geom_smooth(aes(group = sample_day), method = 'loess', se = TRUE, color = "black") +
  labs(x = 'Days after infection', y = 'Alpha Diversity (Shannon Index)', fill = "Group:") +
  theme_classic() +
  theme(legend.position="top",
        axis.title.x=element_text(vjust=0),
        axis.text.x=element_text(size=8, color="black", hjust=0.5),
        axis.line=element_line(colour ="grey"),
        panel.background = element_rect(fill = "transparent", colour = "grey", size = 1, linetype = 1),
        panel.border=element_rect(fill="transparent", color="grey")) +
  facet_grid(mouse_microbiome ~ group, scales = "free_x", space = "free_x")

print(shannon_plot)


# Modified Chao1 Index Plot with updated legend titles
chao1_plot <- ggplot(merged_data, aes(x = sample_day, y = Chao1, fill = group)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), stroke = 1, alpha = 0.6) +
  geom_boxplot(aes(fill = group), position = position_dodge(width = 0.8), alpha = 0.6) +
  geom_smooth(aes(group = sample_day), method = 'loess', se = TRUE, color = "black") +
  labs(x = 'Days after infection', y = 'Alpha Diversity (Chao1 Index)', fill = "Group:") +
  theme_classic() +
  theme(legend.position="top",
        axis.title.x=element_text(vjust=0),
        axis.text.x=element_text(size=8, color="black", hjust=0.5),
        axis.line=element_line(colour ="grey"),
        panel.background = element_rect(fill = "transparent", colour = "grey", size = 1, linetype = 1),
        panel.border=element_rect(fill="transparent", color="grey")) +
  facet_grid(mouse_microbiome ~ group, scales = "free_x", space = "free_x")

print(chao1_plot)

# Save the plot
ggsave("plots/shannon_plot.jpg", plot = shannon_plot, width = 8, height = 5, dpi = 500)
ggsave("plots/chao1_plot.jpg", plot = chao1_plot, width = 8, height = 5, dpi = 500)

##########################################################

# Specify day comparisons for the significance test (this will need to be updated based on your specific analysis)
day_comparisons <- list(c("3", "6"))

# Add significance annotations for the specified days
boxplot_genus <- boxplot_genus +
  geom_signif(
    aes(x = sample_day, group = sample_day),
    comparisons = day_comparisons,
    map_signif_level = TRUE,
    test = "wilcox.test",
    test.args = list(p.adjust.method = "BH", paired = FALSE),
    step_increase = 0.1,
    textsize = 4,
    vjust = -0.2
  )

# Print the plot
print(boxplot_genus)
