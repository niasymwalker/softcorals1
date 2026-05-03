library(ggplot2)
library(tidyverse)
library(fs)
library(readxl)
library(dplyr)
library(cowplot)
library(readr)
library(vegan)
library(forcats)
library(rstatix)
library(dunn.test)
library(ggtext)

#####Universal plotting variables#####
# Genus fill colors — colorblind-safe 
genus_colors <- c(
  "Green Star Polyp" = "#785EF0",   # purple
  "Pulse Coral"            = "#DC267F",   # pinkish-red
  "Leather Coral"    = "#FE6100"    # orange
)

# Store colors for x-axis labels
store_colors <- c(
  "Tropical Reef Store" = "#EE6677",   # rose/pink
  "Jan's"               = "#CCBB44",   # yellow-olive
  "Arts Underwater"     = "#4477AA",   # blue
  "Bob's Tropical Reef" = "#228833",   # green  
  "CK Fish World"       = "#AA3377",   # purple-magenta
  "Pisces Coral"        = "#66CCEE"    # cyan
)
# Color palette 
my_colors_profile <- c(
  # Cladocopium (C) Tol's sequential blue, dark to light
  "C40c"      = "#002288",
  "C1"        = "#0044aa",
  "C3bt"      = "#1166cc",
  "C71"       = "#3388ee",
  "C107b"     = "#66afff",
  "C1bp"      = "#99caff",
  
  # Durusdinium (D)  Tol's sequential orange-red, dark to light
  "D5"        = "#550000",
  "D4"        = "#880000",
  "D3a"       = "#aa1100",
  "D9b"       = "#cc4411",
  "D17m"      = "#ee7733",
  "D2"        = "#ffaa77",
  "D17j"      = "#ffd5b8",
  
  # Other groups  muted neutrals, distinct from both blue and orange
  "Other (C)" = "#77aadd",   # muted blue-grey — matches C family
  "Other (D)" = "#ddaa77",   # muted tan-orange — matches D family
  "Other (B)" = "#44bb99",   # teal — Tol's safe green, distinct from both
  "Other (A)" = "#99ddbb",   # light teal
  "Other (G)" = "#bbcc33",   # yellow-green — Tol safe, distinct from red/blue
  "Other (F)" = "#aaaa00",   # olive
  "Other (I)" = "#dddddd",   # light grey
  "Other"     = "#bbbbbb"    # mid grey
)

#Store and genus order
store_order <- c("Tropical Reef Store", "Jan's", "Arts Underwater",
                 "Bob's Tropical Reef", "CK Fish World", "Pisces Coral")

genus_order <- c("Green Star Polyp", "Pulse Coral", "Leather Coral")

#####Fig 1 plus data mining#####
its2_data <- read.csv("FINAL_ITS2data.csv", stringsAsFactors = FALSE)

table(its2_data$Experiment_Type)

ITS2data <- its2_data %>%
  filter(Experiment_Type == "Store")

ITS2heatstress <- its2_data %>%
  filter(Experiment_Type %in% c("Heated", "Control"))

# DIV columns start at column 6 onward
div_cols <- names(its2_data)[6:ncol(its2_data)]

ITS2heatstress_long <- ITS2heatstress %>%
  pivot_longer(cols      = all_of(div_cols),
               names_to  = "DIV",
               values_to = "reads") %>%
  filter(reads > 0)   

ITS2data_long <- ITS2data %>%
  pivot_longer(cols      = all_of(div_cols),
               names_to  = "DIV",
               values_to = "reads") %>%
  filter(reads > 0)

profiles <- ITS2data %>%
  mutate(
    Sample = as.factor(Sample),
    sum =rowSums(.[, 6:763], na.rm = TRUE)  # More explicit column selection
  ) %>%
  pivot_longer(
    cols = 6:763,
    names_to = "type_profile",
    values_to = "current_reads"
  ) %>%
  mutate(
    type_profile_prop = current_reads / sum,
    type_profile = as.factor(type_profile)
  ) %>%
  filter(current_reads > 0) %>%  # Remove zeros directly instead of converting to NA then dropping
  # Assign clade based on first letter of type profile
  mutate(
    Clades = case_when(
      grepl("^C", type_profile) ~ "Cladocopium",
      grepl("^D", type_profile) ~ "Durusdinium",
      grepl("^A", type_profile) ~ "Symbiodinium",
      grepl("^G", type_profile) ~ "Gerakladium",
      grepl("^I", type_profile) ~ "Clade I",
      grepl("^B", type_profile) ~ "Breviolum",
      grepl("^F", type_profile) ~ "Fugacium",
      TRUE ~ NA_character_
    )
  )

## Summarize ITS2 data by clade per colony
profiles_clades <- profiles %>%
  group_by(Colony, Clades) %>%
  summarise(
    clade_sum = sum(current_reads),
    Genus = dplyr::first(Genus),
    Store = dplyr::first(Store),
    .groups = "drop"
  ) %>%
  group_by(Colony) %>%
  mutate(
    sum_colony = sum(clade_sum),
    clade_prop = round(clade_sum / sum_colony, 5)
  ) %>%
  ungroup() %>%
  dplyr::select(Colony, Genus, Clades, clade_prop, Store)

# Identify the dominant clade for each colony
maxclades <- profiles_clades %>%
  group_by(Colony) %>%
  slice_max(clade_prop, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Colony, Genus, major_clade = Clades, clade_prop)

# Identify which clade has the highest proportion for each colony
dominant_per_colony <- profiles_clades %>%
  group_by(Colony, Genus, Store) %>%
  slice_max(clade_prop, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(dominant_clade      = Clades,
         dominant_clade_prop = clade_prop)

#Summary by Store
store_summary <- dominant_per_colony %>%
  group_by(Store, dominant_clade) %>%
  summarise(
    n_colonies  = n(),
    mean_prop   = round(mean(dominant_clade_prop), 4),
    se_prop     = round(sd(dominant_clade_prop) / sqrt(n()), 4),
    min_prop    = round(min(dominant_clade_prop), 4),
    max_prop    = round(max(dominant_clade_prop), 4),
    .groups     = "drop"
  ) %>%
  mutate(Store = factor(Store, levels = store_order)) %>%
  arrange(Store, desc(n_colonies))

print(store_summary, n = Inf)

#Summary by Genus
genus_summary <- dominant_per_colony %>%
  group_by(Genus, dominant_clade) %>%
  summarise(
    n_colonies  = n(),
    mean_prop   = round(mean(dominant_clade_prop), 4),
    se_prop     = round(sd(dominant_clade_prop) / sqrt(n()), 4),
    min_prop    = round(min(dominant_clade_prop), 4),
    max_prop    = round(max(dominant_clade_prop), 4),
    .groups     = "drop"
  ) %>%
  mutate(Genus = factor(Genus, levels = genus_order)) %>%
  arrange(Genus, desc(n_colonies))

print(genus_summary, n = Inf)

#Summary by Store x Genus
store_genus_summary <- dominant_per_colony %>%
  group_by(Store, Genus, dominant_clade) %>%
  summarise(
    n_colonies  = n(),
    mean_prop   = round(mean(dominant_clade_prop), 4),
    se_prop     = round(sd(dominant_clade_prop) / sqrt(n()), 4),
    min_prop    = round(min(dominant_clade_prop), 4),
    max_prop    = round(max(dominant_clade_prop), 4),
    .groups     = "drop"
  ) %>%
  mutate(
    Store = factor(Store, levels = store_order),
    Genus = factor(Genus, levels = genus_order)
  ) %>%
  arrange(Store, Genus, desc(n_colonies))

print(store_genus_summary, n = Inf)

#Normalize proportions by total reads per colony
profiles <- profiles %>%
  group_by(Colony) %>%
  mutate(
    colony_total_reads = sum(current_reads),
    prop_norm          = current_reads / colony_total_reads
  ) %>%
  ungroup()

#Clade proportions per colony
profiles_clades_wide <- profiles %>%
  group_by(Colony, Clades) %>%
  summarise(clade_prop = sum(prop_norm), .groups = "drop") %>%
  pivot_wider(names_from = Clades, values_from = clade_prop, values_fill = 0)

if (!"Cladocopium" %in% names(profiles_clades_wide))  profiles_clades_wide$Cladocopium  <- 0
if (!"Durusdinium"  %in% names(profiles_clades_wide))  profiles_clades_wide$Durusdinium   <- 0

profiles_clades_wide <- profiles_clades_wide %>%
  rename(Cladocopium_prop = Cladocopium,
         Durusdinium_prop = Durusdinium)

# Collapse rare DIVs into "Other (Clade)"
total_reads_all <- sum(profiles$current_reads)

top_divs <- profiles %>%
  group_by(type_profile) %>%
  summarise(total = sum(current_reads), .groups = "drop") %>%
  filter(total / total_reads_all >= 0.005) %>%
  pull(type_profile)

profiles2 <- profiles %>%
  mutate(
    clade_short = case_when(
      Clades == "Cladocopium"  ~ "C",
      Clades == "Durusdinium"  ~ "D",
      Clades == "Symbiodinium" ~ "A",
      Clades == "Breviolum"    ~ "B",
      Clades == "Gerakladium"  ~ "G",
      Clades == "Fugacium"     ~ "F",
      Clades == "Clade I"      ~ "I",
      TRUE                     ~ "Other"
    ),
    type_profile_plot = if_else(
      type_profile %in% top_divs,
      type_profile,
      paste0("Other (", clade_short, ")")
    )
  ) %>%
  group_by(Colony, Genus, Store, type_profile_plot, Clades) %>%
  summarise(plot_prop = sum(prop_norm), .groups = "drop") %>%
  left_join(profiles_clades_wide %>% select(Colony, Cladocopium_prop, Durusdinium_prop),
            by = "Colony")

colony_order_df <- profiles2 %>%
  distinct(Colony, Genus, Store) %>%
  mutate(
    Store = factor(Store, levels = store_order),
    Genus = factor(Genus, levels = genus_order)
  ) %>%
  arrange(Store, Genus, Colony) %>%
  mutate(Colony = as.character(Colony))

colony_levels <- colony_order_df$Colony
profiles2$Colony <- factor(profiles2$Colony, levels = colony_levels)

# Legend factor levels: named C → named D → all Other groups at bottom
all_plot_levels <- unique(as.character(profiles2$type_profile_plot))

named_c   <- sort(grep("^C", all_plot_levels, value = TRUE))   # e.g. C40c, C1, ...
named_d   <- sort(grep("^D", all_plot_levels, value = TRUE))   # e.g. D17j, D17m, ...
named_b   <- sort(grep("^B", all_plot_levels, value = TRUE))   # Breviolum (if present)
other_grp <- sort(grep("^Other", all_plot_levels, value = TRUE)) # all Other(...) last
remaining  <- setdiff(all_plot_levels, c(named_c, named_d, named_b, other_grp))

profile_levels <- c(named_c, named_d, named_b, remaining, other_grp)
profiles2$type_profile_plot <- factor(profiles2$type_profile_plot, levels = profile_levels)

label_colors <- colony_order_df %>%
  pull(Store) %>%
  as.character() %>%
  { store_colors[.] }

# Store divider positions
store_breaks <- colony_order_df %>%
  mutate(x_pos = row_number()) %>%
  group_by(Store) %>%
  summarise(x_end = max(x_pos), .groups = "drop") %>%
  arrange(x_end)

divider_positions <- store_breaks$x_end[-nrow(store_breaks)] + 0.5

# Build caption with color-coded store key 
store_caption <- paste0(
  "\u25a0 Tropical Reef Store  \u25a0 Arts Underwater  \u25a0 Bob's Tropical Reef  ",
  "\u25a0 CK Fish World  \u25a0 Pisces Coral\n"
)

# Build the actual plot
stacked_barplot <- ggplot(profiles2, aes(x = Colony, y = plot_prop)) +
  geom_col(aes(fill = type_profile_plot), colour = "grey70", linewidth = 0.04) +
  geom_vline(xintercept = divider_positions, color = "grey30",
             linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(values = my_colors_profile, name = "ITS2 Type") +
  scale_y_continuous(
    limits = c(0, 1.001),
    breaks = c(0, 0.33, 0.67, 1),
    labels = c("0", "0.33", "0.67", "1.00")
  ) +
  labs(
    x       = "Colony",
    y       = "Symbiodiniaceae Proportion",
    caption = store_caption
  ) +
  guides(
    fill = guide_legend(ncol = 1, title.theme = element_text(angle = 0, size = 10))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.x    = element_text(size = 13),
    axis.title.y    = element_text(size = 13),
    axis.text.x     = element_text(
      angle  = 90,
      vjust  = 0.5,
      hjust  = 1,
      size   = 11,
      color  = label_colors
    ),
    axis.ticks.x    = element_blank(),
    axis.text.y     = element_text(size = 10),
    text            = element_text(size = 11),
    panel.grid      = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text     = element_text(size = 10),
    plot.title      = element_text(size = 10, face = "bold"),
    plot.caption    = element_text(size = 10, color = "grey30", hjust = 0)
  )

#ggtext for colored caption (install if needed: install.packages("ggtext"))
if (requireNamespace("ggtext", quietly = TRUE)) {
  library(ggtext)
  colored_caption <- paste0(
    "<span style='color:#EE6677'>&#9632;</span> Tropical Reef Store  ",
    "<span style='color:#CCBB44'>&#9632;</span> Jan's  ",
    "<span style='color:#4477AA'>&#9632;</span> Arts Underwater  ",
    "<span style='color:#228833'>&#9632;</span> Bob's Tropical Reef  ",
    "<span style='color:#AA3377'>&#9632;</span> CK Fish World  ",
    "<span style='color:#66CCEE'>&#9632;</span> Pisces Coral<br>"
  )
  F1 <- stacked_barplot +
    labs(caption = colored_caption) +
    theme(plot.caption = ggtext::element_markdown(size = 10, hjust = 0, color = "grey30"))
} else {
  message("Install 'ggtext' for colored caption squares: install.packages('ggtext')")
  # plain text fallback already set above
}

print(F1)

##Fisher's exact test for symbiont proportion by store
store_sym <- matrix(c(7,1,2,0,1,2,1,3,4,2,5,0),ncol=2,nrow=6,dimnames=list(c("TRS", "Jans", "Arts", "BTR", "CKFish","Pisces"),c("Cladocopium","Durusdinium")),byrow=TRUE)
fisher.test(store_sym)
pairwise_fisher_test(as.matrix(store_sym))

#Fisher's exact test for symbiont proportion by genus
genus_sym <- matrix(c(8,2,3,5,9,1),ncol=2,nrow=3,dimnames=list(c("GSP", "PC", "LC"),c("Cladocopium","Durusdinium")),byrow=TRUE)
fisher.test(genus_sym)
pairwise_fisher_test(as.matrix(genus_sym))

quartz(w=7.2, h=4)
plot_grid(F1)
quartz.save("./Fig1_final.pdf", type="pdf")

#####Fig 2A this is the PCA#####
# Build clade proportion matrix
profiles_filtered <- profiles %>%
  filter(Clades %in% c("Cladocopium", "Durusdinium"))

clade_mat <- profiles_filtered %>%
  group_by(Colony, Clades) %>%
  summarise(reads = sum(current_reads), .groups = "drop") %>%
  group_by(Colony) %>%
  mutate(prop = reads / sum(reads)) %>%
  ungroup() %>%
  select(Colony, Clades, prop) %>%
  pivot_wider(names_from = Clades, values_from = prop, values_fill = 0) %>%
  column_to_rownames("Colony")
  
# Metadata (same row order as matrix) 
meta <- profiles_filtered %>%
  distinct(Colony, Store, Genus) %>%
  arrange(match(Colony, rownames(clade_mat)))

stopifnot(all(meta$Colony == rownames(clade_mat)))

# Run PCA 
# Hellinger transformation first — standard for proportion/compositional data
clade_hell <- decostand(clade_mat, method = "hellinger")
pca_result <- rda(clade_hell)

# Re-extract scores using pca_result
site_scores <- scores(pca_result, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Colony") %>%
  left_join(meta, by = "Colony")

species_scores <- scores(pca_result, display = "species") %>%
  as.data.frame() %>%
  rownames_to_column("Clade")

eig     <- eigenvals(pca_result)
pct_var <- round(100 * eig / sum(eig), 1)

# ENVFIT
set.seed(123)
env_fit <- envfit(pca_result,
                  env          = meta %>% select(Store, Genus),
                  permutations = 999,
                  na.rm        = TRUE)

print(env_fit)

# Extract centroids for all levels
env_centroids <- as.data.frame(scores(env_fit, display = "factors")) %>%
  rownames_to_column("Label") %>%
  mutate(
    Variable = case_when(
      Label %in% unique(meta$Store) ~ "Store",
      TRUE ~ "Other"
    ),
    # Pull p-value per level for filtering
    pval = env_fit$factors$pvals[Label]
  )

# Separate store and genus centroids
store_centroids <- env_centroids %>% filter(Variable == "Store")

genus_shapes <- c(
  "Green Star Polyp" = 16,
  "Pulse Coral" = 17,
  "Leather Coral" = 15
)

# Scale factor to make envfit arrows visually comparable to site scores
arrow_scale <- 0.8

F2A <- ggplot(site_scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "grey70", linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(color = Store, shape = Genus), size = 3.5, alpha = 0.9) +
  geom_segment(data = species_scores,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow     = arrow(length = unit(0.25, "cm"), type = "closed"),
               color     = "grey30",
               linewidth = 0.6) +
  geom_text(data = species_scores,
            aes(x = PC1 * 1.15, y = PC2 * 1.15, label = Clade),
            size     = 4,
            color    = "grey20",
            fontface = "italic") +
  geom_segment(data = store_centroids,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale,
                   yend = PC2 * arrow_scale),
               arrow = arrow(length = unit(0.25, "cm"), type = "open"),
               color = "grey50",
               linewidth = 0.5,
               linetype  = "dashed") +
  ggrepel::geom_text_repel(data = store_centroids,
                           aes(x     = PC1 * arrow_scale,
                               y     = PC2 * arrow_scale,
                               label = Label),
                           size        = 3.2,
                           color       = "grey30",
                           fontface    = "bold",
                           box.padding = 0.3) +
  scale_color_manual(values = store_colors, name = "Store") +
  scale_shape_manual(values = genus_shapes, name = "Genus") +
  labs(
    x     = paste0("PC1 (", pct_var[1], "% variance)"),
    y     = paste0("PC2 (", pct_var[2], "% variance)"),
  ) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 11, face = "bold")
  )

print(F2A)

# Average proportion of non-Cladocopium and non-Durusdinium ITS2 types per colony
other_props <- profiles2 %>%
  filter(!Clades %in% c("Cladocopium", "Durusdinium")) %>%
  group_by(Colony) %>%
  summarise(other_prop = sum(plot_prop), .groups = "drop")

# Summary across all colonies
other_summary <- other_props %>%
  summarise(
    mean_other= round(mean(other_prop), 4),
    sd_other  = round(sd(other_prop),   4),
    min_other= round(min(other_prop),  4),
    max_other = round(max(other_prop),  4),
    pct_mean  = scales::percent(mean(other_prop), accuracy = 0.01)
  )

print(other_summary)

#breakdown by genus
profiles2 %>%
  filter(!Clades %in% c("Cladocopium", "Durusdinium")) %>%
  group_by(Colony, Genus) %>%
  summarise(other_prop = sum(plot_prop), .groups = "drop") %>%
  group_by(Genus) %>%
  summarise(
    n = n(),
    mean_other = round(mean(other_prop), 4),
    sd_other  = round(sd(other_prop),  4),
    max_other = round(max(other_prop), 4),
    .groups = "drop"
  )

#####Fig 2B - shannon diversity by store####
# Calculate Shannon diversity per colony
mat <- profiles %>%
  group_by(Colony, type_profile) %>%
  summarise(reads = sum(current_reads), .groups = "drop") %>%
  group_by(Colony) %>%
  mutate(prop = reads / sum(reads)) %>%
  ungroup() %>%
  select(Colony, type_profile, prop) %>%
  pivot_wider(names_from = type_profile, values_from = prop, values_fill = 0) %>%
  column_to_rownames("Colony")

shannon_df <- diversity(mat, index = "shannon") %>%
  enframe(name = "Colony", value = "Shannon")

meta <- profiles %>%
  distinct(Colony, Store, Genus)

shannon_df <- shannon_df %>%
  left_join(meta, by = "Colony") %>%
  mutate(Store = factor(Store, levels = store_order)) %>%
  mutate(Genus = factor(Genus, levels = genus_order))

# Sample size labels
n_labels <- shannon_df %>%
  group_by(Store) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(Store, "\n(n = ", n, ")"))

store_labels <- setNames(n_labels$label, as.character(n_labels$Store))

F2B <- ggplot(shannon_df, aes(x = Store, y = Shannon, fill = Store, color = Store)) +
  geom_boxplot(width = 0.4, alpha = 0.25, linewidth = 0.5,
               outlier.shape = NA, color = "grey20") +
  geom_jitter(aes(shape = Genus), width = 0.1, size = 3, alpha = 0.9) +
  scale_fill_manual(values  = store_colors) +
  scale_color_manual(values = store_colors) +
  scale_shape_manual(values = c("Green Star Polyp" = 16,
                                "Pulse Coral"            = 17,
                                "Leather Coral"    = 15),
                     name = "Genus") +
  scale_x_discrete(labels = store_labels) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = seq(0, 6, 0.5),
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x       = NULL,
    y       = "Shannon Diversity (H')"
  ) +
  guides(
    fill  = "none",
    color = "none",
    shape = guide_legend(override.aes = list(size = 3, color = "grey30"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y       = element_text(size = 13),
    axis.text.x        = element_text(size = 10, angle = 35, hjust = 1,
                                      color = store_colors[store_order]),
    axis.text.y        = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    legend.title       = element_text(size = 11),
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.4, "cm")
  )

print(F2B)

# Summary table
shannon_df %>%
  group_by(Store) %>%
  summarise(
    n       = n(),
    mean_H  = round(mean(Shannon), 3),
    sd_H    = round(sd(Shannon),   3),
    min_H   = round(min(Shannon),  3),
    max_H   = round(max(Shannon),  3),
    .groups = "drop"
  ) %>%
  print()

#Kruskal-Wallis by Store 
kw_store <- kruskal.test(Shannon ~ Store, data = shannon_df)
pairwise.wilcox.test(shannon_df$Shannon, shannon_df$Store, p.adjust.method = "BH")

print(kw_store)

#####Fig 2C shannon diversity by genus#####
# Shannon Alpha Diversity by Genus
# Sample size labels for x-axis
n_labels_genus <- shannon_df %>%
  group_by(Genus) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(Genus, "\n(n = ", n, ")"))

genus_labels <- setNames(n_labels_genus$label, as.character(n_labels_genus$Genus))

F2C <- ggplot(shannon_df,
                            aes(x = Genus, y = Shannon,
                                fill = Genus, color = Genus)) +
  
  geom_boxplot(width = 0.4, alpha = 0.25, linewidth = 0.5,
               outlier.shape = NA, color = "grey20") +
  
  # Points colored by store, shaped by genus
  geom_jitter(aes(color = Store),
              width = 0.1, size = 3, alpha = 0.9) +
  scale_fill_manual(values = genus_colors)  +
  scale_color_manual(values = store_colors, name = "Store") +
  scale_x_discrete(labels = genus_labels) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = seq(0, 6, 0.5),
                     expand = expansion(mult = c(0.02, 0.08))) +
  annotate("segment", x = 1.1, xend = 1.9, y= 2.35, yend=2.35) +
  annotate("text", x=1.5, y=2.4, label="*", size=5) +
  labs(x = NULL,
       y = "Shannon Diversity (H')") +
  guides(
    fill  = "none",
    color = guide_legend(override.aes = list(size = 3, shape = 16))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y       = element_text(size = 13),
    axis.text.x        = element_text(size = 11, angle = 0, hjust = 0.5,
                                      color = genus_colors[genus_order]),
    axis.text.y        = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    legend.title       = element_text(size = 11),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.4, "cm")
  )

print(F2C)

#Kruskal-Wallis by Genus
kw_genus <- kruskal.test(Shannon ~ Genus, data = shannon_df)
pairwise.wilcox.test(shannon_df$Shannon, shannon_df$Genus, p.adjust.method = "hochberg")
print(kw_genus)

# Summary table
shannon_df %>%
  group_by(Genus) %>%
  summarise(
    n      = n(),
    mean_H = round(mean(Shannon), 3),
    sd_H   = round(sd(Shannon),   3),
    min_H  = round(min(Shannon),  3),
    max_H  = round(max(Shannon),  3),
    .groups = "drop"
  ) %>%
  print()

quartz(w=6.5,h=4)
F2A
quartz.save("./Fig2A_final.pdf", type="pdf")
quartz(w=6,h=3)
F2B
quartz.save("./Fig2B_final.pdf", type="pdf")
quartz(w=4,h=3)
F2C
quartz.save("./Fig2C_final.pdf", type="pdf")

quartz(w=11.7,h=6)
F2B <- F2B + ggtitle("B")
F2C <- F2C + ggtitle("C")
F2BC <- plot_grid(F2B, F2C, rel_heights = c(1.2,1), ncol=1)
plot_grid(F2A, F2BC, labels = c("A",""), rel_widths = c(1.2,0.6), nrow=1)
quartz.save("./Fig2_final.pdf", type="pdf")

# Join dominant clade to shannon_df first
shannon_df <- shannon_df %>%
  left_join(
    profiles %>%
      group_by(Colony) %>%
      mutate(total_reads = sum(current_reads)) %>%
      group_by(Colony, Clades) %>%
      summarise(prop = sum(current_reads) / first(total_reads), .groups = "drop") %>%
      group_by(Colony) %>%
      slice_max(prop, n = 1, with_ties = FALSE) %>%
      select(Colony, dominant_clade = Clades),
    by = "Colony"
  )

#---- Split into two groups
clado_shannon <- shannon_df %>% filter(dominant_clade.y == "Cladocopium") %>% pull(Shannon)
durus_shannon  <- shannon_df %>% filter(dominant_clade.y == "Durusdinium")  %>% pull(Shannon)

# ---- Summary stats per group
shannon_df %>%
  group_by(dominant_clade.y) %>%
  summarise(
    n      = n(),
    mean   = round(mean(Shannon), 3),
    median = round(median(Shannon), 3),
    sd     = round(sd(Shannon), 3),
    min    = round(min(Shannon), 3),
    max    = round(max(Shannon), 3)
  )

#Kruskal-Wallis by Store 
kw_store <- kruskal.test(Shannon ~ Store, data = shannon_df)
pairwise.wilcox.test(shannon_df$Shannon, shannon_df$Store, p.adjust.method = "fdr")

print(kw_store)

#Kruskal-Wallis by Genus
kw_genus <- kruskal.test(Shannon ~ Genus, data = shannon_df)
pairwise.wilcox.test(shannon_df$Shannon, shannon_df$Genus, p.adjust.method = "fdr")

print(kw_genus)

# Post-hoc Dunn test (only meaningful if KW is significant 
dunn.test(shannon_df$Shannon, shannon_df$Store,
          method = "bh", kw = TRUE, label = TRUE)


dunn.test(shannon_df$Shannon, shannon_df$Genus,
          method = "bh", kw = TRUE, label = TRUE)

median_genus <- shannon_df %>%
  group_by(Genus) %>%
  summarize(median_value = median(Shannon), na.rm = TRUE) # na.rm = TRUE handles missing values

# View the result
print(median_genus)

#####Fig 3A - thermal stress assay profile + gen analysis####
temp_df <- read_excel('./R_data.xlsx', 2) %>%
  mutate(date_time  = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%S"),
         date       = as.Date(date_time),
         treatment  = str_to_title(treatment))  # "heated" -> "Heated", "control" -> "Control"

#Define objects 
temp_colors    <- c("Heated" = "#CC3311", "Control" = "#4477AA")
temp_linetypes <- c("Heated" = "solid",   "Control" = "solid")

# ---- Plot raw hourly data to show daily ramp cycles
F3A <- ggplot(temp_df,
                       aes(x        = date_time,
                           y        = temperature,
                           color    = treatment,
                           linetype = treatment,
                           group    = treatment)) +
  geom_vline(xintercept = as.POSIXct("2025-06-12 2:30:00"),
             color = "grey40", linewidth = 0.5, linetype = "dotted") +
  geom_line(data = temp_df %>% filter(treatment == "Control"),
            linewidth = 0.9, alpha = 1.0) +
  geom_line(data = temp_df %>% filter(treatment == "Heated"),
            linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(values    = temp_colors,    name = "Treatment") +
  scale_linetype_manual(values = temp_linetypes, name = "Treatment") +
  scale_x_datetime(date_labels = "%b %d", date_breaks = "1 day") +
  scale_y_continuous(limits = c(24, 35),
                     breaks = seq(25, 35, 2)) +
  labs(
    x       = "Date",
    y       = "Temperature (°C)",
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 11),
    axis.text.y      = element_text(size = 11),
    axis.title       = element_text(size = 12),
    panel.grid.minor = element_blank(),
    legend.position = c(0.5, 0.85),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(0.4, "cm"),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.caption     = element_text(size = 10, color = "grey40", hjust = 0)
  ) +
  coord_cartesian(clip = "off")

print(F3A)

#####Fig 3B - heat tolerance Sym stacked bar plot#####
div_cols <- names(ITS2heatstress)[6:ncol(ITS2heatstress)]

hs_long <- ITS2heatstress %>%
  pivot_longer(cols      = all_of(div_cols),
               names_to  = "DIV",
               values_to = "reads") %>%
  filter(reads > 0) %>%
  group_by(Sample) %>%
  mutate(prop = reads / sum(reads)) %>%
  ungroup() %>%
  mutate(Clade = str_extract(DIV, "^[A-Z]"))

# Collapse rare DIVs into "Other (Clade X)" — same approach as store barplot
top_divs <- hs_long %>%
  group_by(DIV) %>%
  summarise(total = sum(reads), .groups = "drop") %>%
  slice_max(total, n = 20) %>%
  pull(DIV)

hs_long <- hs_long %>%
  mutate(DIV_plot = if_else(DIV %in% top_divs, DIV,
                            paste0("Other (", Clade, ")")))

# Ordered factor for x-axis: Colony + Treatment side by side
# Colony order for grouping
colony_order <- c("GSP1","GSP2","GSP3","PC1","PC2","PC3","LC1","LC2","LC3","LC4")

hs_long <- hs_long %>%
  mutate(
    Colony          = factor(Colony, levels = colony_order),
    Experiment_Type = factor(Experiment_Type, levels = c("Control","Heated")),
    Genus           = factor(Genus, levels = genus_order),
    x_label         = paste0(Colony, "\n", Experiment_Type)
  )

# Ordered x_label factor preserving colony + treatment grouping
x_order <- hs_long %>%
  distinct(Colony, Experiment_Type, x_label) %>%
  arrange(Colony, Experiment_Type) %>%
  pull(x_label)

hs_long <- hs_long %>%
  mutate(x_label = factor(x_label, levels = x_order))

# Reuse clade-based palette logic from store barplot
clade_base_colors <- c(
  "C" = "#4477AA",
  "D" = "#EE6677",
  "A" = "#228833",
  "B" = "#CCBB44",
  "F" = "#AA3377",
  "G" = "#66CCEE",
  "H" = "#F1B6DA",
  "I" = "#BBBBBB"
)

div_levels <- hs_long %>%
  group_by(DIV_plot, Clade) %>%
  summarise(total = sum(reads), .groups = "drop") %>%
  arrange(Clade, desc(total))

# Generate shades per clade
make_shades <- function(base_color, n) {
  if (n == 1) return(base_color)
  colorRampPalette(c(base_color, colorspace::lighten(base_color, 0.65)))(n)
}

div_colors <- div_levels %>%
  group_by(Clade) %>%
  mutate(shade = make_shades(clade_base_colors[Clade], n())[row_number()]) %>%
  ungroup() %>%
  select(DIV_plot, shade) %>%
  deframe()

hs_bar <- hs_long %>%
  group_by(Colony, Genus, Experiment_Type, DIV_plot) %>%
  summarise(prop = mean(prop), .groups = "drop") %>%   # mean across replicates
  mutate(
    Colony          = factor(Colony, levels = colony_order),
    Experiment_Type = factor(Experiment_Type, levels = c("Control", "Heated")),
    x_label         = paste0(Colony, "\n", Experiment_Type),
    x_label         = factor(x_label, levels = x_order)
  )

d_proportion_summary <- hs_long %>%
  mutate(is_D = grepl("^D", DIV_plot)) %>%
  group_by(Genus, Colony, Experiment_Type) %>%
  summarise(
    total_reads = sum(reads, na.rm = TRUE),
    d_reads     = sum(reads[is_D], na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(prop_D = d_reads / total_reads) %>%
  select(Genus, Colony, Experiment_Type, prop_D) %>%
  pivot_wider(
    names_from = Experiment_Type, 
    values_from = prop_D
  ) %>%
  mutate(
    perc_diff = (Heated - Control) * 100
  )

print(d_proportion_summary)

mean(d_proportion_summary$perc_diff)
se <- sd(d_proportion_summary$perc_diff) / sqrt(length(d_proportion_summary$perc_diff))

overall_test <- wilcox.test(d_proportion_summary$Heated,d_proportion_summary$Control,paired = TRUE)
overall_test

genus_stats <- d_proportion_summary %>%
  group_by(Genus) %>%
  summarise(
    p_value = wilcox.test(Heated, Control, paired = TRUE)$p.value,
    median_diff = median(perc_diff, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(genus_stats)

genus_summary_stats <- d_proportion_summary %>%
  group_by(Genus) %>%
  summarise(
    mean_diff = mean(perc_diff, na.rm = TRUE),
    sd_diff   = sd(perc_diff, na.rm = TRUE),
    n         = n(),
    # Standard Error calculation:
    se_diff   = sd_diff / sqrt(n),
    .groups   = "drop"
  )

print(genus_summary_stats)

kruskal_result <- kruskal.test(perc_diff ~ Genus, data = d_proportion_summary)
print(kruskal_result)

dunn.test(d_proportion_summary$perc_diff, d_proportion_summary$Genus,
          method = "bh", kw = TRUE, label = TRUE)

print(dunn_result)

genus_dividers <- hs_bar %>%
  distinct(x_label, Genus) %>%
  arrange(x_label) %>%
  mutate(x_pos = as.numeric(x_label)) %>%
  group_by(Genus) %>%
  summarise(x_start = min(x_pos), x_end = max(x_pos), x_mid = mean(x_pos), .groups = "drop")

divider_positions <- genus_dividers$x_end[-nrow(genus_dividers)] + 0.5

F3B <- ggplot(hs_bar,
                     aes(x    = x_label,
                         y    = prop,
                         fill = DIV_plot)) +
  
  geom_vline(xintercept = divider_positions,
             color = "grey50", linewidth = 0.5, linetype = "dashed") +
  
  geom_bar(stat = "identity", position = "fill", width = 0.85, color = "white", linewidth = 0.1) +
  annotate("text",
           x        = genus_dividers$x_mid,
           y        = 1.06,
           label    = c("Green Star Polyp", "Pulse Coral", "Leather Coral"),
           size     = 3.8,
           fontface = "bold",
           color    = "grey20") +
  scale_fill_manual(values = div_colors, name = "ITS2 DIV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.10)),
                     breaks = seq(0, 1, by = 0.25)) +
  scale_x_discrete(labels = function(x) gsub("\n", " ", x)) +
  labs(x = NULL, y = "Relative Abundance") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8, angle = 35, hjust = 1, lineheight = 0.85),
    axis.text.y        = element_text(size = 10),
    axis.title.y       = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 7.5),
    legend.key.size    = unit(0.35, "cm"),
    plot.caption       = element_text(size = 9, color = "grey40", hjust = 0)
  ) +
  coord_cartesian(ylim = c(0, 1.1), clip = "off")

print(F3B)

quartz(w=6.5,h=4)
F3A
quartz.save("./Fig2A_final.pdf", type="pdf")
quartz(w=6.5,h=3)
F3B
quartz.save("./Fig2B_final.pdf", type="pdf")

quartz(w=11.7,h=4)
plot_grid(F3A,F3B, nrow=1, labels=c("A","B"), rel_widths = c(0.8,1),label_size=12)
quartz.save("./Fig3_final.pdf", type="pdf")

#####Fig 4A - heat tolerance vs. coral####
heatstress <- read_excel("R_data.xlsx", 1, na = "NA")

score_variable <- function(x, pattern_scores) {
  case_when(
    grepl(pattern_scores[1], x) ~ 1,  # If first pattern found (e.g., "N" for none) → score = 1
    grepl(pattern_scores[2], x) ~ 2,  # If second pattern found (e.g., "M" for medium) → score = 2
    grepl(pattern_scores[3], x) ~ 3,  # If third pattern found (e.g., "Y" for yes) → score = 3
    TRUE ~ NA_real_                   # If no pattern matches → return NA (as numeric)
  )
}

heatstress1 <- heatstress %>%
  separate(time_point, into = c("trash", "time_point"), sep = " ") %>%
  filter(
    time_point == "10:00:00",
    tank != "1",
    date_sampled == as.Date("2025-06-17")
  ) %>%
  dplyr::select(-trash) %>%
  # Score all health metrics
  mutate(
    # Polyp tentacles extension (N=none, M=medium, Y=yes)
    tentacles_score = score_variable(polyps_tentacles, c("N", "M", "Y")),
    # Polyp stalks extension (N=none, M=medium, Y=yes)
    stalks_score = score_variable(polyps_stalks, c("N", "M", "Y")),
    # Average polyp health
    avg_polyps = case_when(
      !is.na(stalks_score) & !is.na(tentacles_score) ~ (stalks_score + tentacles_score) / 2,
      is.na(stalks_score) & !is.na(tentacles_score) ~ tentacles_score,
      !is.na(stalks_score) & is.na(tentacles_score) ~ stalks_score,
      TRUE ~ NA_real_
    ),
    # Mortality (D=dead, PD=partially dead, ND=not dead)
    mortality_score = case_when(
      grepl("^ND", mortality) ~ 3,  # Check ND before D to avoid misclassification
      grepl("PD", mortality) ~ 2,
      grepl("D", mortality) ~ 1,
      TRUE ~ NA_real_
    ),
    # Bleaching (B=bleached, PB=partially bleached, NB=not bleached)
    bleaching_score = case_when(
      grepl("^NB", bleaching) ~ 3,  # Check NB before B to avoid misclassification
      grepl("PB", bleaching) ~ 2,
      grepl("B", bleaching) ~ 1,
      TRUE ~ NA_real_
    ),
    # Calculate overall tolerance score (sum of all health metrics)
    tolerance_score = rowSums(pick(avg_polyps, mortality_score, bleaching_score), na.rm = TRUE)
  ) %>%
  dplyr::select(coral_replicate, Genus, tank, Colony, original_location, tolerance_score)

#check for any tank differences in heat tolerance!
tank_check <- heatstress1 %>%
  filter(tank != 1) %>%
  select(Colony, tank, tolerance_score) %>%
  group_by(Colony, tank) %>%
  summarise(tolerance_score = mean(tolerance_score, na.rm = TRUE), .groups = "drop") %>%
  mutate(tank = as.factor(tank)) %>% 
  mutate(tank = dplyr::recode(as.character(tank),
                       "2" = "Tank2",
                       "3" = "Tank3")) %>%
  pivot_wider(names_from = tank, values_from = tolerance_score)

wilcox.test(tank_check$Tank2, tank_check$Tank3, paired = TRUE)

# Summary stats per colony
colony_summary <- heatstress1 %>%
  group_by(Colony, Genus, original_location) %>%
  summarise(
    mean_score = mean(tolerance_score, na.rm = TRUE),
    se_score   = sd(tolerance_score, na.rm = TRUE) / sqrt(n()),
    n          = n(),
    .groups    = "drop"
  )

colony_levels <- colony_summary %>%
  mutate(Genus = factor(Genus, levels = genus_order)) %>%
  arrange(Genus, desc(mean_score)) %>%
  pull(Colony)

heatstress1 <- heatstress1 %>%
  mutate(Colony = factor(Colony, levels = colony_levels))

colony_summary <- colony_summary %>%
  mutate(Colony = factor(Colony, levels = colony_levels))

# Store colors for x-axis labels
store_colors <- c(
  "Tropical Reef Fish Store" = "#EE6677",   # red
  "Jan's Tropical Fish Store" = "#CCBB44"   # yellow
)

# Build caption with color-coded store key 
store_caption <- paste0(
  "\u25a0 Tropical Reef Fish Store  \u25a0 Jan's Tropical Fish Store"
)

# X-axis label colors by store
label_colors <- colony_levels %>%
  tibble(Colony = .) %>%
  left_join(colony_summary %>% select(Colony, original_location) %>%
              mutate(Colony = as.character(Colony)), by = "Colony") %>%
  mutate(color = store_colors[original_location]) %>%
  pull(color)

# Genus divider positions
genus_breaks <- tibble(Colony = colony_levels) %>%
  left_join(colony_summary %>% select(Colony, Genus) %>%
              mutate(Colony = as.character(Colony)), by = "Colony") %>%
  mutate(x_pos = row_number(),
         Genus = factor(Genus, levels = genus_order)) %>%
  group_by(Genus) %>%
  summarise(x_end = max(x_pos), x_mid = mean(x_pos), .groups = "drop") %>%
  arrange(x_end)

divider_positions_1 <- genus_breaks$x_end[-nrow(genus_breaks)] + 0.5

colored_caption <- paste0(
  "<span style='color:#EE6677'>■</span> Tropical Reef Fish Store &nbsp;&nbsp;", 
  "<span style='color:#CCBB44'>■</span> Jan's Tropical Fish Store"
)

# Plot
F4A <- ggplot() +
  
  geom_vline(xintercept = as.numeric(divider_positions_1),
             color = "grey50", linewidth = 0.5, linetype = "dashed") +
  annotate("text",
           x        = c("GSP1", "PC3", "LC3"),  # middle colony of each genus group
           y        = 10,
           label    = c("Green Star Polyp", "Pulse Coral", "Leather Coral"),
           size     = 4,
           fontface = "bold",
           color    = "grey20", 
           hjust    = c(0.5, 0.5, 0.8) )+
  geom_jitter(data = heatstress1,
              aes(x = Colony, y = tolerance_score, color = Genus),
              width = 0.12, size = 2, alpha = 0.6) +
  geom_errorbar(data = colony_summary,
                aes(x = Colony,
                    ymin = mean_score - se_score,
                    ymax = mean_score + se_score),
                width = 0.2, linewidth = 0.6, color = "grey20") +
  geom_point(data = colony_summary,
             aes(x = Colony, y = mean_score, fill = Genus),
             shape = 21, size = 4, color = "grey20", stroke = 0.5) +
  scale_color_manual(values = genus_colors, name = "Genus") +
  scale_fill_manual(values  = genus_colors, name = "Genus") +
  scale_y_continuous(
    limits = c(0, 10.5),    # extra headroom for genus labels
    breaks = 1:9,
    labels = 1:9,
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  coord_cartesian(clip = "off") +   # prevents labels being clipped at plot edge 
  labs(
    x       = "Colony",
    y       = "Health Score",
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.x     = element_text(size = 13),
    axis.title.y     = element_text(size = 13),
    axis.text.x      = element_text(size = 10, angle = 35, hjust = 1,
                                    color = label_colors),
    axis.text.y      = element_text(size = 10),
    axis.ticks.x     = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.key.size  = unit(0.4, "cm")
  )

# Try ggtext for colored caption (install if needed: install.packages("ggtext"))
if (requireNamespace("ggtext", quietly = TRUE)) {
  library(ggtext)
  colored_caption <- paste0(
    "<span style='color:#EE6677'>&#9632;</span> Tropical Reef Fish Store",
    "<span style='color:#CCBB44'>&#9632;</span> Jan's Tropical Fish Store"
  )
  F4A <- F4A +
    labs(caption = colored_caption) +
    theme(plot.caption = ggtext::element_markdown(size = 12, hjust = 0, color = "grey30"))
} else {
  message("Install 'ggtext' for colored caption squares: install.packages('ggtext')")
  # plain text fallback already set above
}
print(F4A)

kw_coralheat <- kruskal.test(tolerance_score ~ Genus, data = heatstress1)
dunn.test(heatstress1$tolerance_score, heatstress1$Genus,
          method = "bh", kw = TRUE, label = TRUE)

print(kw_store)

# Summary stats for figure caption — saved to CSV

# Per-colony
colony_stats <- colony_summary %>%
  arrange(Genus, desc(mean_score)) %>%
  select(Colony, Genus, original_location, n, mean_score, se_score) %>%
  mutate(level = "Colony")

# Per-genus
genus_stats <- heatstress1 %>%
  group_by(Genus) %>%
  summarise(
    n          = n(),
    mean_score = mean(tolerance_score, na.rm = TRUE),
    se_score   = sd(tolerance_score,   na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  ) %>%
  mutate(Colony = Genus, original_location = NA, level = "Genus")

# Overall
overall_stats <- heatstress1 %>%
  summarise(
    n          = n(),
    mean_score = mean(tolerance_score, na.rm = TRUE),
    se_score   = sd(tolerance_score,   na.rm = TRUE) / sqrt(n())
  ) %>%
  mutate(Colony = "All", Genus = "All", original_location = NA, level = "Overall")

# Combine and round
all_stats <- bind_rows(colony_stats, genus_stats, overall_stats) %>%
  mutate(across(c(mean_score, se_score), ~ round(.x, 3))) %>%
  select(level, Colony, Genus, original_location, n, mean_score, se_score)

#####Fig 4B - thermal stress vs. symbiont dominance#####
# Classify each colony as Cladocopium or Durusdinium dominant
species_map <- c(
  "A" = "Symbiodinium",
  "B" = "Breviolum",
  "C" = "Cladocopium",
  "D" = "Durusdinium",
  "F" = "Fugacium",
  "G" = "Gerakladium",
  "I" = "Clade I"
)

hs_final <- ITS2heatstress_long %>%
  group_by(Sample) %>%
  mutate(prop = reads / sum(reads, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Colony, Genus, Experiment_Type, DIV) %>%
  summarise(avg_prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    clade_letter = str_extract(DIV, "^[A-Z]"),
    full_species = species_map[clade_letter]
  ) %>%
  group_by(Colony, Experiment_Type) %>%
  mutate(
    clado_prop = sum(avg_prop[clade_letter == "C"], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(Genus, Colony, Experiment_Type, clado_prop) %>%
  mutate(dominant = if_else(clado_prop >= 0.5, "Cladocopium", "Durusdinium")) %>%
  distinct()

#Fix naming convention issues across used datasets
hs_final <- hs_final %>%
  left_join(.,heatstress1,by="Genus", relationship = "many-to-many") %>%
  mutate(Colony = Colony.x) %>%
  filter(Experiment_Type!="Control") %>%
  select(Colony, Genus, Experiment_Type, clado_prop, dominant, tolerance_score) %>%
  distinct()

dominant_colors <- c("Cladocopium" = "#4477AA", "Durusdinium" = "#EE6677")

F4B <- ggplot() +
  geom_boxplot(data = hs_final,
               aes(x   = dominant,
                   y   = tolerance_score,
                   fill = dominant),
               width        = 0.4,
               alpha        = 0.25,
               linewidth    = 0.5,
               outlier.shape = NA,
               color        = "grey20") +
  geom_jitter(data  = hs_final,
              aes(x     = dominant,
                  y     = tolerance_score,
                  color = dominant,
                  shape = Genus),
              width = 0.12, size = 3, alpha = 0.8) +
  scale_fill_manual(values  = dominant_colors, name = "Dominant Symbiont") +
  scale_color_manual(values = dominant_colors, name = "Dominant Symbiont") +
  scale_shape_manual(values = genus_shapes,    name = "Genus") +
  scale_y_continuous(breaks = 1:9,
                     limits = c(0, 10),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Dominant Symbiont",
       y = "Health Score") +
  guides(
    fill  = "none",   # redundant with color
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2,
                         override.aes = list(color = "grey30", size = 3))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text        = element_text(size = 11),
    axis.title       = element_text(size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.key.size  = unit(0.4, "cm")
  )

print(F4B)

kruskal.test(tolerance_score ~ dominant, data = hs_final)

quartz(w=6.5,h=4)
F4A
#quartz.save("./Fig4A_final.pdf", type="pdf")
quartz(w=6,h=3)
F4B
#quartz.save("./Fig4B_final.pdf", type="pdf")

quartz(w=11.5,h=4)
plot_grid(F4A,F4B, nrow=1, labels=c("A","B"), rel_widths = c(1.4,1),label_size=12)
quartz.save("./figs/Fig4_final2.pdf", type="pdf")

#####Fig S1 heat tolerance time series#####
# Load data
heatstress <- read_excel("R_data.xlsx", 1, na="NA")
temp_data  <- read_excel("R_data.xlsx", 2)

# Scoring function (matches your original code)
score_variable <- function(x, pattern_scores) {
  case_when(
    grepl(pattern_scores[1], trimws(x)) ~ 1,
    grepl(pattern_scores[2], trimws(x)) ~ 2,
    grepl(pattern_scores[3], trimws(x)) ~ 3,
    TRUE ~ NA_real_
  )
}

# Calculate tolerance scores across full dataset
heatstress_scored <- heatstress %>%
  mutate(
    polyps_tentacles = trimws(polyps_tentacles),
    polyps_stalks    = trimws(polyps_stalks),
    mortality        = trimws(mortality),
    bleaching        = trimws(bleaching),
    tentacles_score = score_variable(polyps_tentacles, c("N", "M", "Y")),
    stalks_score    = score_variable(polyps_stalks,    c("N", "M", "Y")),
    avg_polyps = case_when(
      !is.na(stalks_score) & !is.na(tentacles_score) ~ (stalks_score + tentacles_score) / 2,
      is.na(stalks_score)  & !is.na(tentacles_score) ~ tentacles_score,
      !is.na(stalks_score) & is.na(tentacles_score)  ~ stalks_score,
      TRUE ~ NA_real_
    ),
    mortality_score = case_when(
      grepl("^ND", mortality) ~ 3,
      grepl("PD",  mortality) ~ 2,
      grepl("^D",  mortality) ~ 1,
      TRUE ~ NA_real_
    ),
    bleaching_score = case_when(
      grepl("^NB", bleaching) ~ 3,
      grepl("PB",  bleaching) ~ 2,
      grepl("^B",  bleaching) ~ 1,
      TRUE ~ NA_real_
    ),
    tolerance_score = rowSums(pick(avg_polyps, mortality_score, bleaching_score),
                              na.rm = TRUE),
    treatment = if_else(tank == 1, "Control", "Heated"),
    datetime = as.POSIXct(paste(date_sampled, time_point),
                          format = "%Y-%m-%d %H:%M:%S")
  )

# Daily mean score per colony x treatment x date
# Use 10:00 timepoint only for daily snapshot (consistent with original analysis)
daily_scores <- heatstress_scored %>%
  filter(format(time_point, "%H:%M:%S") == "10:00:00") %>%
  group_by(Colony, Genus, original_location, treatment, date_sampled) %>%
  summarise(
    mean_score = mean(tolerance_score, na.rm = TRUE),
    se_score   = sd(tolerance_score,   na.rm = TRUE) / sqrt(n()),
    n          = n(),
    .groups    = "drop"
  )

# Daily mean score per genus x treatment x date
daily_genus <- heatstress_scored %>%
  filter(format(time_point, "%H:%M:%S") == "10:00:00") %>%
  group_by(Genus, Colony, treatment, date_sampled) %>%
  summarise(
    mean_score = mean(tolerance_score, na.rm = TRUE),
    se_score   = sd(tolerance_score,   na.rm = TRUE) / sqrt(n()),
    n          = n(),
    .groups    = "drop"
  )

# Temperature overlay (heated treatment only, 10:00 timepoint) 
temp_heated <- temp_data %>%
  filter(treatment == "heated") %>%
  mutate(
    date_sampled = as.Date(date_time),
    time         = format(date_time, "%H:%M:%S")
  ) %>%
  filter(time == "10:00:00") %>%
  select(date_sampled, temperature)

treatment_linetypes <- c("Control" = "dashed", "Heated" = "solid")

daily_genus_filtered <- daily_genus %>%
  filter(!(Colony == "PC2" & date_sampled < as.Date("2025-06-12"))) %>%
  filter(!(Colony == "PC3" & date_sampled < as.Date("2025-06-12"))) %>%
  select(Genus, date_sampled, treatment, mean_score, se_score)
  
S1A <- ggplot(daily_genus_filtered, 
                        aes(x = date_sampled, y = mean_score, 
                            color = Genus, fill = Genus)) +
  geom_vline(xintercept = as.Date("2025-06-12"), 
             color = "black", linewidth = 0.5, linetype = "dashed") +
  stat_summary(aes(group = interaction(Genus, treatment)),
               fun.data = mean_se, 
               geom     = "ribbon", 
               alpha    = 0.15, 
               color    = NA) + # NA prevents a border around the ribbon
  stat_summary(aes(group = interaction(Genus, treatment), linetype = treatment),
               fun       = mean, 
               geom      = "line", 
               linewidth = 0.9) +
  stat_summary(aes(shape = treatment),
               fun  = mean, 
               geom = "point", 
               size = 2.5) +
  scale_color_manual(values = genus_colors, name = "Genus") +
  scale_fill_manual(values  = genus_colors, name = "Genus") +
  scale_linetype_manual(values = treatment_linetypes, name = "Treatment") +
  scale_shape_manual(values = c("Control" = 1, "Heated" = 16), name = "Treatment") +
  scale_x_date(date_labels = "%b %d", date_breaks = "1 day") +
  labs(
    x       = "Date",
    y       = "Mean Health Score",
    caption = "Mean ± SE per Genus. Solid = heated, dashed = control.\nDotted grey line = heated tank temperature (right axis). Score range 1–9."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x             = element_text(angle = 35, hjust = 1, size = 10),
    axis.text.y             = element_text(size = 10),
    panel.grid.minor        = element_blank(),
    legend.position         = "right"
  )

print(S1A) #This is the one for paper! 

last_day_comp <- daily_genus %>%
  mutate(date_dt = as.Date(date_sampled, format = "%m/%d/%Y")) %>%
  filter(date_dt == as.Date("2025-06-17")) %>%
  select(Colony, Genus, treatment, mean_score) %>%
  pivot_wider(names_from = treatment, values_from = mean_score) %>%
  mutate(diff = Heated - Control)

overall_test <- wilcox.test(last_day_comp$Heated, 
                            last_day_comp$Control, 
                            paired = TRUE, 
                            alternative = "less")

last_day_genus_summary <- last_day_comp %>%
  group_by(Genus) %>%
  summarise(
    avg_heated_score = mean(Heated, na.rm = TRUE),
    se_heated        = sd(Heated, na.rm = TRUE) / sqrt(n()),
    n_colonies       = n(),
    .groups = "drop"
  )

print(last_day_genus_summary)

quartz(w=6,h=4)
S1A
quartz.save("./figs/FigS1_final.pdf", type = "pdf")