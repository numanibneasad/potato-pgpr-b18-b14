############################################################
# FINAL PUBLICATION-READY COMBINED FIGURE
# A. Plant establishment response:
#    A1. Emergence (%)
#    A2. Shoot no./tuber
# B. Shoot height dynamics
# C. Leaf number dynamics
# D. Chlorophyll index (SPAD) dynamics
############################################################

library(tidyverse)
library(multcompView)
library(patchwork)

# =========================
# Common settings
# =========================

treat_order <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only",
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

no_path_order <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only"
)

with_path_order <- c(
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

treat_cols <- c(
  "Healthy control"       = "black",
  "PGPR A only"           = "#1F78B4",
  "PGPR B only"           = "#33A02C",
  "PGPR A + B only"       = "#A6CEE3",
  "Pathogen only"         = "#E31A1C",
  "PGPR A + Pathogen"     = "#FDBF6F",
  "PGPR B + Pathogen"     = "#B2DF8A",
  "PGPR A + B + Pathogen" = "#6A3D9A"
)

theme_pub <- theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 17, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    plot.margin = margin(8, 8, 8, 8)
  )

parse_treatment <- function(df) {
  df %>%
    mutate(
      Treatment = factor(as.character(Treatment), levels = treat_order),
      PathogenGroup = ifelse(
        Treatment %in% with_path_order,
        "With Pathogen",
        "No Pathogen"
      ),
      PathogenGroup = factor(PathogenGroup, levels = c("No Pathogen", "With Pathogen")),
      Block = factor(Block)
    )
}

get_tukey_letters <- function(dat, response_col) {
  d <- dat %>% filter(!is.na(.data[[response_col]]))
  
  if (n_distinct(d$Treatment) < 2) {
    return(tibble(Treatment = unique(d$Treatment), Letters = ""))
  }
  
  m <- aov(reformulate(c("Treatment", "Block"), response = response_col), data = d)
  tk <- TukeyHSD(m, "Treatment")
  pvals <- tk$Treatment[, "p adj"]
  names(pvals) <- rownames(tk$Treatment)
  let <- multcompLetters(pvals)$Letters
  
  tibble(
    Treatment = names(let),
    Letters = tolower(gsub("[^A-Za-z]", "", unname(let)))
  )
}

# ============================================================
# A. PLANT ESTABLISHMENT RESPONSE
# ============================================================

shoot_df <- read.delim(
  "shoot.count.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

shoot_df <- shoot_df %>%
  rename(
    Shoot_20Nov = `20-Nov`,
    Shoot_27Nov = `27-Nov`
  ) %>%
  parse_treatment() %>%
  mutate(
    Shoot_20Nov = as.numeric(Shoot_20Nov),
    Shoot_27Nov = as.numeric(Shoot_27Nov),
    Emerged = ifelse(Shoot_20Nov > 0, 1, 0)
  )

# Emergence (%)
emerg_sum <- shoot_df %>%
  group_by(PathogenGroup, Treatment) %>%
  summarise(
    p = mean(Emerged, na.rm = TRUE),
    n = n(),
    Mean = p * 100,
    SE = sqrt(p * (1 - p) / n) * 100,
    .groups = "drop"
  )

letters_emerg <- shoot_df %>%
  group_by(PathogenGroup) %>%
  group_modify(~ get_tukey_letters(.x, "Emerged")) %>%
  ungroup()

emerg_sum <- emerg_sum %>%
  left_join(letters_emerg, by = c("PathogenGroup", "Treatment")) %>%
  mutate(Letters = replace_na(Letters, "")) %>%
  group_by(PathogenGroup) %>%
  mutate(label_y = Mean + SE + 5) %>%
  ungroup()

# Raw shoot number per tuber
shoot_sum <- shoot_df %>%
  group_by(PathogenGroup, Treatment) %>%
  summarise(
    Mean = mean(Shoot_27Nov, na.rm = TRUE),
    SD = sd(Shoot_27Nov, na.rm = TRUE),
    n = sum(!is.na(Shoot_27Nov)),
    SE = SD / sqrt(n),
    .groups = "drop"
  )

letters_shoot <- shoot_df %>%
  group_by(PathogenGroup) %>%
  group_modify(~ get_tukey_letters(.x, "Shoot_27Nov")) %>%
  ungroup()

shoot_sum <- shoot_sum %>%
  left_join(letters_shoot, by = c("PathogenGroup", "Treatment")) %>%
  mutate(Letters = replace_na(Letters, "")) %>%
  group_by(PathogenGroup) %>%
  mutate(label_y = Mean + SE + 0.15 * max(Mean + SE, na.rm = TRUE)) %>%
  ungroup()

plot_establishment <- function(df, ylab, title_text, ylim = NULL) {
  ggplot(df, aes(x = Treatment, y = Mean, fill = Treatment)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.35) +
    geom_errorbar(
      aes(ymin = Mean - SE, ymax = Mean + SE),
      width = 0.18,
      linewidth = 0.65
    ) +
    geom_text(
      aes(label = Letters, y = label_y),
      size = 5,
      fontface = "bold",
      color = "black"
    ) +
    facet_grid(. ~ PathogenGroup, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = treat_cols, breaks = treat_order, drop = FALSE) +
    scale_x_discrete(drop = TRUE) +
    labs(
      title = title_text,
      x = NULL,
      y = ylab,
      fill = "Treatment"
    ) +
    theme_pub +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1, face = "bold")
    ) +
    coord_cartesian(ylim = ylim, clip = "off")
}

pA1 <- plot_establishment(
  emerg_sum,
  ylab = "Emergence (%)",
  title_text = "A. Plant establishment response",
  ylim = c(0, 115)
)

pA2 <- plot_establishment(
  shoot_sum,
  ylab = "Shoot no./tuber",
  title_text = NULL
)

pA <- pA1 / pA2 + plot_layout(heights = c(1, 1))

# ============================================================
# B. SHOOT HEIGHT DYNAMICS
# ============================================================

height_df <- read.delim(
  "Shoot_length.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

height_dates <- c("20-Nov", "27-Nov", "03-Dec", "10-Dec", "17-Dec", "23-Dec")

height_long <- height_df %>%
  pivot_longer(
    cols = all_of(height_dates),
    names_to = "Date",
    values_to = "Height"
  ) %>%
  parse_treatment() %>%
  mutate(
    Date = factor(Date, levels = height_dates, ordered = TRUE),
    Height = as.numeric(Height)
  )

height_sum <- height_long %>%
  group_by(PathogenGroup, Treatment, Date) %>%
  summarise(
    Mean = mean(Height, na.rm = TRUE),
    SD = sd(Height, na.rm = TRUE),
    n = sum(!is.na(Height)),
    SE = SD / sqrt(n),
    .groups = "drop"
  )

pB <- ggplot(height_sum, aes(x = Date, y = Mean, group = Treatment, color = Treatment)) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.10, linewidth = 0.65) +
  facet_wrap(~ PathogenGroup, scales = "free_y") +
  scale_color_manual(values = treat_cols, breaks = treat_order, drop = FALSE) +
  labs(
    title = "B. Shoot height dynamics",
    x = "Date",
    y = "Shoot height (cm)",
    color = "Treatment"
  ) +
  theme_pub +
  theme(legend.position = "none")

# ============================================================
# C. LEAF NUMBER DYNAMICS
# ============================================================

leaf_df <- read.delim(
  "Leaves_number.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

leaf_long <- leaf_df %>%
  pivot_longer(
    cols = all_of(height_dates),
    names_to = "Date",
    values_to = "Leaves"
  ) %>%
  parse_treatment() %>%
  mutate(
    Date = factor(Date, levels = height_dates, ordered = TRUE),
    Leaves = as.numeric(Leaves)
  )

leaf_sum <- leaf_long %>%
  group_by(PathogenGroup, Treatment, Date) %>%
  summarise(
    Mean = mean(Leaves, na.rm = TRUE),
    SD = sd(Leaves, na.rm = TRUE),
    n = sum(!is.na(Leaves)),
    SE = SD / sqrt(n),
    .groups = "drop"
  )

pC <- ggplot(leaf_sum, aes(x = Date, y = Mean, group = Treatment, color = Treatment)) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.10, linewidth = 0.65) +
  facet_wrap(~ PathogenGroup, scales = "free_y") +
  scale_color_manual(values = treat_cols, breaks = treat_order, drop = FALSE) +
  labs(
    title = "C. Leaf number dynamics",
    x = "Date",
    y = "Leaf number",
    color = "Treatment"
  ) +
  theme_pub +
  theme(legend.position = "none")

# ============================================================
# D. CHLOROPHYLL INDEX (SPAD) DYNAMICS
# ============================================================

chl_df <- read.delim(
  "chlorophyll.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

chl_dates <- c("03-Dec", "10-Dec", "23-Dec")

chl_long <- chl_df %>%
  pivot_longer(
    cols = all_of(chl_dates),
    names_to = "Date",
    values_to = "Chlorophyll"
  ) %>%
  parse_treatment() %>%
  mutate(
    Date = factor(Date, levels = chl_dates, ordered = TRUE),
    Chlorophyll = as.numeric(Chlorophyll)
  )

chl_sum <- chl_long %>%
  group_by(PathogenGroup, Treatment, Date) %>%
  summarise(
    Mean = mean(Chlorophyll, na.rm = TRUE),
    SD = sd(Chlorophyll, na.rm = TRUE),
    n = sum(!is.na(Chlorophyll)),
    SE = SD / sqrt(n),
    .groups = "drop"
  )

pD <- ggplot(chl_sum, aes(x = Date, y = Mean, group = Treatment, color = Treatment)) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.10, linewidth = 0.65) +
  facet_wrap(~ PathogenGroup, scales = "free_y") +
  scale_color_manual(values = treat_cols, breaks = treat_order, drop = FALSE) +
  labs(
    title = "D. Chlorophyll index (SPAD) dynamics",
    x = "Date",
    y = "Chlorophyll index (SPAD)",
    color = "Treatment"
  ) +
  theme_pub

# ============================================================
# FINAL PANEL
# ============================================================

final_panel <- 
  (pA | pB) /
  (pC | pD) +
  plot_layout(
    guides = "collect",
    widths = c(1.25, 1),
    heights = c(1.15, 1)
  ) &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15)
  )

print(final_panel)

ggsave(
  "Figure_Establishment_Growth_Physiology_FINAL_corrected.tiff",
  final_panel,
  width = 20,
  height = 16,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  "Figure_Establishment_Growth_Physiology_FINAL_corrected.pdf",
  final_panel,
  width = 20,
  height = 16
)

ggsave(
  "Figure_Establishment_Growth_Physiology_FINAL_corrected.png",
  final_panel,
  width = 20,
  height = 16,
  dpi = 600,
  bg = "white"
)

##########################################################################

############################################################
# PUBLICATION-READY TWO-PANEL FIGURE
# A. Flowering incidence (% pots flowered)
# B. Tuber number per plant
############################################################

library(tidyverse)
library(multcompView)
library(patchwork)

# =========================
# Common settings
# =========================

treat_levels <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only",
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

no_path_levels <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only"
)

with_path_levels <- c(
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

pal <- c(
  "Healthy control"         = "#000000",
  "PGPR A only"             = "#1F78B4",
  "PGPR B only"             = "#33A02C",
  "PGPR A + B only"         = "#A6CEE3",
  "Pathogen only"           = "#E31A1C",
  "PGPR A + Pathogen"       = "#FDBF6F",
  "PGPR B + Pathogen"       = "#B2DF8A",
  "PGPR A + B + Pathogen"   = "#6A3D9A"
)

theme_pub <- theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 17, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    plot.margin = margin(8, 8, 8, 8)
  )

# =========================
# Helper function
# =========================

get_tukey_letters <- function(dat, response_col) {
  d <- dat %>% filter(!is.na(.data[[response_col]]))
  
  if (n_distinct(d$Treatment) < 2) {
    return(tibble(Treatment = unique(d$Treatment), Letters = ""))
  }
  
  m <- aov(reformulate(c("Treatment", "Block"), response = response_col), data = d)
  tk <- TukeyHSD(m, "Treatment")
  
  pvals <- tk$Treatment[, "p adj"]
  names(pvals) <- rownames(tk$Treatment)
  
  let <- multcompLetters(pvals)$Letters
  
  tibble(
    Treatment = names(let),
    Letters = tolower(gsub("[^A-Za-z]", "", unname(let)))
  )
}

# ============================================================
# A. FLOWERING INCIDENCE
# ============================================================

inflo <- read.delim(
  "inflorescence.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

date_cols <- c("17-Dec", "23-Dec")

inflo_long <- inflo %>%
  pivot_longer(
    cols = all_of(date_cols),
    names_to = "Date",
    values_to = "Inflorescences"
  ) %>%
  mutate(
    Block = factor(Block),
    Treatment = factor(as.character(Treatment), levels = treat_levels),
    Date = factor(Date, levels = date_cols),
    Inflorescences = as.numeric(Inflorescences),
    Flowered = ifelse(Inflorescences > 0, 1, 0),
    PathogenGroup = ifelse(
      Treatment %in% with_path_levels,
      "With Pathogen",
      "No Pathogen"
    ),
    PathogenGroup = factor(PathogenGroup, levels = c("No Pathogen", "With Pathogen"))
  )

letters_inflo <- inflo_long %>%
  filter(!is.na(Flowered)) %>%
  group_by(Date, PathogenGroup) %>%
  group_modify(~ get_tukey_letters(.x, "Flowered")) %>%
  ungroup()

sum_inflo <- inflo_long %>%
  filter(!is.na(Flowered)) %>%
  group_by(Date, PathogenGroup, Treatment) %>%
  summarise(
    n = n(),
    p = mean(Flowered, na.rm = TRUE),
    Mean = p * 100,
    SE = sqrt(p * (1 - p) / n) * 100,
    .groups = "drop"
  ) %>%
  left_join(letters_inflo, by = c("Date", "PathogenGroup", "Treatment")) %>%
  mutate(
    Letters = replace_na(Letters, ""),
    x_plot = factor(as.character(Treatment), levels = treat_levels)
  ) %>%
  group_by(Date, PathogenGroup) %>%
  mutate(label_y = pmin(Mean + SE + 8, 108)) %>%
  ungroup()

p_inflo <- ggplot(sum_inflo, aes(x = x_plot, y = Mean, fill = Treatment)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.35) +
  geom_errorbar(
    aes(ymin = pmax(0, Mean - SE), ymax = pmin(100, Mean + SE)),
    width = 0.18,
    linewidth = 0.65
  ) +
  geom_text(
    aes(label = Letters, y = label_y),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  facet_grid(Date ~ PathogenGroup, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, breaks = treat_levels, drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 20)) +
  labs(
    title = "A. Flowering incidence",
    x = NULL,
    y = "Pots flowered (%)",
    fill = "Treatment"
  ) +
  theme_pub +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")

# ============================================================
# B. TUBER NUMBER PER PLANT
# ============================================================

tuber <- read.delim(
  "Tuber.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

tuber <- tuber %>%
  mutate(
    Block = factor(Block),
    Treatment = factor(as.character(Treatment), levels = treat_levels),
    Tubers = as.numeric(Tubers),
    PathogenGroup = ifelse(
      Treatment %in% with_path_levels,
      "With Pathogen",
      "No Pathogen"
    ),
    PathogenGroup = factor(PathogenGroup, levels = c("No Pathogen", "With Pathogen"))
  )

letters_tuber <- tuber %>%
  group_by(PathogenGroup) %>%
  group_modify(~ get_tukey_letters(.x, "Tubers")) %>%
  ungroup()

sum_tuber <- tuber %>%
  group_by(PathogenGroup, Treatment) %>%
  summarise(
    Mean = mean(Tubers, na.rm = TRUE),
    SD = sd(Tubers, na.rm = TRUE),
    n = sum(!is.na(Tubers)),
    SE = SD / sqrt(n),
    .groups = "drop"
  ) %>%
  left_join(letters_tuber, by = c("PathogenGroup", "Treatment")) %>%
  mutate(
    Letters = replace_na(Letters, ""),
    x_plot = factor(as.character(Treatment), levels = treat_levels)
  ) %>%
  group_by(PathogenGroup) %>%
  mutate(label_y = Mean + SE + 0.12 * max(Mean + SE, na.rm = TRUE)) %>%
  ungroup()

p_tuber <- ggplot(sum_tuber, aes(x = x_plot, y = Mean, fill = Treatment)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.35) +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    width = 0.18,
    linewidth = 0.65
  ) +
  geom_text(
    aes(label = Letters, y = label_y),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  facet_grid(. ~ PathogenGroup, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, breaks = treat_levels, drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  labs(
    title = "B. Tuber number",
    x = NULL,
    y = "Tubers per plant",
    fill = "Treatment"
  ) +
  theme_pub +
  coord_cartesian(clip = "off")

# ============================================================
# FINAL TWO-PANEL FIGURE
# ============================================================

final_panel <- p_inflo / p_tuber +
  plot_layout(
    guides = "collect",
    heights = c(1.25, 1)
  ) &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15)
  )

print(final_panel)

ggsave(
  "Figure_Flowering_Tuber_Publication.tiff",
  final_panel,
  width = 16,
  height = 14,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  "Figure_Flowering_Tuber_Publication.pdf",
  final_panel,
  width = 16,
  height = 14
)

ggsave(
  "Figure_Flowering_Tuber_Publication.png",
  final_panel,
  width = 16,
  height = 14,
  dpi = 600,
  bg = "white"
)
##########################################################################
############################################################
# PUBLICATION-READY FIGURE
# Fresh biomass, dry biomass, and WinRHIZO root architecture
# RCBD ANOVA + Tukey HSD letters
############################################################

library(tidyverse)
library(multcompView)
library(patchwork)

# ============================================================
# 1. Common settings
# ============================================================

fresh_file <- "fresh.biomass.txt"
dry_file   <- "dry.biomass.txt"
root_file  <- "winrhizodata.txt"

treat_levels <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only",
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

no_path_levels <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only"
)

with_path_levels <- c(
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

pal <- c(
  "Healthy control"         = "#000000",
  "PGPR A only"             = "#1F78B4",
  "PGPR B only"             = "#33A02C",
  "PGPR A + B only"         = "#A6CEE3",
  "Pathogen only"           = "#E31A1C",
  "PGPR A + Pathogen"       = "#FDBF6F",
  "PGPR B + Pathogen"       = "#B2DF8A",
  "PGPR A + B + Pathogen"   = "#6A3D9A"
)

theme_pub <- theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 17, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 2. Helper functions
# ============================================================

prepare_data <- function(df) {
  df %>%
    mutate(
      ID = as.character(ID),
      Block = stringr::str_extract(ID, "R\\s*\\d+"),
      Block = stringr::str_replace_all(Block, "\\s+", ""),
      Block = factor(Block),
      Treatment = factor(as.character(Treatment), levels = treat_levels),
      PathogenGroup = ifelse(
        Treatment %in% with_path_levels,
        "With Pathogen",
        "No Pathogen"
      ),
      PathogenGroup = factor(PathogenGroup, levels = c("No Pathogen", "With Pathogen"))
    )
}

get_tukey_letters <- function(d) {
  d <- d %>%
    filter(!is.na(Value), !is.na(Treatment), !is.na(Block)) %>%
    mutate(
      Treatment = droplevels(Treatment),
      Block = droplevels(Block)
    )
  
  if (n_distinct(d$Treatment) < 2) {
    return(tibble(Treatment = unique(d$Treatment), Letters = ""))
  }
  
  m <- aov(Value ~ Treatment + Block, data = d)
  tk <- TukeyHSD(m, "Treatment")
  
  pvals <- tk$Treatment[, "p adj"]
  names(pvals) <- rownames(tk$Treatment)
  
  let <- multcompLetters(pvals)$Letters
  
  tibble(
    Treatment = names(let),
    Letters = tolower(gsub("[^A-Za-z]", "", unname(let)))
  )
}

summarise_traits <- function(dat_long) {
  
  letters_df <- dat_long %>%
    group_by(Trait, PathogenGroup) %>%
    group_modify(~ get_tukey_letters(.x)) %>%
    ungroup()
  
  dat_long %>%
    group_by(Trait, TraitLabel, PathogenGroup, Treatment) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      n = sum(!is.na(Value)),
      SE = SD / sqrt(n),
      .groups = "drop"
    ) %>%
    left_join(
      letters_df,
      by = c("Trait", "PathogenGroup", "Treatment")
    ) %>%
    mutate(Letters = replace_na(Letters, "")) %>%
    group_by(TraitLabel, PathogenGroup) %>%
    mutate(
      label_y = Mean + SE + 0.12 * max(Mean + SE, na.rm = TRUE)
    ) %>%
    ungroup()
}

plot_one_trait <- function(sum_df, trait_label, panel_title, y_label, show_legend = FALSE) {
  
  d <- sum_df %>%
    filter(TraitLabel == trait_label) %>%
    mutate(
      x_plot = factor(as.character(Treatment), levels = treat_levels)
    )
  
  ggplot(d, aes(x = x_plot, y = Mean, fill = Treatment)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.35) +
    geom_errorbar(
      aes(ymin = Mean - SE, ymax = Mean + SE),
      width = 0.18,
      linewidth = 0.65
    ) +
    geom_text(
      aes(label = Letters, y = label_y),
      size = 5,
      fontface = "bold",
      color = "black"
    ) +
    facet_grid(. ~ PathogenGroup, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal, breaks = treat_levels, drop = FALSE) +
    scale_x_discrete(drop = TRUE) +
    labs(
      title = panel_title,
      x = NULL,
      y = y_label,
      fill = "Treatment"
    ) +
    theme_pub +
    theme(
      legend.position = ifelse(show_legend, "right", "none")
    ) +
    coord_cartesian(clip = "off")
}

# ============================================================
# 3. Read and prepare fresh biomass data
# ============================================================

fresh_raw <- read.delim(
  fresh_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

fresh <- fresh_raw
names(fresh) <- make.names(names(fresh))

fresh <- prepare_data(fresh)

fresh_traits <- c(
  "total_fresh_weight_gm" = "Total fresh biomass",
  "leaves_fresh_weight_gm" = "Leaf fresh weight",
  "steam_fresh_weight_gm" = "Stem fresh weight",
  "leaf_area_cm2" = "Leaf area"
)

fresh_long <- fresh %>%
  pivot_longer(
    cols = all_of(names(fresh_traits)),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    TraitLabel = recode(Trait, !!!fresh_traits)
  ) %>%
  filter(!is.na(Value))

fresh_sum <- summarise_traits(fresh_long)

# ============================================================
# 4. Read and prepare dry biomass data
# ============================================================

dry_raw <- read.delim(
  dry_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

dry <- dry_raw
names(dry) <- make.names(names(dry))

dry <- prepare_data(dry)

dry_traits <- c(
  "dry_whole_biomass_gm" = "Total dry biomass",
  "leaves_dry_weight_gm" = "Leaf dry weight",
  "steam_dry_weight" = "Stem dry weight",
  "root_dry_weight_gm" = "Root dry weight"
)

dry_long <- dry %>%
  pivot_longer(
    cols = all_of(names(dry_traits)),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    TraitLabel = recode(Trait, !!!dry_traits)
  ) %>%
  filter(!is.na(Value))

dry_sum <- summarise_traits(dry_long)

# ============================================================
# 5. Read and prepare WinRHIZO data
# ============================================================

root_raw <- read.delim(
  root_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

root <- root_raw
names(root) <- make.names(names(root))

root <- prepare_data(root)

root_traits <- c(
  "TotalLength.cm." = "Total root length",
  "SurfArea.cm2." = "Root surface area",
  "AvgDiam.mm." = "Average root diameter",
  "RootVolume.cm3." = "Root volume"
)

root_long <- root %>%
  pivot_longer(
    cols = all_of(names(root_traits)),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    TraitLabel = recode(Trait, !!!root_traits)
  ) %>%
  filter(!is.na(Value))

root_sum <- summarise_traits(root_long)

# ============================================================
# 6. Main publication figure
# Selected non-redundant traits for acceptance/readability
# ============================================================

pA <- plot_one_trait(
  fresh_sum,
  "Total fresh biomass",
  "A. Total fresh biomass",
  "Fresh biomass (g plant\u207B\u00B9)"
)

pB <- plot_one_trait(
  dry_sum,
  "Total dry biomass",
  "B. Total dry biomass",
  "Dry biomass (g plant\u207B\u00B9)"
)

pC <- plot_one_trait(
  dry_sum,
  "Root dry weight",
  "C. Root dry weight",
  "Root dry weight (g plant\u207B\u00B9)"
)

pD <- plot_one_trait(
  fresh_sum,
  "Leaf area",
  "D. Leaf area",
  "Leaf area (cm\u00B2 plant\u207B\u00B9)"
)

pE <- plot_one_trait(
  root_sum,
  "Total root length",
  "E. Total root length",
  "Total root length (cm plant\u207B\u00B9)"
)

pF <- plot_one_trait(
  root_sum,
  "Root surface area",
  "F. Root surface area",
  "Root surface area (cm\u00B2 plant\u207B\u00B9)",
  show_legend = TRUE
)

main_fig <- (pA | pB) /
  (pC | pD) /
  (pE | pF) +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15)
  )

print(main_fig)

ggsave(
  "Figure_Biomass_RootArchitecture_Publication.tiff",
  main_fig,
  width = 20,
  height = 21,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  "Figure_Biomass_RootArchitecture_Publication.pdf",
  main_fig,
  width = 20,
  height = 21
)

ggsave(
  "Figure_Biomass_RootArchitecture_Publication.png",
  main_fig,
  width = 20,
  height = 21,
  dpi = 600,
  bg = "white"
)

# ============================================================
# 7. Optional supplementary figure with all measured traits
# ============================================================

pS1 <- plot_one_trait(fresh_sum, "Leaf fresh weight", "A. Leaf fresh weight", "Leaf fresh weight (g plant\u207B\u00B9)")
pS2 <- plot_one_trait(fresh_sum, "Stem fresh weight", "B. Stem fresh weight", "Stem fresh weight (g plant\u207B\u00B9)")
pS3 <- plot_one_trait(dry_sum, "Leaf dry weight", "C. Leaf dry weight", "Leaf dry weight (g plant\u207B\u00B9)")
pS4 <- plot_one_trait(dry_sum, "Stem dry weight", "D. Stem dry weight", "Stem dry weight (g plant\u207B\u00B9)")
pS5 <- plot_one_trait(root_sum, "Average root diameter", "E. Average root diameter", "Average root diameter (mm)")
pS6 <- plot_one_trait(root_sum, "Root volume", "F. Root volume", "Root volume (cm\u00B3 plant\u207B\u00B9)", show_legend = TRUE)

supp_fig <- (pS1 | pS2) /
  (pS3 | pS4) /
  (pS5 | pS6) +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15)
  )

print(supp_fig)

ggsave(
  "Supplementary_Biomass_RootArchitecture_AllTraits.tiff",
  supp_fig,
  width = 20,
  height = 21,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  "Supplementary_Biomass_RootArchitecture_AllTraits.pdf",
  supp_fig,
  width = 20,
  height = 21
)

# ============================================================
# 8. Save summary tables
# ============================================================

write.csv(fresh_sum, "FreshBiomass_Summary_Tukey.csv", row.names = FALSE)
write.csv(dry_sum, "DryBiomass_Summary_Tukey.csv", row.names = FALSE)
write.csv(root_sum, "WinRHIZO_Summary_Tukey.csv", row.names = FALSE)
##########################################################################
############################################################
# PUBLICATION-READY AUDPC FIGURE
# Verticillium wilt disease progression summarized as AUDPC
# RCBD ANOVA + Tukey HSD letters
#
# Notes:
# - No Pathogen group may have AUDPC = 0 for all treatments.
# - If AUDPC has no variation, Tukey letters are left blank.
# - Tukey letters are calculated only when ANOVA/Tukey is valid.
############################################################

# =========================
# 0) Packages
# =========================

pkgs <- c("tidyverse", "multcompView", "patchwork")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]

if(length(to_install) > 0) {
  install.packages(to_install)
}

library(tidyverse)
library(multcompView)
library(patchwork)

# =========================
# 1) Input file
# =========================

file_path <- "disease scores.txt"

df <- read.delim(
  file_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# =========================
# 2) Treatment order and color palette
# Consistent with previous figures
# =========================

treat_levels <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only",
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

no_path_levels <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only"
)

with_path_levels <- c(
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

pal <- c(
  "Healthy control"         = "#000000",
  "PGPR A only"             = "#1F78B4",
  "PGPR B only"             = "#33A02C",
  "PGPR A + B only"         = "#A6CEE3",
  "Pathogen only"           = "#E31A1C",
  "PGPR A + Pathogen"       = "#FDBF6F",
  "PGPR B + Pathogen"       = "#B2DF8A",
  "PGPR A + B + Pathogen"   = "#6A3D9A"
)

# =========================
# 3) Publication theme
# =========================

theme_pub <- theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 16
    ),
    
    axis.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 14,
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 17
    ),
    
    strip.text = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 15
    ),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    
    plot.margin = margin(10, 10, 10, 10)
  )

# =========================
# 4) Clean and prepare data
# =========================

df <- df %>%
  mutate(
    ID = gsub("\\s+", "", ID),
    Block = factor(Block),
    Treatment = factor(as.character(Treatment), levels = treat_levels),
    PathogenGroup = ifelse(
      Treatment %in% with_path_levels,
      "With Pathogen",
      "No Pathogen"
    ),
    PathogenGroup = factor(
      PathogenGroup,
      levels = c("No Pathogen", "With Pathogen")
    )
  )

# Disease-score date columns
date_cols <- c("03-Dec", "10-Dec", "17-Dec", "23-Dec")

# Convert dates to day offsets for AUDPC calculation
year_used <- 2025

date_vec <- as.Date(
  paste(year_used, date_cols),
  format = "%Y %d-%b"
)

day_offsets <- as.numeric(date_vec - min(date_vec))

# =========================
# 5) Long format
# =========================

dat_long <- df %>%
  pivot_longer(
    cols = all_of(date_cols),
    names_to = "Date",
    values_to = "Score"
  ) %>%
  mutate(
    Date = factor(
      Date,
      levels = date_cols,
      ordered = TRUE
    ),
    
    Day = day_offsets[
      match(as.character(Date), date_cols)
    ],
    
    Score = as.numeric(Score)
  )

# =========================
# 6) AUDPC function
# Trapezoidal method
# =========================

audpc_trapz <- function(day, y) {
  
  keep <- !is.na(day) & !is.na(y)
  day <- day[keep]
  y <- y[keep]
  
  if(length(day) < 2) {
    return(NA_real_)
  }
  
  o <- order(day)
  day <- day[o]
  y <- y[o]
  
  sum(
    diff(day) *
      (head(y, -1) + tail(y, -1)) / 2
  )
}

# =========================
# 7) Calculate AUDPC per experimental unit
# =========================

audpc_df <- dat_long %>%
  group_by(ID, Treatment, Block, PathogenGroup) %>%
  summarise(
    AUDPC = audpc_trapz(Day, Score),
    .groups = "drop"
  )

# Save raw AUDPC per pot
write.csv(
  audpc_df,
  "AUDPC_PerPot.csv",
  row.names = FALSE
)

# =========================
# 8) Safe RCBD ANOVA + Tukey letters
# Handles zero-variance groups safely
# =========================

get_tukey_letters <- function(dat) {
  
  dat <- dat %>%
    filter(!is.na(AUDPC)) %>%
    mutate(
      Treatment = droplevels(Treatment),
      Block = droplevels(Block)
    )
  
  # If fewer than two treatments, no comparison possible
  if(n_distinct(dat$Treatment) < 2) {
    return(
      tibble(
        Treatment = unique(dat$Treatment),
        Letters = ""
      )
    )
  }
  
  # If all AUDPC values are identical, ANOVA/Tukey cannot be calculated
  if(length(unique(dat$AUDPC)) < 2 || var(dat$AUDPC, na.rm = TRUE) == 0) {
    
    message(
      "AUDPC has zero variance for group: ",
      unique(dat$PathogenGroup),
      ". Tukey letters omitted."
    )
    
    return(
      tibble(
        Treatment = levels(dat$Treatment),
        Letters = ""
      )
    )
  }
  
  # RCBD ANOVA
  model <- aov(
    AUDPC ~ Treatment + Block,
    data = dat
  )
  
  cat("\n==============================\n")
  cat("ANOVA for group:", as.character(unique(dat$PathogenGroup)), "\n")
  cat("==============================\n")
  print(anova(model))
  
  aov_tab <- anova(model)
  p_treat <- aov_tab["Treatment", "Pr(>F)"]
  
  # If treatment effect is not significant, all treatments receive same letter
  if(is.na(p_treat) || p_treat >= 0.05) {
    
    message(
      "Treatment effect not significant for group: ",
      unique(dat$PathogenGroup),
      ". All treatments assigned letter 'a'."
    )
    
    return(
      tibble(
        Treatment = levels(dat$Treatment),
        Letters = "a"
      )
    )
  }
  
  # Tukey HSD
  tuk <- TukeyHSD(
    model,
    "Treatment"
  )
  
  cat("\nTukey HSD for group:", as.character(unique(dat$PathogenGroup)), "\n")
  print(tuk)
  
  pvals <- tuk$Treatment[, "p adj"]
  names(pvals) <- rownames(tuk$Treatment)
  
  # Remove NA p-values if any
  pvals <- pvals[!is.na(pvals)]
  
  if(length(pvals) == 0) {
    return(
      tibble(
        Treatment = levels(dat$Treatment),
        Letters = ""
      )
    )
  }
  
  let <- multcompLetters(pvals)$Letters
  
  tibble(
    Treatment = names(let),
    Letters = tolower(
      gsub("[^A-Za-z]", "", unname(let))
    )
  )
}

letters_audpc <- audpc_df %>%
  group_by(PathogenGroup) %>%
  group_modify(~ get_tukey_letters(.x)) %>%
  ungroup()

# =========================
# 9) Summary table for plotting
# =========================

sum_audpc <- audpc_df %>%
  group_by(PathogenGroup, Treatment) %>%
  summarise(
    Mean = mean(AUDPC, na.rm = TRUE),
    SD = sd(AUDPC, na.rm = TRUE),
    n = sum(!is.na(AUDPC)),
    SE = SD / sqrt(n),
    .groups = "drop"
  ) %>%
  left_join(
    letters_audpc,
    by = c("PathogenGroup", "Treatment")
  ) %>%
  mutate(
    Letters = replace_na(Letters, ""),
    x_plot = factor(
      as.character(Treatment),
      levels = treat_levels
    )
  ) %>%
  group_by(PathogenGroup) %>%
  mutate(
    label_y = case_when(
      max(Mean + SE, na.rm = TRUE) == 0 ~ 0.05,
      TRUE ~ Mean + SE + 0.10 * max(Mean + SE, na.rm = TRUE)
    )
  ) %>%
  ungroup()

write.csv(
  sum_audpc,
  "AUDPC_Summary_RCBD_Tukey.csv",
  row.names = FALSE
)

# =========================
# 10) Publication-ready AUDPC plot
# =========================

p_audpc <- ggplot(
  sum_audpc,
  aes(
    x = x_plot,
    y = Mean,
    fill = Treatment
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.35
  ) +
  geom_errorbar(
    aes(
      ymin = Mean - SE,
      ymax = Mean + SE
    ),
    width = 0.20,
    linewidth = 0.70
  ) +
  geom_text(
    aes(
      label = Letters,
      y = label_y
    ),
    size = 5.2,
    fontface = "bold",
    color = "black",
    na.rm = TRUE
  ) +
  facet_grid(
    . ~ PathogenGroup,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = pal,
    breaks = treat_levels,
    drop = FALSE
  ) +
  scale_x_discrete(
    drop = TRUE
  ) +
  labs(
    title = "AUDPC response to PGPR and pathogen treatments",
    x = NULL,
    y = "AUDPC (disease severity × days)",
    fill = "Treatment"
  ) +
  theme_pub +
  coord_cartesian(
    clip = "off"
  )

print(p_audpc)

# =========================
# 11) Save publication-quality outputs
# =========================

ggsave(
  "Figure_AUDPC_Publication.tiff",
  p_audpc,
  width = 12,
  height = 7,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  "Figure_AUDPC_Publication.pdf",
  p_audpc,
  width = 12,
  height = 7
)

ggsave(
  "Figure_AUDPC_Publication.png",
  p_audpc,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white"
)

# =========================
# 12) Suggested caption
# =========================

cat(
  "\nSuggested caption:\n",
  "Figure X. Effects of PGPR treatments on Verticillium wilt progression in potato. ",
  "Disease progression was summarized as the area under the disease progress curve (AUDPC) ",
  "under non-pathogen and Verticillium dahliae-challenged conditions. Values represent mean ± SE ",
  "from biological replicates arranged in an RCBD. Data were analyzed using RCBD ANOVA followed by ",
  "Tukey’s HSD test. Different letters indicate significant differences within the pathogen-challenged group ",
  "(P < 0.05). No statistical letters are shown for the non-pathogen group because AUDPC values were zero across treatments.\n"
)

#############################################################################

############################################################
# STANDALONE SUPPLEMENTARY STATISTICS SCRIPT
# RCBD ANOVA F-values, df, P-values + means, SE, Tukey letters
############################################################

library(tidyverse)
library(multcompView)
library(openxlsx)

# ============================================================
# 1. Treatment settings
# ============================================================

treat_levels <- c(
  "Healthy control",
  "PGPR A only",
  "PGPR B only",
  "PGPR A + B only",
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

with_path_levels <- c(
  "Pathogen only",
  "PGPR A + Pathogen",
  "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

parse_treatment <- function(df) {
  df %>%
    mutate(
      Treatment = factor(as.character(Treatment), levels = treat_levels),
      PathogenGroup = ifelse(
        Treatment %in% with_path_levels,
        "With Pathogen",
        "No Pathogen"
      ),
      PathogenGroup = factor(
        PathogenGroup,
        levels = c("No Pathogen", "With Pathogen")
      ),
      Block = factor(Block)
    )
}

prepare_id_data <- function(df) {
  df %>%
    mutate(
      ID = as.character(ID),
      Block = stringr::str_extract(ID, "R\\s*\\d+"),
      Block = stringr::str_replace_all(Block, "\\s+", ""),
      Block = factor(Block),
      Treatment = factor(as.character(Treatment), levels = treat_levels),
      PathogenGroup = ifelse(
        Treatment %in% with_path_levels,
        "With Pathogen",
        "No Pathogen"
      ),
      PathogenGroup = factor(
        PathogenGroup,
        levels = c("No Pathogen", "With Pathogen")
      )
    )
}

# ============================================================
# 2. Helper functions
# ============================================================

extract_rcbd_anova <- function(dat, response, trait_name, group_name = NA) {
  
  d <- dat %>%
    filter(!is.na(.data[[response]]), !is.na(Treatment), !is.na(Block)) %>%
    mutate(
      Treatment = droplevels(factor(Treatment)),
      Block = droplevels(factor(Block))
    )
  
  if (n_distinct(d$Treatment) < 2 || var(d[[response]], na.rm = TRUE) == 0) {
    return(tibble(
      Trait = trait_name,
      PathogenGroup = group_name,
      Term = "Treatment",
      Df = NA_real_,
      Sum_Sq = NA_real_,
      Mean_Sq = NA_real_,
      F_value = NA_real_,
      P_value = NA_real_,
      Note = "ANOVA not valid: fewer than two treatments or zero variance"
    ))
  }
  
  model <- aov(
    reformulate(c("Treatment", "Block"), response = response),
    data = d
  )
  
  a <- as.data.frame(anova(model))
  a$Term <- rownames(a)
  
  a %>%
    as_tibble() %>%
    filter(Term %in% c("Treatment", "Block", "Residuals")) %>%
    transmute(
      Trait = trait_name,
      PathogenGroup = group_name,
      Term = Term,
      Df = Df,
      Sum_Sq = `Sum Sq`,
      Mean_Sq = `Mean Sq`,
      F_value = `F value`,
      P_value = `Pr(>F)`,
      Note = ""
    )
}

extract_tukey_letters <- function(dat, response, trait_name, group_name = NA) {
  
  d <- dat %>%
    filter(!is.na(.data[[response]]), !is.na(Treatment), !is.na(Block)) %>%
    mutate(
      Treatment = droplevels(factor(Treatment)),
      Block = droplevels(factor(Block))
    )
  
  sum_tab <- d %>%
    group_by(Treatment) %>%
    summarise(
      n = sum(!is.na(.data[[response]])),
      Mean = mean(.data[[response]], na.rm = TRUE),
      SD = sd(.data[[response]], na.rm = TRUE),
      SE = SD / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      Trait = trait_name,
      PathogenGroup = group_name
    )
  
  if (n_distinct(d$Treatment) < 2 || var(d[[response]], na.rm = TRUE) == 0) {
    return(sum_tab %>% mutate(Letters = ""))
  }
  
  model <- aov(
    reformulate(c("Treatment", "Block"), response = response),
    data = d
  )
  
  tuk <- TukeyHSD(model, "Treatment")
  
  pvals <- tuk$Treatment[, "p adj"]
  names(pvals) <- rownames(tuk$Treatment)
  pvals <- pvals[!is.na(pvals)]
  
  if (length(pvals) == 0) {
    return(sum_tab %>% mutate(Letters = ""))
  }
  
  lets <- multcompLetters(pvals)$Letters
  
  letter_tab <- tibble(
    Treatment = names(lets),
    Letters = tolower(gsub("[^A-Za-z]", "", unname(lets)))
  )
  
  sum_tab %>%
    left_join(letter_tab, by = "Treatment") %>%
    mutate(Letters = replace_na(Letters, ""))
}


run_stats_by_group <- function(dat, response, trait_name) {
  
  group_levels <- unique(as.character(dat$PathogenGroup))
  group_levels <- group_levels[!is.na(group_levels)]
  
  anova_out <- map_dfr(group_levels, function(g) {
    d_g <- dat %>% filter(as.character(PathogenGroup) == g)
    extract_rcbd_anova(d_g, response, trait_name, g)
  })
  
  tukey_out <- map_dfr(group_levels, function(g) {
    d_g <- dat %>% filter(as.character(PathogenGroup) == g)
    extract_tukey_letters(d_g, response, trait_name, g)
  })
  
  list(anova = anova_out, tukey = tukey_out)
}



# ============================================================
# 3. Establishment / shoot height / leaf number / SPAD
# ============================================================

shoot_df <- read.delim(
  "shoot.count.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  rename(
    Shoot_20Nov = `20-Nov`,
    Shoot_27Nov = `27-Nov`
  ) %>%
  parse_treatment() %>%
  mutate(
    Shoot_20Nov = as.numeric(Shoot_20Nov),
    Shoot_27Nov = as.numeric(Shoot_27Nov),
    Emerged = ifelse(Shoot_20Nov > 0, 1, 0)
  )

height_dates <- c("20-Nov", "27-Nov", "03-Dec", "10-Dec", "17-Dec", "23-Dec")

height_long <- read.delim(
  "Shoot_length.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  pivot_longer(
    cols = all_of(height_dates),
    names_to = "Date",
    values_to = "Height"
  ) %>%
  parse_treatment() %>%
  mutate(
    Date = factor(Date, levels = height_dates, ordered = TRUE),
    Height = as.numeric(Height)
  )

leaf_long <- read.delim(
  "Leaves_number.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  pivot_longer(
    cols = all_of(height_dates),
    names_to = "Date",
    values_to = "Leaves"
  ) %>%
  parse_treatment() %>%
  mutate(
    Date = factor(Date, levels = height_dates, ordered = TRUE),
    Leaves = as.numeric(Leaves)
  )

chl_dates <- c("03-Dec", "10-Dec", "23-Dec")

chl_long <- read.delim(
  "chlorophyll.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  pivot_longer(
    cols = all_of(chl_dates),
    names_to = "Date",
    values_to = "Chlorophyll"
  ) %>%
  parse_treatment() %>%
  mutate(
    Date = factor(Date, levels = chl_dates, ordered = TRUE),
    Chlorophyll = as.numeric(Chlorophyll)
  )

# ============================================================
# 4. Flowering and tuber data
# ============================================================

inflo_dates <- c("17-Dec", "23-Dec")

inflo_long <- read.delim(
  "inflorescence.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  pivot_longer(
    cols = all_of(inflo_dates),
    names_to = "Date",
    values_to = "Inflorescences"
  ) %>%
  parse_treatment() %>%
  mutate(
    Date = factor(Date, levels = inflo_dates, ordered = TRUE),
    Inflorescences = as.numeric(Inflorescences),
    Flowered = ifelse(Inflorescences > 0, 1, 0)
  )

tuber <- read.delim(
  "Tuber.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  parse_treatment() %>%
  mutate(
    Tubers = as.numeric(Tubers)
  )

# ============================================================
# 5. Fresh biomass, dry biomass, WinRHIZO
# ============================================================

fresh_raw <- read.delim(
  "fresh.biomass.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

names(fresh_raw) <- make.names(names(fresh_raw))

fresh <- prepare_id_data(fresh_raw)

fresh_traits <- c(
  "total_fresh_weight_gm" = "Total fresh biomass",
  "leaves_fresh_weight_gm" = "Leaf fresh weight",
  "steam_fresh_weight_gm" = "Stem fresh weight",
  "leaf_area_cm2" = "Leaf area"
)

fresh_long <- fresh %>%
  pivot_longer(
    cols = all_of(names(fresh_traits)),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    TraitLabel = recode(Trait, !!!fresh_traits)
  ) %>%
  filter(!is.na(Value))

dry_raw <- read.delim(
  "dry.biomass.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

names(dry_raw) <- make.names(names(dry_raw))

dry <- prepare_id_data(dry_raw)

dry_traits <- c(
  "dry_whole_biomass_gm" = "Total dry biomass",
  "leaves_dry_weight_gm" = "Leaf dry weight",
  "steam_dry_weight" = "Stem dry weight",
  "root_dry_weight_gm" = "Root dry weight"
)

dry_long <- dry %>%
  pivot_longer(
    cols = all_of(names(dry_traits)),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    TraitLabel = recode(Trait, !!!dry_traits)
  ) %>%
  filter(!is.na(Value))

root_raw <- read.delim(
  "winrhizodata.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

names(root_raw) <- make.names(names(root_raw))

root <- prepare_id_data(root_raw)

root_traits <- c(
  "TotalLength.cm." = "Total root length",
  "SurfArea.cm2." = "Root surface area",
  "AvgDiam.mm." = "Average root diameter",
  "RootVolume.cm3." = "Root volume"
)

root_long <- root %>%
  pivot_longer(
    cols = all_of(names(root_traits)),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    TraitLabel = recode(Trait, !!!root_traits)
  ) %>%
  filter(!is.na(Value))

# ============================================================
# 6. AUDPC calculation
# ============================================================

audpc_trapz <- function(day, y) {
  
  keep <- !is.na(day) & !is.na(y)
  day <- day[keep]
  y <- y[keep]
  
  if (length(day) < 2) {
    return(NA_real_)
  }
  
  o <- order(day)
  day <- day[o]
  y <- y[o]
  
  sum(
    diff(day) *
      (head(y, -1) + tail(y, -1)) / 2
  )
}

disease_dates <- c("03-Dec", "10-Dec", "17-Dec", "23-Dec")
year_used <- 2025

date_vec <- as.Date(
  paste(year_used, disease_dates),
  format = "%Y %d-%b"
)

day_offsets <- as.numeric(date_vec - min(date_vec))

disease_df <- read.delim(
  "disease scores.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  mutate(
    ID = gsub("\\s+", "", ID),
    Block = factor(Block),
    Treatment = factor(as.character(Treatment), levels = treat_levels),
    PathogenGroup = ifelse(
      Treatment %in% with_path_levels,
      "With Pathogen",
      "No Pathogen"
    ),
    PathogenGroup = factor(
      PathogenGroup,
      levels = c("No Pathogen", "With Pathogen")
    )
  )

disease_long <- disease_df %>%
  pivot_longer(
    cols = all_of(disease_dates),
    names_to = "Date",
    values_to = "Score"
  ) %>%
  mutate(
    Date = factor(Date, levels = disease_dates, ordered = TRUE),
    Day = day_offsets[match(as.character(Date), disease_dates)],
    Score = as.numeric(Score)
  )

audpc_df <- disease_long %>%
  group_by(ID, Treatment, Block, PathogenGroup) %>%
  summarise(
    AUDPC = audpc_trapz(Day, Score),
    .groups = "drop"
  )

write.csv(audpc_df, "AUDPC_PerPot.csv", row.names = FALSE)

# ============================================================
# 7. Run statistics
# ============================================================

est_stats <- list(
  run_stats_by_group(shoot_df, "Emerged", "Emergence"),
  run_stats_by_group(shoot_df, "Shoot_27Nov", "Shoot number per tuber")
)

height_stats <- height_long %>%
  group_by(Date) %>%
  group_split() %>%
  lapply(function(x) {
    run_stats_by_group(x, "Height", paste0("Shoot height_", unique(x$Date)))
  })

leaf_stats <- leaf_long %>%
  group_by(Date) %>%
  group_split() %>%
  lapply(function(x) {
    run_stats_by_group(x, "Leaves", paste0("Leaf number_", unique(x$Date)))
  })

chl_stats <- chl_long %>%
  group_by(Date) %>%
  group_split() %>%
  lapply(function(x) {
    run_stats_by_group(x, "Chlorophyll", paste0("SPAD_", unique(x$Date)))
  })

flower_stats <- inflo_long %>%
  group_by(Date) %>%
  group_split() %>%
  lapply(function(x) {
    run_stats_by_group(x, "Flowered", paste0("Flowering incidence_", unique(x$Date)))
  })

tuber_stats <- list(
  run_stats_by_group(tuber, "Tubers", "Tuber number")
)

fresh_stats <- fresh_long %>%
  group_by(TraitLabel) %>%
  group_split() %>%
  lapply(function(x) {
    run_stats_by_group(x, "Value", unique(x$TraitLabel))
  })

dry_stats <- dry_long %>%
  group_by(TraitLabel) %>%
  group_split() %>%
  lapply(function(x) {
    run_stats_by_group(x, "Value", unique(x$TraitLabel))
  })

root_stats <- root_long %>%
  group_by(TraitLabel) %>%
  group_split() %>%
  lapply(function(x) {
    run_stats_by_group(x, "Value", unique(x$TraitLabel))
  })

audpc_stats <- list(
  run_stats_by_group(audpc_df, "AUDPC", "AUDPC")
)

# ============================================================
# 8. Combine and export
# ============================================================

all_stats <- c(
  est_stats,
  height_stats,
  leaf_stats,
  chl_stats,
  flower_stats,
  tuber_stats,
  fresh_stats,
  dry_stats,
  root_stats,
  audpc_stats
)

anova_supplement <- bind_rows(lapply(all_stats, `[[`, "anova")) %>%
  mutate(
    F_value = round(F_value, 3),
    P_value = signif(P_value, 4)
  )

tukey_supplement <- bind_rows(lapply(all_stats, `[[`, "tukey")) %>%
  mutate(
    Mean = round(Mean, 3),
    SD = round(SD, 3),
    SE = round(SE, 3)
  )

write.csv(
  anova_supplement,
  "Supplementary_Table_ANOVA_Fvalues.csv",
  row.names = FALSE
)

write.csv(
  tukey_supplement,
  "Supplementary_Table_Means_SE_TukeyLetters.csv",
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "ANOVA_F_values")
writeData(wb, "ANOVA_F_values", anova_supplement)

addWorksheet(wb, "Means_SE_Tukey")
writeData(wb, "Means_SE_Tukey", tukey_supplement)

saveWorkbook(
  wb,
  "Supplementary_Statistical_Tables_RCBD_ANOVA_Tukey.xlsx",
  overwrite = TRUE
)

# ============================================================
# 9. Key values for manuscript Results text
# ============================================================

key_traits <- c(
  "Shoot height_23-Dec",
  "Leaf number_23-Dec",
  "SPAD_23-Dec",
  "Tuber number",
  "Total fresh biomass",
  "Total dry biomass",
  "Root dry weight",
  "Leaf area",
  "Total root length",
  "Root surface area",
  "AUDPC"
)

key_anova <- anova_supplement %>%
  filter(
    Trait %in% key_traits,
    Term == "Treatment"
  ) %>%
  arrange(Trait, PathogenGroup)

write.csv(
  key_anova,
  "Key_ANOVA_Values_for_Manuscript_Text.csv",
  row.names = FALSE
)

print(key_anova)

cat("\nDone. Files created:\n")
cat("1. Supplementary_Table_ANOVA_Fvalues.csv\n")
cat("2. Supplementary_Table_Means_SE_TukeyLetters.csv\n")
cat("3. Supplementary_Statistical_Tables_RCBD_ANOVA_Tukey.xlsx\n")
cat("4. Key_ANOVA_Values_for_Manuscript_Text.csv\n")
cat("5. AUDPC_PerPot.csv\n")

