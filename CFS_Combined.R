# ============================================================
# CFS assay: B18 + B14 combined analysis
# Germination (%) + Radicle length
# Bars = Mean ± 95% CI
# Tukey HSD letters
# One 4-panel supplementary figure
# ============================================================

# ---- 0) Packages ----
req <- c("tidyverse", "car", "multcompView", "patchwork")
to_install <- req[!req %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

library(tidyverse)
library(car)
library(multcompView)
library(patchwork)

# ---- 1) User settings ----
files <- tibble(
  Strain = c("B18", "B14"),
  infile = c("CFS_B18.txt", "CFS_B14.txt"),
  total_seeds = c(25, 20)   # change if needed
)

out_prefix <- "CFS_B18_B14_combined"

# ---- 2) Read and clean function ----
read_cfs <- function(infile, strain, total_seeds) {
  
  read.delim(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(
      Strain = strain,
      Treatment = trimws(Treatment),
      Replication = as.factor(Replication),
      Germination_count = as.numeric(as.character(Germination_count)),
      Radicle_avg_cm = suppressWarnings(as.numeric(as.character(Radicle_avg_cm))),
      
      Product = case_when(
        str_detect(Treatment, regex("^CFS", ignore_case = TRUE)) ~ "CFS",
        TRUE ~ "Control"
      ),
      
      Concentration_num = str_match(
        Treatment,
        "\\((\\s*[0-9]+\\.?[0-9]*)\\s*%\\)"
      )[, 2],
      
      Concentration_num = as.numeric(str_replace_all(Concentration_num, "\\s+", "")),
      Concentration_num = ifelse(is.na(Concentration_num), 0, Concentration_num),
      
      Germination_pct = 100 * Germination_count / total_seeds
    )
}

# ---- 3) Import both datasets ----
df <- pmap_dfr(
  files,
  function(Strain, infile, total_seeds) {
    read_cfs(infile, Strain, total_seeds)
  }
) %>%
  mutate(
    Concentration = factor(Concentration_num, levels = c(0, 2.5, 5, 10)),
    Product = factor(Product, levels = c("Control", "CFS")),
    Strain = factor(Strain, levels = c("B18", "B14")),
    Trt = interaction(Product, Concentration, sep = ":", drop = TRUE)
  )

cat("\nParsed treatment levels:\n")
print(df %>% count(Strain, Product, Concentration, Trt))

# ---- 4) Summary function: mean, SD, SE, 95% CI ----
summary_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  mean_x <- ifelse(n == 0, NA_real_, mean(x))
  sd_x   <- ifelse(n < 2, NA_real_, sd(x))
  se_x   <- ifelse(n < 2, NA_real_, sd_x / sqrt(n))
  ci_x   <- ifelse(n < 2, 0, qt(0.975, df = n - 1) * se_x)
  
  tibble(
    n = n,
    mean = mean_x,
    sd = sd_x,
    se = se_x,
    ci95 = ci_x
  )
}

# ---- 5) Tukey letters function ----
get_letters <- function(data, response) {
  
  data <- data %>% filter(is.finite(.data[[response]]))
  
  form <- as.formula(paste(response, "~ Trt"))
  model <- aov(form, data = data)
  
  cat("\n==============================\n")
  cat("ANOVA for:", unique(data$Strain), "-", response, "\n")
  print(summary(model))
  
  cat("\nLevene test:\n")
  print(car::leveneTest(form, data = data))
  
  tuk <- TukeyHSD(model)
  cat("\nTukey HSD:\n")
  print(tuk)
  
  letters <- multcompLetters4(model, tuk)$Trt$Letters
  
  list(
    model = model,
    tukey = tuk,
    letters = letters
  )
}

# ---- 6) Make summary for one strain ----
make_summary <- function(data, response) {
  
  out <- get_letters(data, response)
  
  data %>%
    filter(is.finite(.data[[response]])) %>%
    group_by(Strain, Product, Concentration, Trt) %>%
    summarise(summary_ci(.data[[response]]), .groups = "drop") %>%
    mutate(
      letters = out$letters[as.character(Trt)],
      response = response
    )
}

# ---- 7) Analyze B18 and B14 separately ----
sum_germ <- df %>%
  split(.$Strain) %>%
  map_dfr(~ make_summary(.x, "Germination_pct"))

df_rad <- df %>% filter(is.finite(Radicle_avg_cm))

sum_rad <- df_rad %>%
  split(.$Strain) %>%
  map_dfr(~ make_summary(.x, "Radicle_avg_cm"))

# ---- 8) Prepare y positions ----
sum_germ <- sum_germ %>%
  mutate(
    ymin = pmax(0, mean - ci95),
    ymax = pmin(100, mean + ci95),
    y_lab = pmin(108, ymax + 3)
  )

sum_rad <- sum_rad %>%
  mutate(
    ymin = pmax(0, mean - ci95),
    ymax = mean + ci95,
    y_lab = ymax + 0.25
  )

# ---- 9) Plot function ----
plot_bar <- function(data, ylab, title_text, ylim_top = NULL, show_legend = TRUE) {
  
  p <- ggplot(data, aes(x = Concentration, y = mean, fill = Product)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      linewidth = 0.4
    ) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      position = position_dodge(width = 0.8),
      width = 0.2,
      linewidth = 0.5
    ) +
    geom_text(
      aes(label = letters, y = y_lab),
      position = position_dodge(width = 0.8),
      size = 5
    ) +
    labs(
      title = title_text,
      x = "Concentration (%)",
      y = ylab,
      fill = "Product"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12)
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (!is.null(ylim_top)) {
    p <- p + coord_cartesian(ylim = c(0, ylim_top))
  }
  
  p
}

# ---- 10) Create four panels ----
p1 <- plot_bar(
  filter(sum_germ, Strain == "B18"),
  ylab = "Germination (%)",
  title_text = "A. B18 germination (Mean ± 95% CI)",
  ylim_top = 110,
  show_legend = FALSE
)

p2 <- plot_bar(
  filter(sum_rad, Strain == "B18"),
  ylab = "Radicle length (cm)",
  title_text = "B. B18 radicle length (Mean ± 95% CI)",
  show_legend = TRUE
)

p3 <- plot_bar(
  filter(sum_germ, Strain == "B14"),
  ylab = "Germination (%)",
  title_text = "C. B14 germination (Mean ± 95% CI)",
  ylim_top = 110,
  show_legend = FALSE
)

p4 <- plot_bar(
  filter(sum_rad, Strain == "B14"),
  ylab = "Radicle length (cm)",
  title_text = "D. B14 radicle length (Mean ± 95% CI)",
  show_legend = TRUE
)

# ---- 11) Patch all figures into one panel ----
combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

print(combined_plot)

# ---- 12) Save outputs ----
ggsave(
  paste0(out_prefix, "_4panel_Mean_CI_Tukey.png"),
  combined_plot,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  paste0(out_prefix, "_4panel_Mean_CI_Tukey.pdf"),
  combined_plot,
  width = 14,
  height = 10
)

write.csv(sum_germ, paste0(out_prefix, "_Germination_summary.csv"), row.names = FALSE)
write.csv(sum_rad, paste0(out_prefix, "_Radicle_summary.csv"), row.names = FALSE)

cat("\nSaved files:\n")
cat(paste0(out_prefix, "_4panel_Mean_CI_Tukey.png\n"))
cat(paste0(out_prefix, "_4panel_Mean_CI_Tukey.pdf\n"))
cat(paste0(out_prefix, "_Germination_summary.csv\n"))
cat(paste0(out_prefix, "_Radicle_summary.csv\n"))


# ============================================================
# ADDITIONAL SECTION: Export ANOVA summaries for CFS assays
# Germination (%) and Radicle length
# ============================================================

extract_cfs_anova <- function(data, response_name, response_label) {
  
  data <- data %>%
    filter(is.finite(.data[[response_name]])) %>%
    mutate(
      Strain = factor(Strain),
      Trt = factor(Trt)
    )
  
  # Run separately for each strain
  anova_out <- data %>%
    split(.$Strain) %>%
    map_dfr(function(d) {
      
      strain_name <- unique(as.character(d$Strain))
      
      model <- aov(as.formula(paste(response_name, "~ Trt")), data = d)
      aov_tab <- summary(model)[[1]]
      
      tibble(
        Strain = strain_name,
        Response = response_label,
        Term = rownames(aov_tab),
        Df = aov_tab[, "Df"],
        Sum_Sq = aov_tab[, "Sum Sq"],
        Mean_Sq = aov_tab[, "Mean Sq"],
        F_value = aov_tab[, "F value"],
        P_value = aov_tab[, "Pr(>F)"]
      )
    }) %>%
    mutate(
      F_value = round(F_value, 3),
      P_value = signif(P_value, 4)
    )
  
  return(anova_out)
}

extract_cfs_levene <- function(data, response_name, response_label) {
  
  data <- data %>%
    filter(is.finite(.data[[response_name]])) %>%
    mutate(
      Strain = factor(Strain),
      Trt = factor(Trt)
    )
  
  levene_out <- data %>%
    split(.$Strain) %>%
    map_dfr(function(d) {
      
      strain_name <- unique(as.character(d$Strain))
      
      form <- as.formula(paste(response_name, "~ Trt"))
      lev <- car::leveneTest(form, data = d)
      lev_df <- as.data.frame(lev)
      lev_df$Term <- rownames(lev_df)
      
      tibble(
        Strain = strain_name,
        Response = response_label,
        Term = lev_df$Term,
        Df = lev_df$Df,
        F_value = round(lev_df$`F value`, 3),
        P_value = signif(lev_df$`Pr(>F)`, 4)
      )
    })
  
  return(levene_out)
}

# ---- ANOVA summaries ----
anova_cfs_germ <- extract_cfs_anova(
  df,
  response_name = "Germination_pct",
  response_label = "Germination percentage"
)

anova_cfs_radicle <- extract_cfs_anova(
  df_rad,
  response_name = "Radicle_avg_cm",
  response_label = "Radicle length"
)

anova_cfs_all <- bind_rows(
  anova_cfs_germ,
  anova_cfs_radicle
)

# ---- Levene summaries ----
levene_cfs_germ <- extract_cfs_levene(
  df,
  response_name = "Germination_pct",
  response_label = "Germination percentage"
)

levene_cfs_radicle <- extract_cfs_levene(
  df_rad,
  response_name = "Radicle_avg_cm",
  response_label = "Radicle length"
)

levene_cfs_all <- bind_rows(
  levene_cfs_germ,
  levene_cfs_radicle
)

# ---- Export CSV files ----
write.csv(
  anova_cfs_all,
  paste0(out_prefix, "_ANOVA_summary.csv"),
  row.names = FALSE
)

write.csv(
  levene_cfs_all,
  paste0(out_prefix, "_Levene_summary.csv"),
  row.names = FALSE
)

# ---- Print key Treatment effects only ----
key_cfs_anova <- anova_cfs_all %>%
  filter(Term == "Trt") %>%
  select(Strain, Response, Df, F_value, P_value)

write.csv(
  key_cfs_anova,
  paste0(out_prefix, "_Key_ANOVA_for_Manuscript.csv"),
  row.names = FALSE
)

cat("\nCFS ANOVA summaries saved:\n")
cat(paste0(out_prefix, "_ANOVA_summary.csv\n"))
cat(paste0(out_prefix, "_Levene_summary.csv\n"))
cat(paste0(out_prefix, "_Key_ANOVA_for_Manuscript.csv\n"))

cat("\nKey ANOVA values for manuscript:\n")
print(key_cfs_anova)




