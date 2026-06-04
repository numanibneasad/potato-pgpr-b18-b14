# ============================================================
# Combined seed germination analysis: B14 + B18
# Germination percentage + Maguire Germination Index (MGI)
# Mean ± 95% CI + Tukey letters
#
# Figure labels:
# - Strain labels: P. sp. strain_B18 / P. sp. strain_B14
# - Legend title: CFU/mL
# - Control shown simply as "Control"
# - Caption explains Control = sterile 10 mM MgSO4 only
# ============================================================

# ---- 0) Packages ----
req <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "stringr",
  "janitor", "patchwork", "multcompView", "car", "purrr"
)

to_install <- req[!req %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(janitor)
library(patchwork)
library(multcompView)
library(car)
library(purrr)

# ---- 1) Settings ----
out_dir <- "seed_germination_combined_output"
dir.create(out_dir, showWarnings = FALSE)

N_SEEDS <- 25

files <- tibble(
  file = c("seed_B14.txt", "seed_ B18.txt"),
  strain = c("P. sp. strain_B14", "P. sp. strain_B18")
)

cfu_levels <- c("Control", "10^5", "10^6", "10^7")

cfu_labels <- c(
  "Control" = "Control",
  "10^5" = expression(10^5),
  "10^6" = expression(10^6),
  "10^7" = expression(10^7)
)

# ---- 2) Read + clean function ----
read_seed_file <- function(file, strain_name) {
  
  read_tsv(file, show_col_types = FALSE) %>%
    clean_names() %>%
    select(any_of(c("treatment", "concentration", "day_1", "day_2", "day_3"))) %>%
    mutate(
      strain = strain_name,
      treatment = str_squish(as.character(treatment)),
      concentration = str_squish(as.character(concentration)),
      day_1 = as.numeric(day_1),
      day_2 = as.numeric(day_2),
      day_3 = as.numeric(day_3)
    ) %>%
    filter(!is.na(day_1), !is.na(day_2), !is.na(day_3)) %>%
    mutate(
      cfu_ml = case_when(
        concentration %in% c("0", "00") ~ "Control",
        concentration == "OD1" ~ "10^5",
        concentration == "OD2" ~ "10^6",
        concentration == "OD3" ~ "10^7",
        TRUE ~ concentration
      ),
      cfu_ml = factor(cfu_ml, levels = cfu_levels),
      strain = factor(
        strain,
        levels = c("P. sp. strain_B18", "P. sp. strain_B14")
      )
    )
}

# ---- 3) Combine datasets ----
df <- map2_dfr(files$file, files$strain, read_seed_file)

cat("\nParsed data:\n")
print(df %>% count(strain, cfu_ml))

# ---- 4) Calculate germination % and MGI ----
df_calc <- df %>%
  mutate(
    germ_pct_d1 = (day_1 / N_SEEDS) * 100,
    germ_pct_d2 = (day_2 / N_SEEDS) * 100,
    germ_pct_d3 = (day_3 / N_SEEDS) * 100,
    
    n1 = pmax(0, day_1),
    n2 = pmax(0, day_2 - day_1),
    n3 = pmax(0, day_3 - day_2),
    
    maguire_gri = (n1 / 1) + (n2 / 2) + (n3 / 3)
  )

# ---- 5) Long germination data ----
germ_long <- df_calc %>%
  select(strain, cfu_ml, germ_pct_d1, germ_pct_d2, germ_pct_d3) %>%
  pivot_longer(
    cols = starts_with("germ_pct"),
    names_to = "day",
    values_to = "germination_pct"
  ) %>%
  mutate(
    day = factor(
      day,
      levels = c("germ_pct_d1", "germ_pct_d2", "germ_pct_d3"),
      labels = c("Day 1", "Day 2", "Day 3")
    )
  )

# ---- 6) Summary function: mean ± 95% CI ----
summary_ci <- function(data, response) {
  data %>%
    summarise(
      n = n(),
      mean = mean(.data[[response]], na.rm = TRUE),
      sd = sd(.data[[response]], na.rm = TRUE),
      se = sd / sqrt(n),
      ci95 = qt(0.975, df = n - 1) * se,
      lower = mean - ci95,
      upper = mean + ci95,
      .groups = "drop"
    )
}

sum_germ <- germ_long %>%
  group_by(strain, cfu_ml, day) %>%
  summary_ci("germination_pct") %>%
  mutate(
    lower = pmax(0, lower),
    upper = pmin(100, upper)
  )

sum_gri <- df_calc %>%
  group_by(strain, cfu_ml) %>%
  summary_ci("maguire_gri")

# ---- 7) Tukey letters for germination % within each strain and day ----
get_germ_letters <- function(data) {
  
  fit <- aov(germination_pct ~ cfu_ml, data = data)
  tuk <- TukeyHSD(fit)
  letters <- multcompLetters4(fit, tuk)$cfu_ml$Letters
  
  tibble(
    cfu_ml = names(letters),
    letters = unname(letters)
  )
}

letters_germ <- germ_long %>%
  group_by(strain, day) %>%
  group_modify(~ get_germ_letters(.x)) %>%
  ungroup() %>%
  mutate(cfu_ml = factor(cfu_ml, levels = cfu_levels))

sum_germ <- sum_germ %>%
  left_join(letters_germ, by = c("strain", "day", "cfu_ml")) %>%
  group_by(strain) %>%
  mutate(
    y_letters = pmin(104, upper + 4)
  ) %>%
  ungroup()


# ============================================================
# ANOVA for Germination Percentage
# ============================================================

extract_germ_anova <- function(data){
  
  fit <- aov(germination_pct ~ cfu_ml, data = data)
  
  aov_tab <- summary(fit)[[1]]
  
  tibble(
    Strain = unique(as.character(data$strain)),
    Day = unique(as.character(data$day)),
    DF_Treatment = aov_tab["cfu_ml","Df"],
    DF_Residual = aov_tab["Residuals","Df"],
    F_value = round(aov_tab["cfu_ml","F value"], 2),
    P_value = signif(aov_tab["cfu_ml","Pr(>F)"], 3)
  )
}

anova_germ <- germ_long %>%
  group_by(strain, day) %>%
  group_split() %>%
  map_dfr(extract_germ_anova)

print(anova_germ)

write_csv(
  anova_germ,
  file.path(out_dir, "ANOVA_Germination_summary.csv")
)


# ---- 8) Tukey letters for MGI within each strain ----
get_gri_letters <- function(data) {
  
  fit <- aov(maguire_gri ~ cfu_ml, data = data)
  
  cat("\n=============================\n")
  cat("MGI ANOVA for:", unique(data$strain), "\n")
  print(summary(fit))
  
  cat("\nLevene test:\n")
  print(car::leveneTest(maguire_gri ~ cfu_ml, data = data))
  
  tuk <- TukeyHSD(fit)
  print(tuk)
  
  letters <- multcompLetters4(fit, tuk)$cfu_ml$Letters
  
  tibble(
    cfu_ml = names(letters),
    letters = unname(letters)
  )
}

letters_gri <- df_calc %>%
  group_by(strain) %>%
  group_modify(~ get_gri_letters(.x)) %>%
  ungroup() %>%
  mutate(cfu_ml = factor(cfu_ml, levels = cfu_levels))



# ============================================================
# Extract ANOVA statistics for MGI
# ============================================================

extract_mgi_anova <- function(data){
  
  fit <- aov(maguire_gri ~ cfu_ml, data = data)
  
  aov_tab <- summary(fit)[[1]]
  
  tibble(
    Strain = unique(as.character(data$strain)),
    Response = "Maguire Germination Index",
    DF_Treatment = aov_tab["cfu_ml","Df"],
    DF_Residual = aov_tab["Residuals","Df"],
    F_value = round(aov_tab["cfu_ml","F value"], 2),
    P_value = signif(aov_tab["cfu_ml","Pr(>F)"], 3)
  )
}

anova_mgi <- df_calc %>%
  group_by(strain) %>%
  group_modify(~{
    
    fit <- aov(maguire_gri ~ cfu_ml, data = .x)
    
    aov_tab <- summary(fit)[[1]]
    
    tibble(
      DF_Treatment = aov_tab["cfu_ml","Df"],
      DF_Residual = aov_tab["Residuals","Df"],
      F_value = round(aov_tab["cfu_ml","F value"],2),
      P_value = signif(aov_tab["cfu_ml","Pr(>F)"],3)
    )
    
  }) %>%
  ungroup()

anova_mgi


anova_mgi <- df_calc %>%
  group_by(strain) %>%
  group_split() %>%
  map_dfr(extract_mgi_anova)

print(anova_mgi)

write_csv(
  anova_mgi,
  file.path(out_dir, "ANOVA_MGI_summary.csv")
)

sum_gri <- sum_gri %>%
  left_join(letters_gri, by = c("strain", "cfu_ml")) %>%
  group_by(strain) %>%
  mutate(
    y_letters = upper + 0.08 * max(upper, na.rm = TRUE)
  ) %>%
  ungroup()

# ---- 9) Publication-style theme ----
theme_pub <- theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 19),
    plot.title = element_text(size = 22, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    plot.margin = margin(8, 8, 8, 8)
  )

# ---- 10) Germination percentage line graph ----
p_germ <- ggplot(
  sum_germ,
  aes(
    x = day,
    y = mean,
    group = cfu_ml,
    color = cfu_ml,
    shape = cfu_ml
  )
) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 4.4) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.12,
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = letters, y = y_letters),
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~ strain, nrow = 1) +
  scale_y_continuous(
    limits = c(0, 106),
    breaks = seq(0, 100, 20)
  ) +
  scale_color_discrete(
    name = "CFU/mL",
    labels = cfu_labels
  ) +
  scale_shape_discrete(
    name = "CFU/mL",
    labels = cfu_labels
  ) +
  labs(
    title = "A. Germination percentage",
    x = "Germination time",
    y = "Germination (%)"
  ) +
  theme_pub

# ---- 11) MGI bar graph ----
p_gri <- ggplot(
  sum_gri,
  aes(x = cfu_ml, y = mean, fill = cfu_ml)
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.20,
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = letters, y = y_letters),
    size = 5.5,
    fontface = "bold"
  ) +
  facet_wrap(~ strain, nrow = 1) +
  scale_fill_discrete(
    name = "CFU/mL",
    labels = cfu_labels
  ) +
  scale_x_discrete(
    labels = cfu_labels
  ) +
  labs(
    title = "B. Maguire germination index",
    x = "CFU/mL",
    y = "MGI"
  ) +
  theme_pub

# ---- 12) Combine plots ----
combined_plot <- p_germ / p_gri +
  plot_layout(heights = c(1.15, 1), guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)

# ---- 13) Save outputs ----
ggsave(
  file.path(out_dir, "combined_germination_MGI_B14_B18.png"),
  combined_plot,
  width = 15,
  height = 11,
  dpi = 300
)

ggsave(
  file.path(out_dir, "combined_germination_MGI_B14_B18.tiff"),
  combined_plot,
  width = 15,
  height = 11,
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(out_dir, "combined_germination_MGI_B14_B18.pdf"),
  combined_plot,
  width = 15,
  height = 11
)

write_csv(sum_germ, file.path(out_dir, "summary_germination_percentage.csv"))
write_csv(sum_gri, file.path(out_dir, "summary_maguire_GRI.csv"))

cat("\nDone. Outputs saved in:", out_dir, "\n")
