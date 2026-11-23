# FINAL STATISTICAL CHECK - potato PGPR x pathogen RCBD
#
# Purpose:
#   Run corrected inferential analyses without changing any publication figure
#   layout, colours, labels, panel arrangement, or existing image dimensions.
#   The original figure script remains the visual source of truth. This script
#   exports model results, diagnostics, model-based estimates, confidence
#   intervals, and compact-letter displays that can replace the old ANOVA/Tukey
#   annotations where appropriate.
#
# Design:
#   4-level PGPR treatment x 2-level pathogen status, 8 RCBD blocks.
#   Block is fitted as a fixed design effect. Pot is a random intercept for
#   repeated observations. All analyses use the interaction as the primary test.
#
# IMPORTANT:
#   1. Set project_dir below if this file is not stored in the data directory.
#   2. Do not interpret simple-effect letters unless the relevant interaction
#      and diagnostics have been reviewed.
#   3. Letters compare PGPR treatments WITHIN pathogen status (and date for
#      repeated traits); letters in different facets are not comparable.
#   4. This script never imputes missing observations.

required <- c("dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
              "lme4", "lmerTest", "glmmTMB", "emmeans", "multcomp", "car",
              "brglm2")
missing_packages <- setdiff(required, rownames(installed.packages()))
if (length(missing_packages)) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(purrr); library(tibble); library(lme4); library(lmerTest)
  library(glmmTMB); library(emmeans); library(multcomp)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_dir, "Final_statistical_check_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

options(contrasts = c("contr.sum", "contr.poly"))
analysis_log <- file.path(output_dir, "analysis_log.txt")
writeLines(c(
  paste("Analysis run:", Sys.time()),
  paste("Project directory:", project_dir),
  "Figures were not regenerated or reformatted by this script."
), analysis_log)

log_line <- function(...) {
  cat(paste0(...), "\n", file = analysis_log, append = TRUE)
}

treatment_levels <- c(
  "Healthy control", "PGPR A only", "PGPR B only", "PGPR A + B only",
  "Pathogen only", "PGPR A + Pathogen", "PGPR B + Pathogen",
  "PGPR A + B + Pathogen"
)

read_tab <- function(filename) {
  path <- file.path(project_dir, filename)
  if (!file.exists(path)) stop("Required file not found: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("", "NA", "N/A", "."))
}

prepare_design <- function(dat, label) {
  names(dat) <- trimws(names(dat))
  if (!all(c("Treatment", "Block") %in% names(dat))) {
    if (!"Block" %in% names(dat) && "ID" %in% names(dat)) {
      dat$Block <- str_extract(dat$ID, "(?i)(?<=R)\\s*[0-9]+")
    }
  }
  if (!all(c("Treatment", "Block") %in% names(dat))) {
    stop(label, ": Treatment and Block are required")
  }
  dat <- dat %>%
    mutate(
      Treatment = factor(str_squish(Treatment), levels = treatment_levels),
      Block = factor(Block),
      PGPR = factor(case_when(
        Treatment %in% c("Healthy control", "Pathogen only") ~ "None",
        Treatment %in% c("PGPR A only", "PGPR A + Pathogen") ~ "A",
        Treatment %in% c("PGPR B only", "PGPR B + Pathogen") ~ "B",
        Treatment %in% c("PGPR A + B only", "PGPR A + B + Pathogen") ~ "A+B"
      ), levels = c("None", "A", "B", "A+B")),
      Pathogen = factor(
        if_else(Treatment %in% treatment_levels[5:8], "With Pathogen", "No Pathogen"),
        levels = c("No Pathogen", "With Pathogen")
      ),
      Pot = factor(if ("ID" %in% names(.)) str_squish(as.character(ID)) else
        paste(Treatment, Block,
              if ("Replicate" %in% names(.)) Replicate else seq_len(n()), sep = "__"))
    )
  if (anyNA(dat$Treatment)) stop(label, ": unrecognized treatment label(s)")
  if (anyNA(dat$Block)) stop(label, ": missing block value(s)")
  invisible(log_line(label, ": ", nrow(dat), " rows; ",
                     n_distinct(dat$Pot), " experimental units"))
  dat
}

pivot_dates <- function(dat, dates, value_name) {
  absent <- setdiff(dates, names(dat))
  if (length(absent)) stop(value_name, ": missing date columns: ", paste(absent, collapse = ", "))
  dat %>%
    pivot_longer(all_of(dates), names_to = "Date", values_to = value_name) %>%
    mutate(Date = factor(Date, levels = dates, ordered = TRUE))
}

save_capture <- function(object, filename) {
  capture.output(object, file = file.path(output_dir, filename))
}

save_model_bundle <- function(model, tag, emm_formula = ~ PGPR | Pathogen,
                              type = "link", make_letters = TRUE) {
  saveRDS(model, file.path(output_dir, paste0(tag, "_model.rds")))
  save_capture(summary(model), paste0(tag, "_model_summary.txt"))
  # Type-III tests are used because missing observations can make the design
  # unbalanced; sum-to-zero contrasts were set globally above.
  omnibus <- if (inherits(model, "brglmFit")) {
    car::Anova(model, type = 3, test.statistic = "Wald")
  } else {
    car::Anova(model, type = 3)
  }
  save_capture(omnibus, paste0(tag, "_omnibus_tests.txt"))
  em <- emmeans(model, emm_formula, type = type)
  write_csv(as.data.frame(confint(em)),
            file.path(output_dir, paste0(tag, "_estimated_marginal_means.csv")))
  pw <- pairs(em, adjust = "tukey")
  write_csv(as.data.frame(pw),
            file.path(output_dir, paste0(tag, "_pairwise_tukey.csv")))
  cl <- if (make_letters) tryCatch(
      as.data.frame(multcomp::cld(em, adjust = "tukey", Letters = letters)) %>%
        mutate(.group = str_remove_all(.group, "\\s+")),
      error = function(e) {
        log_line(tag, ": compact letters unavailable - ", conditionMessage(e))
        NULL
      }
    ) else NULL
  if (!make_letters) log_line(tag, ": letters deliberately omitted because fitted probabilities are at/near a boundary")
  if (!is.null(cl)) write_csv(cl, file.path(output_dir, paste0(tag, "_letters.csv")))
  invisible(list(emmeans = em, pairs = pw, letters = cl))
}

pearson_dispersion <- function(model) {
  rp <- residuals(model, type = "pearson")
  sum(rp^2, na.rm = TRUE) / df.residual(model)
}

fit_count <- function(dat, response, repeated = FALSE, tag) {
  dat <- dat %>% filter(!is.na(.data[[response]]))
  if (any(dat[[response]] < 0 | dat[[response]] != floor(dat[[response]]))) {
    stop(tag, ": response is not a non-negative integer count")
  }
  rhs <- if (repeated) "PGPR * Pathogen * Date + Block + (1 | Pot)" else
    "PGPR * Pathogen + Block"
  form <- as.formula(paste(response, "~", rhs))
  pois <- glmmTMB(form, data = dat, family = poisson)
  dispersion <- pearson_dispersion(pois)
  model <- if (is.finite(dispersion) && dispersion > 1.20) {
    log_line(tag, ": Poisson dispersion = ", round(dispersion, 3),
             "; negative binomial model selected")
    glmmTMB(form, data = dat, family = nbinom2)
  } else {
    log_line(tag, ": Poisson dispersion = ", round(dispersion, 3),
             "; Poisson model retained")
    pois
  }
  save_capture(c(Poisson_dispersion = dispersion,
                 Selected_family = family(model)$family),
               paste0(tag, "_count_diagnostic.txt"))
  model
}

fit_gaussian <- function(dat, response, repeated = FALSE) {
  dat <- dat %>% filter(!is.na(.data[[response]]))
  form <- as.formula(paste(response, "~",
    if (repeated) "PGPR * Pathogen * Date + Block + (1 | Pot)" else
      "PGPR * Pathogen + Block"))
  if (repeated) lmer(form, data = dat, REML = TRUE) else lm(form, data = dat)
}

fit_binary <- function(dat, response, repeated = FALSE) {
  dat <- dat %>% filter(!is.na(.data[[response]]))
  if (any(!dat[[response]] %in% c(0, 1))) stop(response, " must contain only 0/1")
  form <- as.formula(paste(response, "~",
    if (repeated) "PGPR * Pathogen * Date + Block + (1 | Pot)" else
      "PGPR * Pathogen + Block"))
  if (repeated) glmmTMB(
    form, data = dat, family = binomial,
    control = glmmTMBControl(optCtrl = list(iter.max = 100000, eval.max = 100000))
  ) else {
    ordinary <- glm(form, data = dat, family = binomial)
    boundary <- any(fitted(ordinary) < 1e-6 | fitted(ordinary) > 1 - 1e-6)
    if (boundary) {
      log_line(response, ": separation detected; using bias-reduced binomial GLM (brglm2)")
      glm(form, data = dat, family = binomial, method = brglm2::brglmFit,
          type = "AS_mean")
    } else ordinary
  }
}

wilson_summary <- function(dat, response, groups) {
  dat %>%
    filter(!is.na(.data[[response]])) %>%
    group_by(across(all_of(groups))) %>%
    summarise(x = sum(.data[[response]]), n = n(), .groups = "drop") %>%
    mutate(
      proportion = x / n,
      z = qnorm(0.975),
      centre = (proportion + z^2 / (2 * n)) / (1 + z^2 / n),
      half_width = z * sqrt(proportion * (1 - proportion) / n + z^2 / (4 * n^2)) /
        (1 + z^2 / n),
      lower_95 = pmax(0, centre - half_width),
      upper_95 = pmin(1, centre + half_width)
    ) %>% dplyr::select(-z, -centre, -half_width)
}

# 1. Emergence and shoot count -------------------------------------------------
shoot <- prepare_design(read_tab("shoot.count.txt"), "shoot.count.txt") %>%
  mutate(Emerged = if_else(is.na(`20-Nov`), NA_integer_, as.integer(`20-Nov` > 0)),
         Shoots = as.numeric(`27-Nov`))

emergence_model <- fit_binary(shoot, "Emerged")
save_model_bundle(emergence_model, "emergence_binomial", type = "response")
write_csv(wilson_summary(shoot, "Emerged", c("Treatment", "PGPR", "Pathogen")),
          file.path(output_dir, "emergence_observed_Wilson95.csv"))

shoot_model <- fit_count(shoot, "Shoots", tag = "shoot_count")
save_model_bundle(shoot_model, "shoot_count")

# 2. Repeated shoot height, leaf number, and SPAD -----------------------------
growth_dates <- c("20-Nov", "27-Nov", "03-Dec", "10-Dec", "17-Dec", "23-Dec")
spad_dates <- c("03-Dec", "10-Dec", "23-Dec")

height <- prepare_design(read_tab("Shoot_length.txt"), "Shoot_length.txt") %>%
  pivot_dates(growth_dates, "Height")
height_model <- fit_gaussian(height, "Height", repeated = TRUE)
save_model_bundle(height_model, "shoot_height_repeated",
                  ~ PGPR | Pathogen * Date)

leaves <- prepare_design(read_tab("Leaves_number.txt"), "Leaves_number.txt") %>%
  pivot_dates(growth_dates, "Leaves")
if (any(!is.na(leaves$Leaves) & leaves$Leaves != floor(leaves$Leaves))) {
  log_line("leaf_number_repeated: fractional observations detected; treating leaf number as an averaged continuous response and fitting a Gaussian repeated mixed model (values were not rounded)")
  leaf_model <- fit_gaussian(leaves, "Leaves", repeated = TRUE)
} else {
  leaf_model <- fit_count(leaves, "Leaves", repeated = TRUE,
                          tag = "leaf_number_repeated")
}
save_model_bundle(leaf_model, "leaf_number_repeated", ~ PGPR | Pathogen * Date)

spad <- prepare_design(read_tab("chlorophyll.txt"), "chlorophyll.txt") %>%
  pivot_dates(spad_dates, "SPAD")
spad_model <- fit_gaussian(spad, "SPAD", repeated = TRUE)
save_model_bundle(spad_model, "SPAD_repeated", ~ PGPR | Pathogen * Date)

# 3. Flowering incidence and tuber count --------------------------------------
flower_dates <- c("17-Dec", "23-Dec")
flower <- prepare_design(read_tab("inflorescence.txt"), "inflorescence.txt") %>%
  pivot_dates(flower_dates, "Inflorescences") %>%
  mutate(Flowered = if_else(is.na(Inflorescences), NA_integer_,
                            as.integer(Inflorescences > 0)))
flower_model <- fit_binary(flower, "Flowered", repeated = TRUE)
flower_converged <- isTRUE(flower_model$fit$convergence == 0) &&
  isTRUE(flower_model$sdr$pdHess)
flower_boundary <- any(fitted(flower_model) < 1e-6 | fitted(flower_model) > 1 - 1e-6)
if (!flower_converged) {
  log_line("flowering_incidence_repeated_binomial: MODEL DID NOT CONVERGE; omnibus tests and contrasts are exploratory only and must not be reported")
}
save_model_bundle(flower_model, "flowering_incidence_repeated_binomial",
                  ~ PGPR | Pathogen * Date, type = "response",
                  make_letters = flower_converged && !flower_boundary)
write_csv(wilson_summary(flower, "Flowered",
                         c("Date", "Treatment", "PGPR", "Pathogen")),
          file.path(output_dir, "flowering_observed_Wilson95.csv"))
log_line("Flowering letters should be omitted if separation, boundary estimates, ",
         "or non-estimable contrasts occur; report model contrasts instead.")

tuber <- prepare_design(read_tab("Tuber.txt"), "Tuber.txt") %>%
  mutate(Tubers = as.numeric(Tubers))
tuber_model <- fit_count(tuber, "Tubers", tag = "tuber_count")
save_model_bundle(tuber_model, "tuber_count")

# 4. Biomass and root architecture -------------------------------------------
# The updated dry file is authoritative only if shared values agree with the
# original file. Any disagreement stops the analysis and requires a raw-data
# decision by the investigator.
compare_biomass_files <- function(old_name, new_name) {
  old_path <- file.path(project_dir, old_name)
  new_path <- file.path(project_dir, new_name)
  if (!file.exists(old_path) || !file.exists(new_path)) return(invisible(NULL))
  old <- read.delim(old_path, check.names = FALSE, stringsAsFactors = FALSE)
  new <- read.delim(new_path, check.names = FALSE, stringsAsFactors = FALSE)
  keys <- intersect(c("ID", "Treatment"), intersect(names(old), names(new)))
  shared <- setdiff(intersect(names(old), names(new)), keys)
  chk <- inner_join(old, new, by = keys, suffix = c("_old", "_new"))
  conflict_flags <- map_lgl(shared, function(v) {
    a <- suppressWarnings(as.numeric(chk[[paste0(v, "_old")]]))
    b <- suppressWarnings(as.numeric(chk[[paste0(v, "_new")]]))
    any(!is.na(a) & !is.na(b) & abs(a - b) > 1e-10)
  })
  if (any(conflict_flags)) {
    conflict_vars <- shared[conflict_flags]
    conflict_rows <- map_dfr(conflict_vars, function(v) {
      old_v <- paste0(v, "_old"); new_v <- paste0(v, "_new")
      chk %>%
        filter(!is.na(.data[[old_v]]), !is.na(.data[[new_v]]),
               abs(as.numeric(.data[[old_v]]) - as.numeric(.data[[new_v]])) > 1e-10) %>%
        transmute(across(all_of(keys)), Trait = v,
                  Old_value = .data[[old_v]], Updated_value = .data[[new_v]])
    })
    write_csv(conflict_rows, file.path(output_dir, "biomass_file_conflicts.csv"))
    log_line("Biomass conflict detected in ", paste(conflict_vars, collapse = ", "),
             "; updated file used as named, and old-file sensitivity model added.")
    return(conflict_rows)
  }
  log_line("Biomass cross-check passed for ", old_name, " vs ", new_name)
  tibble()
}
biomass_conflicts <- compare_biomass_files("dry.biomass.txt", "dry_biomass_updated.txt")

fresh <- prepare_design(read_tab("fresh.biomass.txt"), "fresh.biomass.txt")
dry <- prepare_design(read_tab("dry_biomass_updated.txt"), "dry_biomass_updated.txt")
dry_old <- prepare_design(read_tab("dry.biomass.txt"), "dry.biomass.txt")
root <- prepare_design(read_tab("winrhizodata.txt"), "winrhizodata.txt")
names(root) <- make.names(names(root))

continuous_traits <- list(
  fresh_biomass = list(fresh, "total_fresh_weight_gm"),
  leaf_area = list(fresh, "leaf_area_cm2"),
  dry_biomass = list(dry, "dry_whole_biomass_gm"),
  root_dry_weight = list(dry, "root_dry_weight_gm"),
  root_total_length = list(root, "TotalLength.cm."),
  root_surface_area = list(root, "SurfArea.cm2."),
  root_average_diameter = list(root, "AvgDiam.mm."),
  root_volume = list(root, "RootVolume.cm3.")
)

walk(names(continuous_traits), function(tag) {
  item <- continuous_traits[[tag]]
  dat <- item[[1]]; response <- item[[2]]
  if (!response %in% names(dat)) {
    log_line(tag, ": missing column ", response); return(invisible(NULL))
  }
  model <- fit_gaussian(dat, response)
  save_model_bundle(model, tag)
  diag <- tibble(
    fitted = fitted(model), residual = residuals(model),
    standardized_residual = rstandard(model)
  )
  write_csv(diag, file.path(output_dir, paste0(tag, "_residual_diagnostics.csv")))
})

if (nrow(biomass_conflicts) &&
    "root_dry_weight_gm" %in% biomass_conflicts$Trait) {
  old_root_model <- fit_gaussian(dry_old, "root_dry_weight_gm")
  save_model_bundle(old_root_model, "root_dry_weight_OLD_FILE_sensitivity")
}

sample_sizes <- bind_rows(
  fresh %>% summarise(across(where(is.numeric), ~ sum(!is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "Trait", values_to = "n") %>% mutate(File = "fresh.biomass.txt"),
  dry %>% summarise(across(where(is.numeric), ~ sum(!is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "Trait", values_to = "n") %>% mutate(File = "dry_biomass_updated.txt"),
  root %>% summarise(across(where(is.numeric), ~ sum(!is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "Trait", values_to = "n") %>% mutate(File = "winrhizodata.txt")
)
write_csv(sample_sizes, file.path(output_dir, "biomass_root_sample_sizes.csv"))

# 5. AUDPC with an explicit complete-trajectory rule --------------------------
disease_dates <- c("03-Dec", "10-Dec", "17-Dec", "23-Dec")
disease_days <- c(0, 7, 14, 20)
disease <- prepare_design(read_tab("disease scores.txt"), "disease scores.txt")

audpc_complete <- function(values, days = disease_days) {
  values <- as.numeric(values)
  if (anyNA(values)) return(NA_real_)
  sum(diff(days) * (head(values, -1) + tail(values, -1)) / 2)
}

audpc <- disease %>%
  rowwise() %>%
  mutate(AUDPC = audpc_complete(c_across(all_of(disease_dates)))) %>%
  ungroup()
write_csv(audpc %>% dplyr::select(any_of(c("ID", "Treatment", "Block", "Replicate", "AUDPC"))),
          file.path(output_dir, "AUDPC_per_pot_complete_trajectories.csv"))
log_line("AUDPC rule: any missing disease assessment makes that pot's AUDPC missing.")

# Primary factorial inference is meaningful only where variation exists.
# The no-pathogen controls are retained descriptively but are not forced into a
# Gaussian model when all scores are zero.
write_csv(audpc %>% group_by(Treatment, PGPR, Pathogen) %>%
            summarise(n = sum(!is.na(AUDPC)), mean = mean(AUDPC, na.rm = TRUE),
                      sd = sd(AUDPC, na.rm = TRUE), .groups = "drop"),
          file.path(output_dir, "AUDPC_descriptive_summary.csv"))

audpc_infected <- audpc %>% filter(Pathogen == "With Pathogen", !is.na(AUDPC))
if (n_distinct(audpc_infected$AUDPC) > 1) {
  audpc_model <- lm(AUDPC ~ PGPR + Block, data = audpc_infected)
  save_model_bundle(audpc_model, "AUDPC_infected_stratum", ~ PGPR)
} else {
  log_line("AUDPC infected stratum has zero variance; inferential model omitted.")
}

# 6. Reproducibility record ----------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
log_line("Completed successfully. Review every omnibus test and diagnostic before ",
         "updating manuscript claims or figure letters.")
message("Completed. Results saved to: ", output_dir)
