library(readr)
library(dplyr)
library(purrr)
library(broom)
library(stringr)
library(tibble)
library(sensemakr)
library(tidyr)

treatment   <- "high_religiosity_75" 
outcome     <- "num_children"      
controls    <- c("age","age_squared","education","income","married","female","urban",
                 "economic_hardship_index","most_people_trusted","left_right_scale",
                 "well_being_index")   
benchmark   <- "married"        
alpha_level <- 0.05

pooled_path  <- "data/sensitivity_data/pooled_religion_data.csv"
country_dir  <- "data/sensitivity_data/countries/religion/"
out_dir      <- "data/sensitivity_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_cols <- c(treatment, outcome, controls)

label_from_file <- function(fp, df) {
  # Try country name, else B_COUNTRY code, else filename stem
  if ("country" %in% names(df) && length(df$country) > 0) {
    as.character(df$country[1])
  } else if ("B_COUNTRY" %in% names(df) && length(df$B_COUNTRY) > 0) {
    as.character(df$B_COUNTRY[1])
  } else {
    tools::file_path_sans_ext(basename(fp))
  }
}

run_sense_on_df <- function(df, unit_label) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    return(list(
      summary = tibble(
        unit = unit_label,
        n = nrow(df),
        coef = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_,
        partial_r2_treat = NA_real_, rv_q1 = NA_real_, rv_q1_alpha05 = NA_real_,
        note = paste0("Missing columns: ", paste(missing, collapse = ", "))
      ),
      sense = NULL
    ))
  }

  work <- df |>
    dplyr::select(dplyr::all_of(c(outcome, treatment, controls))) |>
    tidyr::drop_na()

  tr_counts <- table(work[[treatment]])
  if (length(tr_counts) < 2 || min(tr_counts) < 30) {
    return(list(
      summary = tibble(
        unit = unit_label,
        n = nrow(work),
        coef = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_,
        partial_r2_treat = NA_real_, rv_q1 = NA_real_, rv_q1_alpha05 = NA_real_,
        note = "Insufficient treatment variation after NA drop"
      ),
      sense = NULL
    ))
  }

  fml   <- as.formula(paste(outcome, "~", treatment, "+", paste(controls, collapse = " + ")))
  model <- lm(fml, data = work)

  tid  <- broom::tidy(model)
  rowt <- dplyr::filter(tid, term == treatment)
  if (nrow(rowt) != 1) {
    return(list(
      summary = tibble(
        unit = unit_label,
        n = nrow(work),
        coef = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_,
        partial_r2_treat = NA_real_, rv_q1 = NA_real_, rv_q1_alpha05 = NA_real_,
        note = "Treatment term not found in model"
      ),
      sense = NULL
    ))
  }

  sense <- sensemakr::sensemakr(
    model,
    treatment = treatment,
    benchmark_covariates = benchmark,
    kd = 1, q = 1, reduce = TRUE
  )

  pr2  <- as.numeric(sensemakr::partial_r2(model, covariates = treatment))
  rv1  <- as.numeric(sensemakr::robustness_value(model, covariates = treatment, q = 1))
  rv1a <- as.numeric(sensemakr::robustness_value(model, covariates = treatment, q = 1, alpha = alpha_level))

  list(
    summary = tibble(
      unit  = unit_label,
      n     = nrow(work),
      coef  = rowt$estimate,
      se    = rowt$std.error,
      t     = rowt$statistic,
      p     = rowt$p.value,
      partial_r2_treat = pr2,
      rv_q1           = rv1,
      rv_q1_alpha05   = rv1a,
      note  = NA_character_
    ),
    sense = sense
  )
}

safe_run_file <- function(fp) {
  out <- try({
    df <- readr::read_csv(fp, show_col_types = FALSE)
    lab <- label_from_file(fp, df)
    res <- run_sense_on_df(df, lab)
    res$summary %>% dplyr::mutate(source_file = basename(fp), .before = 1)
  }, silent = TRUE)

  if (inherits(out, "try-error")) {
    tibble(
      source_file = basename(fp),
      unit = tools::file_path_sans_ext(basename(fp)),
      n = NA_integer_,
      coef = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_,
      partial_r2_treat = NA_real_, rv_q1 = NA_real_, rv_q1_alpha05 = NA_real_,
      note = as.character(attr(out, "condition")$message)
    )
  } else out
}


# Pooled run 
pooled_res <- NULL
if (file.exists(pooled_path)) {
  message("Running pooled sensitivity on: ", pooled_path)
  pooled_df  <- readr::read_csv(pooled_path, show_col_types = FALSE)
  pooled_out <- run_sense_on_df(pooled_df, unit_label = "Pooled")

  # summary row
  pooled_res <- pooled_out$summary %>%
    dplyr::mutate(source_file = basename(pooled_path), .before = 1)

  if (!is.null(pooled_out$sense)) {
    pooled_plot_path <- file.path(out_dir, "pooled_sense_plot_religion.png")
    png(pooled_plot_path, width = 800, height = 600)
    plot(pooled_out$sense)
    dev.off()
    message("Saved pooled plot to: ", pooled_plot_path)
  } else {
  message("No sense object for pooled — no plot saved.")
}
} else {
  message("No pooled.csv found at ", pooled_path, " — skipping pooled run.")
}
# Country runs 
country_files <- list.files(country_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(country_files) == 0) {
  warning("No country CSV files found in ", country_dir)
  country_results <- tibble()
} else {
  country_results <- map_dfr(country_files, safe_run_file)
}

# combine and save
final_results <- bind_rows(pooled_res, country_results)

readr::write_csv(final_results, file.path(out_dir, "sensemakr_religion_summary_table.csv"))
print(final_results %>% arrange(unit))


#### MEDIA RUN #####
treatment   <- "heavy_media_use" 
outcome     <- "num_children"      
controls    <- c("age","age_squared","education","income","married","female","urban",
                 "economic_hardship_index","most_people_trusted","left_right_scale", 
                 "work_duty_society", "work_comes_first", "religion_importance", "service_attendance")   
benchmark   <- "married"        
alpha_level <- 0.05

pooled_path  <- "data/sensitivity_data/pooled_media_data.csv"
country_dir  <- "data/sensitivity_data/countries/media/"
out_dir      <- "data/sensitivity_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_cols <- c(treatment, outcome, controls)

# Pooled run 
pooled_res <- NULL
if (file.exists(pooled_path)) {
  message("Running pooled sensitivity on: ", pooled_path)
  pooled_df  <- readr::read_csv(pooled_path, show_col_types = FALSE)
  pooled_out <- run_sense_on_df(pooled_df, unit_label = "Pooled")

  # summary row
  pooled_res <- pooled_out$summary %>%
    dplyr::mutate(source_file = basename(pooled_path), .before = 1)

  if (!is.null(pooled_out$sense)) {
    pooled_plot_path <- file.path(out_dir, "pooled_sense_plot_media.png")
    png(pooled_plot_path, width = 800, height = 600)
    plot(pooled_out$sense)
    dev.off()
    message("Saved pooled plot to: ", pooled_plot_path)
  } else {
  message("No sense object for pooled — no plot saved.")
}
} else {
  message("No pooled.csv found at ", pooled_path, " — skipping pooled run.")
}
# Country runs 
country_files <- list.files(country_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(country_files) == 0) {
  warning("No country CSV files found in ", country_dir)
  country_results <- tibble()
} else {
  country_results <- map_dfr(country_files, safe_run_file)
}

# combine and save
final_results <- bind_rows(pooled_res, country_results)

readr::write_csv(final_results, file.path(out_dir, "sensemakr_media_summary_table.csv"))
print(final_results %>% arrange(unit))


##### GENDER RUN #####
treatment   <- "high_traditional_gender" 
outcome     <- "num_children"      
controls    <- c("age","age_squared","education","income","married","female","urban",
                 "economic_hardship_index","most_people_trusted","left_right_scale", 
                 "well_being_index", "religion_importance", "service_attendance")   
benchmark   <- "married"        
alpha_level <- 0.05

pooled_path  <- "data/sensitivity_data/pooled_gender_data.csv"
country_dir  <- "data/sensitivity_data/countries/gender/"
out_dir      <- "data/sensitivity_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_cols <- c(treatment, outcome, controls)

# Pooled run 
pooled_res <- NULL
if (file.exists(pooled_path)) {
  message("Running pooled sensitivity on: ", pooled_path)
  pooled_df  <- readr::read_csv(pooled_path, show_col_types = FALSE)
  pooled_out <- run_sense_on_df(pooled_df, unit_label = "Pooled")

  # summary row
  pooled_res <- pooled_out$summary %>%
    dplyr::mutate(source_file = basename(pooled_path), .before = 1)

  if (!is.null(pooled_out$sense)) {
    pooled_plot_path <- file.path(out_dir, "pooled_sense_plot_gender.png")
    png(pooled_plot_path, width = 800, height = 600)
    plot(pooled_out$sense)
    dev.off()
    message("Saved pooled plot to: ", pooled_plot_path)
  } else {
  message("No sense object for pooled — no plot saved.")
}
} else {
  message("No pooled.csv found at ", pooled_path, " — skipping pooled run.")
}
# Country runs 
country_files <- list.files(country_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(country_files) == 0) {
  warning("No country CSV files found in ", country_dir)
  country_results <- tibble()
} else {
  country_results <- map_dfr(country_files, safe_run_file)
}

# combine and save
final_results <- bind_rows(pooled_res, country_results)

readr::write_csv(final_results, file.path(out_dir, "sensemakr_gender_summary_table.csv"))


print(final_results %>% arrange(unit))


### HIGH WELL BEING RUN ####
treatment   <- "high_well_being" 
outcome     <- "num_children"      
controls    <- c("age","age_squared","education","income","married","female","urban",
                 "economic_hardship_index","most_people_trusted","left_right_scale", "religion_importance",
                 "service_attendance", "work_duty_society", "work_comes_first", "traditional_gender_index")   
benchmark   <- "married"        
alpha_level <- 0.05

pooled_path  <- "data/sensitivity_data/pooled_happy_data.csv"
country_dir  <- "data/sensitivity_data/countries/happy/"
out_dir      <- "data/sensitivity_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_cols <- c(treatment, outcome, controls)

# Pooled run 
pooled_res <- NULL
if (file.exists(pooled_path)) {
  message("Running pooled sensitivity on: ", pooled_path)
  pooled_df  <- readr::read_csv(pooled_path, show_col_types = FALSE)
  pooled_out <- run_sense_on_df(pooled_df, unit_label = "Pooled")

  # summary row
  pooled_res <- pooled_out$summary %>%
    dplyr::mutate(source_file = basename(pooled_path), .before = 1)

  if (!is.null(pooled_out$sense)) {
    pooled_plot_path <- file.path(out_dir, "pooled_sense_plot_happy.png")
    png(pooled_plot_path, width = 800, height = 600)
    plot(pooled_out$sense)
    dev.off()
    message("Saved pooled plot to: ", pooled_plot_path)
  } else {
  message("No sense object for pooled — no plot saved.")
}
} else {
  message("No pooled.csv found at ", pooled_path, " — skipping pooled run.")
}
# Country runs 
country_files <- list.files(country_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(country_files) == 0) {
  warning("No country CSV files found in ", country_dir)
  country_results <- tibble()
} else {
  country_results <- map_dfr(country_files, safe_run_file)
}

# combine and save
final_results <- bind_rows(pooled_res, country_results)

readr::write_csv(final_results, file.path(out_dir, "sensemakr_happy_summary_table.csv"))
print(final_results %>% arrange(unit))

