# Longitudinal alpha-diversity and community-change analysis

# Purpose:
#   1) Calculate Shannon alpha diversity from rarefied ASV count data.
#   2) Plot temporal Shannon diversity and longitudinal community change.
#   3) Test temporal effects using mixed-effects models with sampling event as
#      a categorical predictor.

# Notes:
#   - Shannon diversity is calculated after rarefaction to the minimum
#     sequencing depth across samples.
#   - Longitudinal community change is based on a precomputed Aitchison
#     distance matrix and is not rarefied.
#   - Figures show observed data and LOESS mean curves; model fits are used
#     only for statistical testing.

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lubridate)
  library(vegan)
  library(RColorBrewer)
  library(ragg)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(scales)
  library(grid)
})


# File paths
asv_file  <- file.path("data", "ASV_16s_final.csv")
meta_file <- file.path("data", "metadata_16s_final.xlsx")
dist_file <- file.path("data", "dist_mat.xlsx")


# Output folder
out_dir <- file.path("results", "longitudinal_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# Helper functions
std_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  num <- suppressWarnings(as.integer(gsub("^a0*", "", x)))
  ifelse(is.na(num), NA_character_, paste0("a", num))
}

parse_sampling_date <- function(x) {
  x_date <- suppressWarnings(as.Date(x))
  x_chr  <- as.character(x)
  x_date2 <- suppressWarnings(dmy(x_chr))
  x_date3 <- suppressWarnings(ymd(x_chr))
  coalesce(x_date, x_date2, x_date3)
}

std_site_code <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

add_sampling_event <- function(dat) {
  event_levels <- sort(unique(dat$Sampling.date[!is.na(dat$Sampling.date)]))

  dat %>%
    mutate(
      Sampling.event = factor(
        Sampling.date,
        levels = event_levels,
        labels = paste0("T", seq_along(event_levels))
      )
    )
}

safe_write_csv <- function(x, file, row.names = FALSE, ...) {
  # write.csv() cannot handle list-columns. This wrapper removes them.
  x <- as.data.frame(x)
  list_cols <- names(x)[vapply(x, is.list, logical(1))]

  if (length(list_cols) > 0) {
    message("Removing list-columns before writing CSV: ", paste(list_cols, collapse = ", "))
    x <- x[, !names(x) %in% list_cols, drop = FALSE]
  }

  utils::write.csv(x, file = file, row.names = row.names, ...)
}

make_loess_ci_curve <- function(dat, x_col, y_col, group_col = "Area",
                                span = 0.9, n_grid = 300) {
  # Creates fitted LOESS mean curves and 95% CI boundary lines per group.
  # CI values are drawn as lines rather than ribbons to avoid overlapping fields.
  dat <- dat %>%
    mutate(
      .x_date = as.Date(.data[[x_col]]),
      .x_num  = as.numeric(.x_date),
      .y      = as.numeric(.data[[y_col]]),
      .group  = .data[[group_col]]
    ) %>%
    filter(!is.na(.x_date), is.finite(.x_num), is.finite(.y), !is.na(.group))

  dat %>%
    group_by(.group) %>%
    group_modify(~ {
      d <- .x %>% arrange(.x_num)

      if (nrow(d) < 4 || length(unique(d$.x_num)) < 4) {
        return(tibble())
      }

      x_grid <- seq(min(d$.x_num), max(d$.x_num), length.out = n_grid)
      fit <- loess(
        .y ~ .x_num,
        data = d,
        span = span,
        degree = 2,
        control = loess.control(surface = "direct")
      )

      pred <- predict(fit, newdata = data.frame(.x_num = x_grid), se = TRUE)
      tcrit <- qt(0.975, df = max(1, pred$df))

      tibble(
        Sampling.date = as.Date(x_grid, origin = "1970-01-01"),
        fit = as.numeric(pred$fit),
        ci_lower = as.numeric(pred$fit - tcrit * pred$se.fit),
        ci_upper = as.numeric(pred$fit + tcrit * pred$se.fit)
      )
    }) %>%
    ungroup() %>%
    rename(!!group_col := .group)
}

extract_lrt_p <- function(lrt) {
  lrt_df <- as.data.frame(lrt)

  if ("Pr(>Chisq)" %in% names(lrt_df)) {
    return(lrt_df$`Pr(>Chisq)`[2])
  }

  NA_real_
}

make_lrt_row <- function(lrt, response, test, model_reduced, model_full) {
  lrt_df <- as.data.frame(lrt)

  tibble(
    Response = response,
    Test = test,
    Reduced_model = model_reduced,
    Full_model = model_full,
    df_reduced = lrt_df$Df[1],
    df_full = lrt_df$Df[2],
    AIC_reduced = lrt_df$AIC[1],
    AIC_full = lrt_df$AIC[2],
    BIC_reduced = lrt_df$BIC[1],
    BIC_full = lrt_df$BIC[2],
    logLik_reduced = lrt_df$logLik[1],
    logLik_full = lrt_df$logLik[2],
    Chisq = lrt_df$Chisq[2],
    Chi_df = lrt_df$`Chi Df`[2],
    p.value = extract_lrt_p(lrt)
  )
}

#Plot settings
english_month_labels <- scales::label_date(format = "%b %Y", locale = "en")

area_colors <- brewer.pal(4, "Set1")
names(area_colors) <- c(
  "Cold_room",
  "Processing_finished_product",
  "Raw_material_processing",
  "Receiving_area"
)

area_labels <- c(
  "Cold_room" = "Cold room",
  "Processing_finished_product" = "Processing finished product",
  "Raw_material_processing" = "Raw material processing",
  "Receiving_area" = "Receiving area"
)

# Load and prepare data

# ASV table: rows = ASVs, columns = samples.
asv <- read.csv(asv_file, check.names = FALSE)

asv_ids <- asv[[1]]
asv_mat <- asv[, -1, drop = FALSE]

rownames(asv_mat) <- asv_ids
colnames(asv_mat) <- std_id(colnames(asv_mat))
asv_mat <- as.matrix(asv_mat)
mode(asv_mat) <- "numeric"

# Metadata
meta <- read_excel(meta_file)
meta$SampleID <- std_id(meta[[1]])

meta_df <- meta %>%
  select(-1) %>%
  mutate(
    Sampling.date = parse_sampling_date(Sampling.date),
    Sampling.site.detail = trimws(as.character(Sampling.site.detail)),
    Sampling.Site.Code = std_site_code(Sampling.Site.Code),
    Facility = factor(Facility),
    Area = factor(
      Area,
      levels = c(
        "Cold_room",
        "Processing_finished_product",
        "Raw_material_processing",
        "Receiving_area"
      )
    ),
    Species = factor(Species)
  ) %>%
  column_to_rownames("SampleID")

# Match samples between ASV table and metadata.
common_samples <- intersect(colnames(asv_mat), rownames(meta_df))

if (length(common_samples) == 0) {
  stop("No shared samples found between ASV table and metadata.")
}

asv_mat <- asv_mat[, common_samples, drop = FALSE]
meta_df <- meta_df[common_samples, , drop = FALSE]

#Part A: Alpha diversity — Shannon index from rarefied ASV counts

# vegan::rrarefy() expects samples as rows and ASVs as columns.
asv_counts_for_alpha <- t(asv_mat)

seq_depth_alpha <- rowSums(asv_counts_for_alpha, na.rm = TRUE)
min_depth_alpha <- min(seq_depth_alpha)

cat("\n------------------------------\n")
cat("RAREFACTION FOR SHANNON DIVERSITY\n")
cat("------------------------------\n")
cat("Minimum sequencing depth used for rarefaction:", min_depth_alpha, "reads per sample\n")
cat("Number of samples retained:", sum(seq_depth_alpha >= min_depth_alpha), "of", length(seq_depth_alpha), "\n\n")

set.seed(123)
asv_rarefied_alpha <- vegan::rrarefy(asv_counts_for_alpha, sample = min_depth_alpha)

# Save rarefaction documentation.
seq_depth_summary_alpha <- tibble(
  SampleID = rownames(asv_counts_for_alpha),
  Sequencing_depth = seq_depth_alpha,
  Rarefaction_depth = min_depth_alpha
)

safe_write_csv(
  seq_depth_summary_alpha,
  file.path(out_dir, "Shannon_sequencing_depth_and_rarefaction_depth.csv")
)

safe_write_csv(
  asv_rarefied_alpha %>%
    as.data.frame() %>%
    rownames_to_column("SampleID"),
  file.path(out_dir, "ASV_counts_rarefied_for_Shannon.csv")
)

alpha <- data.frame(
  SampleID = rownames(asv_rarefied_alpha),
  Shannon = vegan::diversity(asv_rarefied_alpha, index = "shannon")
) %>%
  left_join(meta_df %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  filter(
    !is.na(Sampling.date),
    !is.na(Area),
    !is.na(Facility),
    !is.na(Sampling.Site.Code)
  ) %>%
  add_sampling_event() %>%
  mutate(time_num = as.numeric(Sampling.date - min(Sampling.date, na.rm = TRUE)))

alpha_site_all <- alpha %>%
  group_by(Facility, Area, Sampling.Site.Code, Sampling.date, Sampling.event) %>%
  summarise(
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Sampling.date = as.Date(Sampling.date),
    Shannon_mean = as.numeric(Shannon_mean),
    site_id = interaction(Facility, Sampling.Site.Code, drop = TRUE)
  ) %>%
  filter(!is.na(Sampling.date), is.finite(Shannon_mean)) %>%
  arrange(Facility, Area, Sampling.Site.Code, Sampling.date)

alpha_fac_area <- alpha %>%
  group_by(Facility, Area, Sampling.date, Sampling.event) %>%
  summarise(
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(Sampling.date), is.finite(Shannon_mean))

alpha_area_all <- alpha_fac_area %>%
  group_by(Area, Sampling.date, Sampling.event) %>%
  summarise(
    Shannon_mean = mean(Shannon_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Area) %>%
  mutate(n_points = n()) %>%
  ungroup()

alpha_area_curve_data <- alpha_site_all %>%
  group_by(Area) %>%
  mutate(n_points = n()) %>%
  ungroup()

alpha_loess_ci <- make_loess_ci_curve(
  dat = alpha_area_curve_data %>% filter(n_points >= 4),
  x_col = "Sampling.date",
  y_col = "Shannon_mean",
  group_col = "Area",
  span = 0.9,
  n_grid = 300
)

p_alpha <- ggplot() +
  geom_point(
    data = alpha_site_all,
    aes(x = Sampling.date, y = Shannon_mean, color = Area),
    size = 1.9,
    alpha = 0.35
  ) +
  geom_line(
    data = alpha_loess_ci,
    aes(x = Sampling.date, y = ci_lower, color = Area, group = Area),
    linewidth = 0.8,
    alpha = 0.15
  ) +
  geom_line(
    data = alpha_loess_ci,
    aes(x = Sampling.date, y = ci_upper, color = Area, group = Area),
    linewidth = 0.8,
    alpha = 0.15
  ) +
  geom_line(
    data = alpha_loess_ci,
    aes(x = Sampling.date, y = fit, color = Area, group = Area),
    linewidth = 2.4,
    alpha = 1
  ) +
  scale_color_manual(values = area_colors, labels = area_labels) +
  scale_x_date(
    date_breaks = "1 month",
    labels = english_month_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_classic(base_size = 12) +
  theme(
    line = element_line(lineend = "round", linejoin = "round"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12.5),
    legend.key.size = grid::unit(0.65, "cm")
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Mean alpha diversity (Shannon)",
    color = "Area"
  )

print(p_alpha)

ggsave(
  filename = file.path(out_dir, "Figure_alpha_shannon_site_points_mean_curve_95CI_boundary_lines.pdf"),
  plot = p_alpha,
  width = 10,
  height = 5,
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "Figure_alpha_shannon_site_points_mean_curve_95CI_boundary_lines.png"),
  plot = p_alpha,
  width = 10,
  height = 5,
  dpi = 600,
  device = ragg::agg_png
)

#Part B: Longitudinal community change from precomputed Aitchison distances

dist_raw <- read_excel(dist_file)
dist_ids <- std_id(dist_raw[[1]])

dist_mat <- dist_raw[, -1] %>% as.data.frame()
colnames(dist_mat) <- std_id(colnames(dist_mat))

dist_mat <- as.matrix(dist_mat)
mode(dist_mat) <- "numeric"
rownames(dist_mat) <- dist_ids

dist_mat <- dist_mat[!is.na(rownames(dist_mat)), !is.na(colnames(dist_mat)), drop = FALSE]

common_dist_ids <- intersect(rownames(dist_mat), colnames(dist_mat))
dist_mat <- dist_mat[common_dist_ids, common_dist_ids, drop = FALSE]

meta2 <- read_excel(meta_file) %>%
  mutate(
    Sample.ID = std_id(Sample.ID),
    Sampling.date = parse_sampling_date(Sampling.date),
    Sampling.site.detail = trimws(as.character(Sampling.site.detail)),
    Sampling.Site.Code = std_site_code(Sampling.Site.Code),
    Area = factor(
      as.character(Area),
      levels = c(
        "Cold_room",
        "Processing_finished_product",
        "Raw_material_processing",
        "Receiving_area"
      )
    ),
    Facility = factor(Facility)
  ) %>%
  filter(Sample.ID %in% rownames(dist_mat)) %>%
  filter(!is.na(Sampling.date), !is.na(Sampling.Site.Code)) %>%
  add_sampling_event()

event_tbl <- meta2 %>%
  group_by(Sampling.Site.Code, Facility, Area, Sampling.date, Sampling.event) %>%
  summarise(
    sample_ids = list(Sample.ID),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(Sampling.Site.Code, Sampling.date)

calc_prev_distance <- function(current_ids, previous_ids, dist_mat) {
  current_ids <- intersect(current_ids, rownames(dist_mat))
  previous_ids <- intersect(previous_ids, colnames(dist_mat))

  if (length(current_ids) == 0 || length(previous_ids) == 0) {
    return(NA_real_)
  }

  submat <- dist_mat[current_ids, previous_ids, drop = FALSE]
  mean(submat, na.rm = TRUE)
}

event_prev <- event_tbl %>%
  group_by(Facility, Sampling.Site.Code) %>%
  arrange(Sampling.date, .by_group = TRUE) %>%
  mutate(
    prev_sample_ids = lag(sample_ids),
    prev_date = lag(Sampling.date),
    days_since_prev = as.numeric(Sampling.date - prev_date),
    dist_from_previous = purrr::map2_dbl(
      sample_ids,
      prev_sample_ids,
      ~ calc_prev_distance(.x, .y, dist_mat)
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(dist_from_previous)) %>%
  mutate(
    site_id = interaction(Facility, Sampling.Site.Code, drop = TRUE),
    time_num = as.numeric(Sampling.date - min(Sampling.date, na.rm = TRUE))
  )

beta_fac_area <- event_prev %>%
  group_by(Facility, Area, Sampling.date, Sampling.event) %>%
  summarise(
    mean_dist = mean(dist_from_previous, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(Sampling.date), is.finite(mean_dist))

trend_area <- beta_fac_area %>%
  group_by(Area, Sampling.date, Sampling.event) %>%
  summarise(
    mean_dist = mean(mean_dist, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Area) %>%
  mutate(n_points = n()) %>%
  ungroup()

beta_area_curve_data <- event_prev %>%
  group_by(Area) %>%
  mutate(n_points = n()) %>%
  ungroup()

beta_loess_ci <- make_loess_ci_curve(
  dat = beta_area_curve_data %>% filter(n_points >= 4),
  x_col = "Sampling.date",
  y_col = "dist_from_previous",
  group_col = "Area",
  span = 0.9,
  n_grid = 300
)

p_prev <- ggplot() +
  geom_point(
    data = event_prev,
    aes(x = Sampling.date, y = dist_from_previous, color = Area),
    size = 1.9,
    alpha = 0.35
  ) +
  geom_line(
    data = beta_loess_ci,
    aes(x = Sampling.date, y = ci_lower, color = Area, group = Area),
    linewidth = 0.8,
    alpha = 0.15
  ) +
  geom_line(
    data = beta_loess_ci,
    aes(x = Sampling.date, y = ci_upper, color = Area, group = Area),
    linewidth = 0.8,
    alpha = 0.15
  ) +
  geom_line(
    data = beta_loess_ci,
    aes(x = Sampling.date, y = fit, color = Area, group = Area),
    linewidth = 2.2,
    alpha = 1
  ) +
  scale_color_manual(values = area_colors, labels = area_labels) +
  scale_x_date(
    date_breaks = "1 month",
    labels = english_month_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_classic(base_size = 12) +
  theme(
    line = element_line(lineend = "round", linejoin = "round"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12.5),
    legend.key.size = grid::unit(0.65, "cm")
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Distance from previous event",
    color = "Area"
  )

print(p_prev)

ggsave(
  filename = file.path(out_dir, "Figure_distance_from_previous_event_site_points_mean_curve_95CI_boundary_lines.pdf"),
  plot = p_prev,
  width = 10,
  height = 5,
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "Figure_distance_from_previous_event_site_points_mean_curve_95CI_boundary_lines.png"),
  plot = p_prev,
  width = 10,
  height = 5,
  dpi = 600,
  device = ragg::agg_png
)

# Save plotted summary data.
safe_write_csv(alpha_area_all, file.path(out_dir, "Figure_alpha_shannon_area_mean_points.csv"))
safe_write_csv(alpha_area_curve_data, file.path(out_dir, "Figure_alpha_shannon_site_level_curve_data.csv"))
safe_write_csv(alpha_loess_ci, file.path(out_dir, "Figure_alpha_shannon_loess_95CI_boundary_lines.csv"))
safe_write_csv(trend_area, file.path(out_dir, "Figure_distance_from_previous_area_mean_points.csv"))
safe_write_csv(beta_area_curve_data, file.path(out_dir, "Figure_distance_from_previous_site_level_curve_data.csv"))
safe_write_csv(beta_loess_ci, file.path(out_dir, "Figure_distance_from_previous_loess_95CI_boundary_lines.csv"))

# Reduced statistical testing

# Four tests are retained for each response:
#   1) Sampling.event effect
#   2) Area effect adjusted for Sampling.event
#   3) Area x Sampling.event interaction
#   4) Non-linearity test comparing categorical Sampling.event with linear time

cat("\n==============================\n")
cat("PART C: REDUCED STATISTICAL TESTING\n")
cat("==============================\n\n")

# Shannon diversity------------------------------------------------------

alpha_model_data <- alpha %>%
  filter(
    !is.na(Sampling.event),
    !is.na(time_num),
    !is.na(Area),
    !is.na(Facility),
    !is.na(Sampling.Site.Code),
    !is.na(Species),
    !is.na(Shannon)
  ) %>%
  mutate(
    site_id = interaction(Facility, Sampling.Site.Code, drop = TRUE),
    Area = droplevels(Area),
    Sampling.event = droplevels(Sampling.event),
    Facility = droplevels(Facility),
    Species = droplevels(Species)
  )

cat("\n------------------------------\n")
cat("SHANNON DIVERSITY\n")
cat("------------------------------\n\n")

m_alpha_no_time <- lmer(
  Shannon ~ Area + Species + (1 | Facility) + (1 | site_id),
  data = alpha_model_data,
  REML = FALSE
)

m_alpha_event <- lmer(
  Shannon ~ Area + Sampling.event + Species + (1 | Facility) + (1 | site_id),
  data = alpha_model_data,
  REML = FALSE
)

m_alpha_no_area <- lmer(
  Shannon ~ Sampling.event + Species + (1 | Facility) + (1 | site_id),
  data = alpha_model_data,
  REML = FALSE
)

m_alpha_interaction <- lmer(
  Shannon ~ Area * Sampling.event + Species + (1 | Facility) + (1 | site_id),
  data = alpha_model_data,
  REML = FALSE
)

m_alpha_linear <- lmer(
  Shannon ~ Area + time_num + Species + (1 | Facility) + (1 | site_id),
  data = alpha_model_data,
  REML = FALSE
)

lrt_alpha_event <- anova(m_alpha_no_time, m_alpha_event)
cat("1) Sampling.event effect: temporal variation across sampling events\n")
print(lrt_alpha_event)

lrt_alpha_area <- anova(m_alpha_no_area, m_alpha_event)
cat("\n2) Area effect adjusted for Sampling.event\n")
print(lrt_alpha_area)

lrt_alpha_interaction <- anova(m_alpha_event, m_alpha_interaction)
cat("\n3) Area x Sampling.event interaction\n")
print(lrt_alpha_interaction)

lrt_alpha_nonlinear <- anova(m_alpha_linear, m_alpha_event)
cat("\n4) Non-linearity test: categorical Sampling.event vs linear time_num\n")
print(lrt_alpha_nonlinear)

alpha_results_reduced <- bind_rows(
  make_lrt_row(
    lrt_alpha_event,
    response = "Shannon diversity",
    test = "Sampling.event effect",
    model_reduced = "Shannon ~ Area + Species + random effects",
    model_full = "Shannon ~ Area + Sampling.event + Species + random effects"
  ),
  make_lrt_row(
    lrt_alpha_area,
    response = "Shannon diversity",
    test = "Area effect adjusted for Sampling.event",
    model_reduced = "Shannon ~ Sampling.event + Species + random effects",
    model_full = "Shannon ~ Area + Sampling.event + Species + random effects"
  ),
  make_lrt_row(
    lrt_alpha_interaction,
    response = "Shannon diversity",
    test = "Area x Sampling.event interaction",
    model_reduced = "Shannon ~ Area + Sampling.event + Species + random effects",
    model_full = "Shannon ~ Area * Sampling.event + Species + random effects"
  ),
  make_lrt_row(
    lrt_alpha_nonlinear,
    response = "Shannon diversity",
    test = "Non-linearity: Sampling.event factor vs linear time_num",
    model_reduced = "Shannon ~ Area + time_num + Species + random effects",
    model_full = "Shannon ~ Area + Sampling.event + Species + random effects"
  )
) %>%
  mutate(
    q.value = p.adjust(p.value, method = "fdr"),
    singular_reduced_or_full = c(
      isSingular(m_alpha_no_time) | isSingular(m_alpha_event),
      isSingular(m_alpha_no_area) | isSingular(m_alpha_event),
      isSingular(m_alpha_event) | isSingular(m_alpha_interaction),
      isSingular(m_alpha_linear) | isSingular(m_alpha_event)
    )
  )

safe_write_csv(
  alpha_results_reduced,
  file.path(out_dir, "Shannon_reduced_4tests_LRT_FDR.csv")
)

# Distance from previous event

beta_model_data <- event_prev %>%
  filter(
    !is.na(Sampling.event),
    !is.na(time_num),
    !is.na(Area),
    !is.na(Facility),
    !is.na(Sampling.Site.Code),
    !is.na(dist_from_previous)
  ) %>%
  mutate(
    site_id = interaction(Facility, Sampling.Site.Code, drop = TRUE),
    Area = droplevels(Area),
    Sampling.event = droplevels(Sampling.event),
    Facility = droplevels(Facility)
  )

cat("\n------------------------------\n")
cat("DISTANCE FROM PREVIOUS EVENT\n")
cat("------------------------------\n\n")

m_beta_no_time <- lmer(
  dist_from_previous ~ Area + (1 | Facility) + (1 | site_id),
  data = beta_model_data,
  REML = FALSE
)

m_beta_event <- lmer(
  dist_from_previous ~ Area + Sampling.event + (1 | Facility) + (1 | site_id),
  data = beta_model_data,
  REML = FALSE
)

m_beta_no_area <- lmer(
  dist_from_previous ~ Sampling.event + (1 | Facility) + (1 | site_id),
  data = beta_model_data,
  REML = FALSE
)

m_beta_interaction <- lmer(
  dist_from_previous ~ Area * Sampling.event + (1 | Facility) + (1 | site_id),
  data = beta_model_data,
  REML = FALSE
)

m_beta_linear <- lmer(
  dist_from_previous ~ Area + time_num + (1 | Facility) + (1 | site_id),
  data = beta_model_data,
  REML = FALSE
)

lrt_beta_event <- anova(m_beta_no_time, m_beta_event)
cat("1) Sampling.event effect: temporal variation across sampling events\n")
print(lrt_beta_event)

lrt_beta_area <- anova(m_beta_no_area, m_beta_event)
cat("\n2) Area effect adjusted for Sampling.event\n")
print(lrt_beta_area)

lrt_beta_interaction <- anova(m_beta_event, m_beta_interaction)
cat("\n3) Area x Sampling.event interaction\n")
print(lrt_beta_interaction)

lrt_beta_nonlinear <- anova(m_beta_linear, m_beta_event)
cat("\n4) Non-linearity test: categorical Sampling.event vs linear time_num\n")
print(lrt_beta_nonlinear)

beta_results_reduced <- bind_rows(
  make_lrt_row(
    lrt_beta_event,
    response = "Distance from previous event",
    test = "Sampling.event effect",
    model_reduced = "distance ~ Area + random effects",
    model_full = "distance ~ Area + Sampling.event + random effects"
  ),
  make_lrt_row(
    lrt_beta_area,
    response = "Distance from previous event",
    test = "Area effect adjusted for Sampling.event",
    model_reduced = "distance ~ Sampling.event + random effects",
    model_full = "distance ~ Area + Sampling.event + random effects"
  ),
  make_lrt_row(
    lrt_beta_interaction,
    response = "Distance from previous event",
    test = "Area x Sampling.event interaction",
    model_reduced = "distance ~ Area + Sampling.event + random effects",
    model_full = "distance ~ Area * Sampling.event + random effects"
  ),
  make_lrt_row(
    lrt_beta_nonlinear,
    response = "Distance from previous event",
    test = "Non-linearity: Sampling.event factor vs linear time_num",
    model_reduced = "distance ~ Area + time_num + random effects",
    model_full = "distance ~ Area + Sampling.event + random effects"
  )
) %>%
  mutate(
    q.value = p.adjust(p.value, method = "fdr"),
    singular_reduced_or_full = c(
      isSingular(m_beta_no_time) | isSingular(m_beta_event),
      isSingular(m_beta_no_area) | isSingular(m_beta_event),
      isSingular(m_beta_event) | isSingular(m_beta_interaction),
      isSingular(m_beta_linear) | isSingular(m_beta_event)
    )
  )

safe_write_csv(
  beta_results_reduced,
  file.path(out_dir, "Distance_reduced_4tests_LRT_FDR.csv")
)

# Combined overview
combined_reduced_results <- bind_rows(alpha_results_reduced, beta_results_reduced) %>%
  group_by(Response) %>%
  mutate(q.value_within_response = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

safe_write_csv(
  combined_reduced_results,
  file.path(out_dir, "Combined_reduced_4tests_LRT_FDR.csv")
)

cat("\n==============================\n")
cat("Reduced 4-test analysis completed.\n")
cat("Outputs saved to:\n")
cat(out_dir, "\n")
cat("==============================\n")
