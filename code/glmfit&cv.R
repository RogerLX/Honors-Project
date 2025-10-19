suppressPackageStartupMessages({
  library(lubridate)
  library(tidyr)
  library(MASS)
  library(broom)
  library(ggplot2)
  library(readr)
  library(dplyr)
})

if (!exists("df")) {
  data_path <- file.path("patient_clean.csv")
  if (!file.exists(data_path)) {
    stop("`df` must exist in the environment or `patient_clean.csv` must be present in the working directory.")
  }
  message("Reading data from ", data_path)
  df <- read_csv(data_path, show_col_types = FALSE)
}

# ---------------------------------------------------------------------------
# 1) Melbourne metro filter
# ---------------------------------------------------------------------------
df_melb <- df %>%
  filter(!is.na(sceneaddress), !is.na(hospitalname)) %>%
  filter(
    between(long_racf, 144.5, 145.5),
    between(lat_racf,  -38.2, -37.4),
    between(long_hosp, 144.5, 145.5),
    between(lat_hosp,  -38.2, -37.4)
  ) %>%
  mutate(
    from     = toupper(trimws(sceneaddress)),
    to       = toupper(trimws(hospitalname)),
    date     = make_date(year, month, day),
    covid    = if_else(date < as.Date("2020-03-15"), "Pre", "During"),
    day_type = if_else(wday(date) %in% c(1, 7), "Weekend", "Weekday")
  )

# ---------------------------------------------------------------------------
# 2) Observed transfers
# ---------------------------------------------------------------------------
agg_melb <- df_melb %>%
  count(from, to, covid, day_type, name = "count")

print(count(agg_melb, to))

# ---------------------------------------------------------------------------
# 3) Add zeros (expand all combinations)
# ---------------------------------------------------------------------------
all_combos <- expand_grid(
  from     = unique(df_melb$from),
  to       = unique(df_melb$to),
  covid    = c("Pre", "During"),
  day_type = c("Weekday", "Weekend")
)

agg_full <- all_combos %>%
  left_join(agg_melb, by = c("from", "to", "covid", "day_type")) %>%
  mutate(count = if_else(is.na(count), 0L, count))

# ---------------------------------------------------------------------------
# 4) Find top hospitals and RACFs
# ---------------------------------------------------------------------------
q_rmh         <- 0.98  # keep RACFs in the top 2% by transfers -> RMH
q_total       <- 0.98  # keep RACFs in the top 2% by total transfers (all hospitals)
q_top5        <- 0.98  # keep RACFs in the top 2% by transfers to top-5 hospitals
min_top5_hosp <- 2     # must send to at least this many of the top-5 hospitals

rmh_name <- "ROYAL MELBOURNE HOSPITAL - CITY CAMPUS"

# -- Top-5 hospitals by total volume ------------------------------------------------
hospital_counts <- agg_full %>%
  group_by(to) %>%
  summarise(total_transfers = sum(count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_transfers))

top5_hospitals <- head(hospital_counts$to, 5)
hosp_keep <- union(top5_hospitals, rmh_name)

cat("Top 5 hospitals:\n"); print(top5_hospitals)

# -- RACFs with many transfers to RMH -----------------------------------------------
racfs_to_rmh <- agg_full %>%
  filter(to == rmh_name) %>%
  group_by(from) %>%
  summarise(total_to_rmh = sum(count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_to_rmh))

thr_rmh <- quantile(racfs_to_rmh$total_to_rmh, q_rmh, na.rm = TRUE)
racfs_many_to_rmh <- racfs_to_rmh %>%
  filter(total_to_rmh >= thr_rmh) %>%
  pull(from) %>%
  unique()

cat("\nRACFs in top", (1 - q_rmh) * 100, "% for transfers -> RMH:", length(racfs_many_to_rmh), "\n")

# -- RACFs with large total volume ---------------------------------------------------
racf_totals <- agg_full %>%
  group_by(from) %>%
  summarise(total_from = sum(count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_from))

thr_total <- quantile(racf_totals$total_from, q_total, na.rm = TRUE)
racfs_large_total <- racf_totals %>%
  filter(total_from >= thr_total) %>%
  pull(from) %>%
  unique()

cat("RACFs in top", (1 - q_total) * 100, "% by total volume:", length(racfs_large_total), "\n")

# -- RACFs with strong ties to multiple top-5 hospitals ------------------------------
racfs_mult_top5 <- agg_full %>%
  filter(to %in% top5_hospitals) %>%
  group_by(from) %>%
  summarise(
    total_to_top5 = sum(count, na.rm = TRUE),
    n_top5        = n_distinct(to),
    .groups       = "drop"
  ) %>%
  arrange(desc(total_to_top5))

thr_top5 <- quantile(racfs_mult_top5$total_to_top5, q_top5, na.rm = TRUE)
racfs_top5_multi <- racfs_mult_top5 %>%
  filter(total_to_top5 >= thr_top5, n_top5 >= min_top5_hosp) %>%
  pull(from) %>%
  unique()

cat(
  "RACFs with transfers to >=", min_top5_hosp,
  "of top-5 hospitals and above median top5 volume:",
  length(racfs_top5_multi), "\n"
)

# -- Final RACF selection -----------------------------------------------------------
key_racfs <- unique(c(racfs_many_to_rmh, racfs_large_total, racfs_top5_multi))

cat("\nTotal selected RACFs:", length(key_racfs), "\n")
print(head(key_racfs))

# -- Reduced dataset for modelling --------------------------------------------------
sub_df <- agg_full %>%
  filter(from %in% key_racfs, to %in% hosp_keep) %>%
  mutate(
    to       = factor(to),
    from     = factor(from),
    covid    = factor(covid, levels = c("Pre", "During")),
    day_type = factor(day_type, levels = c("Weekday", "Weekend"))
  )

cat(
  "\nReduced dataset rows:", nrow(sub_df),
  "| unique RACFs:", n_distinct(sub_df$from),
  "| hospitals kept:", paste(hosp_keep, collapse = ", "), "\n"
)

# ---------------------------------------------------------------------------
# 5) Fit GLMs
# ---------------------------------------------------------------------------
fit_pois <- glm(
  count ~ covid + day_type + from + to,
  family = poisson(),
  data   = sub_df
)

fit_nb <- glm.nb(
  count ~ covid * from + day_type + to,
  data = sub_df
)

if (!requireNamespace("pscl", quietly = TRUE)) {
  stop("Package `pscl` is required for the zero-inflated Poisson model.")
}

library(pscl)

fit_zip <- zeroinfl(
  count ~ covid + day_type + from + to,
  data = sub_df,
  dist = "poisson"
)

cat("\n--- Poisson GLM summary ---\n")
print(summary(fit_pois))

cat("\n--- Negative binomial GLM summary ---\n")
print(summary(fit_nb))

cat("\n--- Zero-inflated Poisson summary ---\n")
print(summary(fit_zip))

# ---------------------------------------------------------------------------
# 6) Model comparison and diagnostics
# ---------------------------------------------------------------------------
model_stats <- tibble(
  model  = c("Poisson", "Negative binomial", "Zero-inflated Poisson"),
  aic    = c(AIC(fit_pois), AIC(fit_nb), AIC(fit_zip)),
  bic    = c(BIC(fit_pois), BIC(fit_nb), BIC(fit_zip)),
  logLik = c(
    as.numeric(logLik(fit_pois)),
    as.numeric(logLik(fit_nb)),
    as.numeric(logLik(fit_zip))
  ),
  nobs   = c(
    nobs(fit_pois),
    nobs(fit_nb),
    nrow(sub_df)  # zeroinfl lacks an nobs() method; row-count gives the sample size
  )
)

cat("\n--- Model comparison ---\n")
print(model_stats)

augment_model <- function(data, fit, model_name) {
  tibble(
    model         = model_name,
    count         = data$count,
    covid         = data$covid,
    day_type      = data$day_type,
    fitted        = as.numeric(predict(fit, type = "response")),
    resid_pearson = as.numeric(residuals(fit, type = "pearson"))
  )
}

audited_models <- bind_rows(
  augment_model(sub_df, fit_pois, "Poisson"),
  augment_model(sub_df, fit_nb, "Negative binomial"),
  augment_model(sub_df, fit_zip, "Zero-inflated Poisson")
)

plot_prediction_comparison <- function(data, colour_var) {
  ggplot(data, aes(x = count, y = fitted, colour = .data[[colour_var]])) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_x_sqrt() +
    scale_y_sqrt() +
    facet_wrap(~model, scales = "free") +
    labs(
      x = "Observed transfers",
      y = "Predicted transfers",
      colour = colour_var,
      title = paste("Predicted vs observed counts by", colour_var)
    )
}

plot_residuals <- function(data, colour_var) {
  ggplot(data, aes(x = fitted, y = resid_pearson, colour = .data[[colour_var]])) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.6) +
    facet_wrap(~model, scales = "free") +
    labs(
      x = "Predicted transfers",
      y = "Pearson residuals",
      colour = colour_var,
      title = paste("Residual diagnostics by", colour_var)
    )
}

pred_obs_covid_plot  <- plot_prediction_comparison(audited_models, "covid")
pred_obs_day_plot    <- plot_prediction_comparison(audited_models, "day_type")
residuals_covid_plot <- plot_residuals(audited_models, "covid")
residuals_day_plot   <- plot_residuals(audited_models, "day_type")

print(pred_obs_covid_plot)
print(pred_obs_day_plot)
print(residuals_covid_plot)
print(residuals_day_plot)

# ---------------------------------------------------------------------------
# 7) Filter: only pairs with >2 transfers, then NB/Poisson GLM
# ---------------------------------------------------------------------------
sub_gt2 <- sub_df %>%
  filter(count > 2) %>%
  droplevels()

cat(
  "Rows kept:", nrow(sub_gt2), "\n",
  "Unique from:", n_distinct(sub_gt2$from), "\n",
  "Unique to:",   n_distinct(sub_gt2$to),   "\n"
)

fit_nbg2 <- glm.nb(
  count ~ covid + from + day_type + to,
  data = sub_gt2
)

fit_poig2 <- glm(
  count ~ covid + day_type + from + to,
  family = poisson(),
  data   = sub_gt2
)

cat("\n--- Restricted negative binomial GLM summary ---\n")
print(summary(fit_nbg2))

cat("\n--- Restricted Poisson GLM summary ---\n")
print(summary(fit_poig2))

restricted_stats <- tibble(
  model  = c("Restricted negative binomial", "Restricted Poisson"),
  aic    = c(AIC(fit_nbg2), AIC(fit_poig2)),
  bic    = c(BIC(fit_nbg2), BIC(fit_poig2)),
  logLik = c(as.numeric(logLik(fit_nbg2)), as.numeric(logLik(fit_poig2))),
  nobs   = c(nobs(fit_nbg2), nobs(fit_poig2))
)

cat("\n--- Restricted model comparison ---\n")
print(restricted_stats)

restricted_aug <- bind_rows(
  augment_model(sub_gt2, fit_nbg2, "Restricted negative binomial"),
  augment_model(sub_gt2, fit_poig2, "Restricted Poisson")
)

restricted_pred_covid_plot <- plot_prediction_comparison(restricted_aug, "covid")
restricted_pred_day_plot   <- plot_prediction_comparison(restricted_aug, "day_type")
restricted_resid_covid_plot <- plot_residuals(restricted_aug, "covid")
restricted_resid_day_plot   <- plot_residuals(restricted_aug, "day_type")

print(restricted_pred_covid_plot)
print(restricted_pred_day_plot)
print(restricted_resid_covid_plot)
print(restricted_resid_day_plot)

# ---------------------------------------------------------------------------
# 8) Weekly out-of-sample cross-validation (restricted NB model)
# ---------------------------------------------------------------------------
# Build a week-level dataset that keeps the same factor structure as the
# restricted model (count > 2, selected RACFs/hospitals, covid, day_type).
weekly_counts <- df_melb %>%
  filter(from %in% key_racfs, to %in% hosp_keep) %>%
  mutate(
    week     = floor_date(date, unit = "week", week_start = 1L),
    covid    = factor(covid, levels = levels(sub_gt2$covid)),
    day_type = factor(day_type, levels = levels(sub_gt2$day_type)),
    from     = factor(from, levels = levels(sub_gt2$from)),
    to       = factor(to,   levels = levels(sub_gt2$to))
  ) %>%
  count(from, to, week, covid, day_type, name = "count")

weekly_gt2 <- weekly_counts %>%
  semi_join(
    sub_gt2 %>% distinct(from, to, covid, day_type),
    by = c("from", "to", "covid", "day_type")
  ) %>%
  arrange(week, from, to, covid, day_type)

if (nrow(weekly_gt2) == 0) {
  warning("No weekly observations available for the restricted dataset.")
} else {
  perform_weekly_cv <- function(data, formula) {
    weeks <- sort(unique(data$week))
    
    metrics_list <- vector("list", length(weeks))
    preds_list   <- vector("list", length(weeks))
    
    for (i in seq_along(weeks)) {
      wk <- weeks[[i]]
      train <- data %>% filter(week != wk)
      test  <- data %>% filter(week == wk)
      
      if (nrow(test) == 0L || nrow(train) == 0L) {
        metrics_list[[i]] <- tibble(
          week    = wk,
          n_train = nrow(train),
          n_test  = nrow(test),
          rmse    = NA_real_,
          mae     = NA_real_,
          log_score = NA_real_,
          theta     = NA_real_,
          status    = "insufficient data"
        )
        preds_list[[i]] <- tibble()
        next
      }
      
      fit <- tryCatch(
        glm.nb(formula, data = train),
        error = function(e) e
      )
      
      if (inherits(fit, "error")) {
        metrics_list[[i]] <- tibble(
          week    = wk,
          n_train = nrow(train),
          n_test  = nrow(test),
          rmse    = NA_real_,
          mae     = NA_real_,
          log_score = NA_real_,
          theta     = NA_real_,
          status    = paste("fit error:", conditionMessage(fit))
        )
        preds_list[[i]] <- tibble()
        next
      }
      
      link_preds <- tryCatch(
        predict(fit, newdata = test, type = "link", se.fit = TRUE),
        error = function(e) e
      )
      
      if (inherits(link_preds, "error")) {
        metrics_list[[i]] <- tibble(
          week    = wk,
          n_train = nrow(train),
          n_test  = nrow(test),
          rmse    = NA_real_,
          mae     = NA_real_,
          log_score = NA_real_,
          theta     = unname(fit$theta),
          status    = paste("prediction error:", conditionMessage(link_preds))
        )
        preds_list[[i]] <- tibble()
        next
      }
      
      preds <- exp(link_preds$fit)
      se_link <- as.numeric(link_preds$se.fit)
      crit <- qnorm(0.975)
      preds_lower <- exp(link_preds$fit - crit * se_link)
      preds_upper <- exp(link_preds$fit + crit * se_link)
      
      preds <- as.numeric(preds)
      preds_lower <- as.numeric(pmax(preds_lower, 0))
      preds_upper <- as.numeric(pmax(preds_upper, 0))
      residuals <- test$count - preds
      
      metrics_list[[i]] <- tibble(
        week      = wk,
        n_train   = nrow(train),
        n_test    = nrow(test),
        rmse      = sqrt(mean(residuals^2)),
        mae       = mean(abs(residuals)),
        log_score = mean(dnbinom(test$count, mu = preds, size = fit$theta, log = TRUE)),
        theta     = unname(fit$theta),
        status    = "ok"
      )
      
      preds_list[[i]] <- test %>%
        mutate(
          week       = wk,
          predicted  = preds,
          predicted_lower = preds_lower,
          predicted_upper = preds_upper,
          residual   = count - predicted
        )
    }
    
    list(
      metrics      = bind_rows(metrics_list),
      predictions  = bind_rows(preds_list)
    )
  }
  
  cv_out <- perform_weekly_cv(
    weekly_gt2,
    count ~ covid + from + day_type + to
  )
  
  weekly_cv_metrics <- cv_out$metrics %>%
    mutate(week = as.Date(week))
  
  weekly_cv_summary <- weekly_cv_metrics %>%
    filter(status == "ok") %>%
    summarise(
      folds    = n(),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_mae  = mean(mae, na.rm = TRUE),
      mean_log_score = mean(log_score, na.rm = TRUE),
      mean_theta     = mean(theta, na.rm = TRUE)
    )
  
  cat("\n--- Weekly cross-validation metrics (restricted NB) ---\n")
  print(weekly_cv_metrics)
  
  cat("\n--- Weekly cross-validation summary (restricted NB) ---\n")
  print(weekly_cv_summary)
  
  weekly_cv_predictions <- cv_out$predictions %>%
    mutate(week = as.Date(week))
  
  cat("\nSample of out-of-sample predictions:\n")
  print(head(weekly_cv_predictions))
  
  if (nrow(weekly_cv_metrics) > 0) {
    metric_plot_data <- weekly_cv_metrics %>%
      filter(status == "ok") %>%
      dplyr::select(week, rmse, mae, log_score, theta) %>%
      pivot_longer(
        cols = c(rmse, mae, log_score, theta),
        names_to = "metric",
        values_to = "value"
      )
    
    cv_metric_plot <- ggplot(metric_plot_data, aes(x = week, y = value, colour = metric)) +
      geom_line() +
      geom_point() +
      facet_wrap(~metric, scales = "free_y") +
      labs(
        title = "Leave-one-week-out CV diagnostics",
        x = "Held-out week",
        y = "Metric value"
      ) +
      theme(legend.position = "none")
    
    print(cv_metric_plot)
  }
  
  if (nrow(weekly_cv_predictions) > 0) {
    weekly_totals <- weekly_cv_predictions %>%
      group_by(week) %>%
      summarise(
        observed_total  = sum(count, na.rm = TRUE),
        predicted_total = sum(predicted, na.rm = TRUE),
        predicted_lower_total = sum(predicted_lower, na.rm = TRUE),
        predicted_upper_total = sum(predicted_upper, na.rm = TRUE),
        .groups = "drop"
      )
    
    cv_weekly_totals_plot <- ggplot(weekly_totals, aes(x = week)) +
      geom_ribbon(
        aes(x = week, ymin = predicted_lower_total, ymax = predicted_upper_total),
        fill = "#d95f02",
        alpha = 0.2,
        inherit.aes = FALSE
      ) +
      geom_line(aes(y = observed_total, colour = "Observed")) +
      geom_point(aes(y = observed_total, colour = "Observed")) +
      geom_line(aes(y = predicted_total, colour = "Predicted")) +
      geom_point(aes(y = predicted_total, colour = "Predicted")) +
      scale_colour_manual(values = c("Observed" = "#1b9e77", "Predicted" = "#d95f02")) +
      labs(
        title = "Observed vs predicted weekly totals (CV folds)",
        x = "Week",
        y = "Transfers",
        colour = "Series"
      )
    
    print(cv_weekly_totals_plot)
    
    cv_residual_scatter <- ggplot(weekly_cv_predictions, aes(x = predicted, y = residual)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(alpha = 0.6) +
      labs(
        title = "Cross-validated residuals by prediction",
        x = "Predicted transfers",
        y = "Residual (observed - predicted)"
      )
    
    print(cv_residual_scatter)
    
    cv_prediction_interval_plot <- ggplot(weekly_cv_predictions, aes(x = count, y = predicted)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbar(aes(ymin = predicted_lower, ymax = predicted_upper), width = 0.2, alpha = 0.5) +
      geom_point(alpha = 0.7, colour = "#1b9e77") +
      labs(
        title = "Observed vs predicted counts with 95% confidence intervals",
        x = "Observed transfers",
        y = "Predicted transfers"
      )
    
    print(cv_prediction_interval_plot)
  }
}