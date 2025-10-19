# =========================
# 0) Packages
# =========================
# install.packages(c("readr","dplyr","lubridate","ggplot2","geosphere","broom","MASS"))
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(geosphere)  # distHaversine
library(broom)
library(MASS)       # glm.nb

# =========================
# 1) Load + minimal clean
# =========================
df <- read_csv("patient_clean.csv", show_col_types = FALSE)

df_clean <- df %>%
  filter(toupper(trimws(transportflag)) == "Y") %>%
  filter(!is.na(sceneaddress), !is.na(hospitalname)) %>%
  mutate(
    from   = toupper(trimws(sceneaddress)),
    to     = toupper(trimws(hospitalname)),
    date   = make_date(year, month, day),
    period = if_else(date < as.Date("2020-03-15"), "Pre", "During"),
    day_type = if_else(wday(date) %in% c(1,7), "Weekend", "Weekday"),
    distance_km = distHaversine(
      cbind(long_racf, lat_racf),
      cbind(long_hosp, lat_hosp)
    ) / 1000
  ) %>%
  filter(is.finite(distance_km), distance_km > 0)

# =========================
# 2) Collapse to pair-level (edge list)
# =========================
edges <- df_clean %>%
  group_by(from, to) %>%
  summarise(
    n_transfers = n(),                         # count per pair
    distance_km = median(distance_km, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  filter(distance_km > 0)

cat(sprintf("Pairs: %d | median distance = %.1f km\n",
            nrow(edges), median(edges$distance_km)))

# =========================
# 3) Distance vs Count — Plots (log–log)
# =========================

# (A) Scatter with loess; log–log scales
p_scatter <- ggplot(edges, aes(distance_km, n_transfers)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_x_log10() + scale_y_log10() +
  labs(
    x = "Distance between RACF and Hospital (km, log scale)",
    y = "Transfers per pair (log scale)",
    title = "Transfers vs Distance (pair-level, log–log)"
  ) +
  theme_minimal()
print(p_scatter)

# (B) Binned counts by distance using LOG-SPACED bins
# =========================
# Re-bin distances in log scale
# =========================

# Define log-spaced bins for distance
breaks <- exp(seq(log(min(edges$distance_km, na.rm = TRUE)),
                  log(max(edges$distance_km, na.rm = TRUE)),
                  length.out = 25))  # ~25 bins

edges_bins <- edges %>%
  mutate(bin = cut(distance_km, breaks = breaks, include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    total_transfers = sum(n_transfers),
    n_pairs = n(),
    bin_min = min(distance_km),
    bin_max = max(distance_km),
    .groups = "drop"
  ) %>%
  mutate(bin_mid_km = sqrt(bin_min * bin_max))   # geometric mean as midpoint

# =========================
# Plot transfer volume by log-distance bins
# =========================
p_bins <- ggplot(edges_bins, aes(x = bin_mid_km, y = total_transfers)) +
  geom_col(width = 0.15 * edges_bins$bin_mid_km) +  # width ~15% of mid-point
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Distance (km, log scale)",
    y = "Total transfers in bin (log scale)",
    title = "Transfer volume by log-distance bins (log-log)"
  ) +
  theme_minimal()

print(p_bins)

# =========================
# 4) Correlation diagnostics
# =========================

# Spearman (monotone)
cor_spear <- cor.test(edges$distance_km, edges$n_transfers, method = "spearman")
print(cor_spear)

# Pearson on logs (power-law)
log_df <- edges %>%
  mutate(
    log_dist = log(distance_km),
    log_n    = log(pmax(n_transfers, 1))  # avoid log(0)
  )
cor_pearson_log <- cor.test(log_df$log_dist, log_df$log_n, method = "pearson")
print(cor_pearson_log)

# =========================
# 5) Regression models
# =========================

# 5a) Log–log OLS: log(n) ~ log(distance)   => n ~ distance^beta
fit_loglog <- lm(log_n ~ log_dist, data = log_df)
cat("\n--- Log–log OLS (n ~ distance^beta) ---\n")
print(glance(fit_loglog))
print(tidy(fit_loglog))

# 5b) Quasi-Poisson on counts
fit_qpois <- glm(n_transfers ~ log(distance_km), data = edges, family = quasipoisson())
cat("\n--- Quasi-Poisson: n ~ log(distance) ---\n")
print(summary(fit_qpois))

# 5c) Negative Binomial on counts
fit_nb <- glm.nb(n_transfers ~ log(distance_km), data = edges)
cat("\n--- Negative Binomial: n ~ log(distance) ---\n")
print(summary(fit_nb))

# Predicted curve overlay (NegBin) on log–log axes
pred_grid <- data.frame(distance_km = exp(seq(log(min(edges$distance_km)),
                                              log(max(edges$distance_km)),
                                              length.out = 300)))
pred_grid$fit_nb <- predict(fit_nb, newdata = pred_grid, type = "response")

p_nb <- ggplot(edges, aes(distance_km, n_transfers)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_line(data = pred_grid, aes(distance_km, fit_nb), linewidth = 1.1) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Distance (km, log scale)", y = "Transfers per pair (log scale)",
       title = "Negative Binomial fit (log–log): Transfers ~ Distance") +
  theme_minimal()
print(p_nb)

# =========================
# 6) Pre vs During COVID (log–log plots)
# =========================
edges_period <- df_clean %>%
  group_by(period, from, to) %>%
  summarise(n_transfers = n(),
            distance_km = median(distance_km, na.rm=TRUE),
            .groups = "drop") %>%
  filter(distance_km > 0)

p_period <- ggplot(edges_period, aes(distance_km, n_transfers, color = period)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Distance (km, log scale)", y = "Transfers per pair (log scale)",
       title = "Transfers vs Distance by COVID period (log–log)") +
  theme_minimal()
print(p_period)

# Interaction on log scale
edges_period_log <- edges_period %>%
  mutate(log_dist = log(distance_km), log_n = log(pmax(n_transfers,1)))
fit_interact <- lm(log_n ~ log_dist * period, data = edges_period_log)
cat("\n--- Log–log OLS with interaction (period) ---\n")
print(anova(fit_interact))
print(tidy(fit_interact))

# =========================
# 7) Weekday vs Weekend (log–log plots)
# =========================
edges_day <- df_clean %>%
  group_by(day_type, from, to) %>%
  summarise(n_transfers = n(),
            distance_km = median(distance_km, na.rm=TRUE),
            .groups = "drop") %>%
  filter(distance_km > 0)

p_day <- ggplot(edges_day, aes(distance_km, n_transfers, color = day_type)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Distance (km, log scale)", y = "Transfers per pair (log scale)",
       title = "Transfers vs Distance: Weekday vs Weekend (log–log)") +
  theme_minimal()
print(p_day)

edges_day_log <- edges_day %>%
  mutate(log_dist = log(distance_km), log_n = log(pmax(n_transfers,1)))
fit_day_interact <- lm(log_n ~ log_dist * day_type, data = edges_day_log)
cat("\n--- Log–log OLS with interaction (day type) ---\n")
print(anova(fit_day_interact))
print(tidy(fit_day_interact))
