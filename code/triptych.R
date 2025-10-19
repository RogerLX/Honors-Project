# =========================
# Libraries
# =========================
library(dplyr)
library(lubridate)
library(Matrix)
library(readr)
library(superheat)

# =========================
# Load & clean
# =========================
df <- read_csv("patient_clean.csv")

df_clean <- df %>%
  filter(toupper(trimws(transportflag)) == "Y") %>%
  filter(!is.na(sceneaddress), !is.na(hospitalname)) %>%
  mutate(
    date         = make_date(year, month, day),
    day_type     = if_else(wday(date) %in% c(1, 7), "Weekend", "Weekday"),
    covid_period = if_else(date < as.Date("2020-03-15"), "Pre", "During"),
    from         = toupper(trimws(sceneaddress)),
    to           = toupper(trimws(hospitalname))
  )

row_names <- sort(unique(df_clean$from))  # RACFs
col_names <- sort(unique(df_clean$to))    # Hospitals

# =========================
# Helpers
# =========================
get_days <- function(dat, filter_expr) {
  dat %>% filter(!!filter_expr) %>% pull(date) %>% n_distinct()
}

create_adj_matrix <- function(data, n_days, row_names, col_names) {
  data %>%
    count(from, to, name = "weight") %>%
    mutate(
      from = factor(from, levels = row_names),
      to   = factor(to,   levels = col_names),
      norm_weight = sqrt(weight / n_days)  # variance-stabilised per-day rate
    ) %>%
    {
      sparseMatrix(
        i = as.integer(.$from),
        j = as.integer(.$to),
        x = .$norm_weight,
        dims = c(length(row_names), length(col_names)),
        dimnames = list(row_names, col_names)
      )
    }
}

top_abs_diff_matrix <- function(M, top_n = 200) {
  dense <- as.matrix(M)
  long  <- as.data.frame(as.table(dense)) %>%
    rename(RACF = Var1, Hospital = Var2, diff = Freq) %>%
    mutate(absdiff = abs(diff)) %>%
    arrange(desc(absdiff)) %>%
    slice(1:top_n)
  
  rows_keep <- unique(long$RACF)
  cols_keep <- unique(long$Hospital)
  dense[rows_keep, cols_keep, drop = FALSE]
}

# A robust superheat plot wrapper (no abbreviations; all labels shown)
plot_sh <- function(mat, title, outfile, col_limits = NULL,
                    width = 2600, height = 2200, res = 220) {
  
  # symmetric colour limits if not provided
  if (is.null(col_limits)) {
    lim <- max(abs(range(mat, na.rm = TRUE)))
    col_limits <- c(-lim, lim)
  }
  
  png(outfile, width = width, height = height, res = res)
  
  superheat(
    mat,
    title = title,
    title.alignment = "center",
    
    # Colour palette & clipping
    heat.pal = colorRampPalette(c("blue", "white", "red"))(201),
    # superheat clips to range of the matrix; to emphasise symmetry,
    # we rescale values beforehand:
    # (do not modify labels—keep full names)
    legend = TRUE,
    
    # Axis labels — SHOW FULL NAMES
    left.label.text.size = 4.0,            # y-axis (RACFs)
    left.label.text.alignment = "right",
    bottom.label.text.size = 3.3,          # x-axis (Hospitals)
    bottom.label.text.angle = 90,          # vertical to avoid overlap
    
    # Keep current ordering (no clustering) since we pre-sort
    pretty.order.rows = FALSE,
    pretty.order.cols = FALSE,
    
    # Thin white grid for readability
    grid.hline.col = "white",
    grid.vline.col = "white"
  )
  
  dev.off()
}

# =========================
# A) Weekday – Weekend
# =========================
n_days_weekday <- get_days(df_clean, quote(day_type == "Weekday"))
n_days_weekend <- get_days(df_clean, quote(day_type == "Weekend"))

A_weekday <- create_adj_matrix(filter(df_clean, day_type == "Weekday"),
                               n_days_weekday, row_names, col_names)
A_weekend <- create_adj_matrix(filter(df_clean, day_type == "Weekend"),
                               n_days_weekend, row_names, col_names)

# Order for this comparison using combined activity
activity_day <- rowSums(A_weekday) + colSums(A_weekday) +
  rowSums(A_weekend) + colSums(A_weekend)
row_order_day <- row_names[order(activity_day[row_names], decreasing = TRUE)]
col_order_day <- col_names[order(activity_day[col_names], decreasing = TRUE)]

A_weekday <- A_weekday[row_order_day, col_order_day]
A_weekend <- A_weekend[row_order_day, col_order_day]

D_daytype <- A_weekday - A_weekend  # Weekday – Weekend

# Top-200 absolute changes for plotting
D_daytype_top <- top_abs_diff_matrix(D_daytype, top_n = 200)

# Symmetric colour limits for this comparison
lim_day <- max(abs(range(D_daytype_top, na.rm = TRUE)))

plot_sh(
  D_daytype_top,
  "Top 200 Transfer Differences: Weekday – Weekend (normalized rate)",
  "sh_weekday_minus_weekend.png",
  col_limits = c(-lim_day, lim_day)
)

# =========================
# B) During-COVID – Pre-COVID
# =========================
n_days_pre    <- get_days(df_clean, quote(covid_period == "Pre"))
n_days_during <- get_days(df_clean, quote(covid_period == "During"))

A_pre    <- create_adj_matrix(filter(df_clean, covid_period == "Pre"),
                              n_days_pre, row_names, col_names)
A_during <- create_adj_matrix(filter(df_clean, covid_period == "During"),
                              n_days_during, row_names, col_names)

# Order for this comparison using combined activity
activity_cov <- rowSums(A_pre) + colSums(A_pre) +
  rowSums(A_during) + colSums(A_during)
row_order_cov <- row_names[order(activity_cov[row_names], decreasing = TRUE)]
col_order_cov <- col_names[order(activity_cov[col_names], decreasing = TRUE)]

A_pre    <- A_pre[row_order_cov, col_order_cov]
A_during <- A_during[row_order_cov, col_order_cov]

D_covid <- A_during - A_pre  # During – Pre

# Top-200 absolute changes for plotting
D_covid_top <- top_abs_diff_matrix(D_covid, top_n = 200)

# Symmetric colour limits for this comparison
lim_covid <- max(abs(range(D_covid_top, na.rm = TRUE)))

plot_sh(
  D_covid_top,
  "Top 200 Transfer Differences: During-COVID – Pre-COVID (normalized rate)",
  "sh_during_minus_pre.png",
  col_limits = c(-lim_covid, lim_covid)
)

# ============
# Done.
# Outputs:
#   sh_weekday_minus_weekend.png
#   sh_during_minus_pre.png
# ============