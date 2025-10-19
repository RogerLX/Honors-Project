# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Load the data
df <- read_csv("patient_clean.csv")

# Clean and prepare data
df_clean <- df %>%
  mutate(
    age = as.numeric(age),
    diagnosis = as.factor(diagnosis),
    dispatch = as.factor(dispatch),
    casenature = as.factor(casenature)
  ) %>%
  drop_na(age)

# Histogram of Age
ggplot(df_clean, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Age", x = "Age", y = "Count") +
  theme_minimal()

# Top 10 Diagnoses
df_clean %>%
  count(diagnosis, sort = TRUE) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = reorder(diagnosis, n), y = n)) +
  geom_bar(stat = "identity", fill = "salmon", color = "black") +
  coord_flip() +
  labs(title = "Top 10 Diagnoses", x = "Diagnosis", y = "Count") +
  theme_minimal()

# Top 10 Dispatch Reasons
df_clean %>%
  count(dispatch, sort = TRUE) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = reorder(dispatch, n), y = n)) +
  geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
  coord_flip() +
  labs(title = "Top 10 Dispatch Reasons", x = "Dispatch Reason", y = "Count") +
  theme_minimal()

# Top 10 Case Natures
df_clean %>%
  count(casenature, sort = TRUE) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = reorder(casenature, n), y = n)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "black") +
  coord_flip() +
  labs(title = "Top 10 Case Natures", x = "Case Nature", y = "Count") +
  theme_minimal()

# Top RACF → Hospital Pairs by Transfer Count
top_pairs <- df_edges %>%
  arrange(desc(weight)) %>%
  select(from, to, weight) %>%
  head(20)

print(top_pairs)

#Aggregate by Facility Name Only (Ignore lat/lon differences)
##Sometimes transfers are split across minor address variations. To fix that, group only by names:
##This gives you a cleaner transfer matrix between facilities

df_pairs_agg <- df %>%
  filter(toupper(trimws(transportflag)) == "Y") %>%
  mutate(
    from = toupper(trimws(sceneaddress)),
    to = toupper(trimws(hospitalname))
  ) %>%
  group_by(from, to) %>%
  summarise(weight = n(), .groups = "drop") %>%
  arrange(desc(weight))

head(df_pairs_agg, 20)




library(dplyr)
library(igraph)





# Install if not done
install.packages(c("sf", "ggplot2", "ggspatial", "rnaturalearth", "rnaturalearthdata", "dplyr"))

# Install required packages if missing
install.packages(c("sf", "ggplot2", "rnaturalearth", "rnaturalearthdata", "dplyr"))

library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# 1. Filter and prepare transfer data
df_edges <- df %>%
  filter(toupper(trimws(transportflag)) == "Y") %>%
  filter(!is.na(hospitalname), !is.na(sceneaddress),
         !is.na(long_hosp), !is.na(lat_hosp),
         !is.na(long_racf), !is.na(lat_racf)) %>%
  mutate(
    from = toupper(trimws(sceneaddress)),
    to = toupper(trimws(hospitalname)),
    long_from = long_racf,
    lat_from = lat_racf,
    long_to = long_hosp,
    lat_to = lat_hosp
  ) %>%
  group_by(from, to, long_from, lat_from, long_to, lat_to) %>%
  summarise(weight = n(), .groups = "drop")

# 2. Create edges as spatial lines
edges_sf <- df_edges %>%
  rowwise() %>%
  mutate(geometry = st_sfc(st_linestring(matrix(c(long_from, long_to, lat_from, lat_to), ncol = 2, byrow = FALSE)), crs = 4326)) %>%
  st_as_sf()

# 3. Create nodes for RACFs and Hospitals
racf_nodes <- df_edges %>%
  select(name = from, long = long_from, lat = lat_from) %>%
  distinct() %>%
  mutate(type = "RACF")

hospital_nodes <- df_edges %>%
  select(name = to, long = long_to, lat = lat_to) %>%
  distinct() %>%
  mutate(type = "Hospital")

nodes_sf <- bind_rows(racf_nodes, hospital_nodes) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

# 4. Load Victoria map
australia <- ne_states(country = "australia", returnclass = "sf")
vic <- australia %>% filter(name == "Victoria")

# 5. Plot the network over Victoria
ggplot() +
  geom_sf(data = vic, fill = "white", color = "black") +
  geom_sf(data = edges_sf, aes(size = weight), color = "#94B447", alpha = 0.3, show.legend = "line") +
  guides(size = guide_legend(override.aes = list(color = "#94B447", alpha = 3))) +
  geom_sf(data = nodes_sf, aes(color = type), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("RACF" = "orange", "Hospital" = "steelblue")) +
  scale_size_continuous(range = c(0.2, 2)) +
  theme_minimal() +
  labs(
    title = "Transfers from RACFs to Hospitals in Victoria",
    color = "Node Type", size = "Transfer Count"
  )
# density+strength+


# Time Varying Comparison 

#1. Add Date and COVID Period Label
library(dplyr)
library(lubridate)

library(dplyr)
library(lubridate)

df <- df %>%
  mutate(
    date = make_date(year, month, day),
    covid_period = case_when(
      date < as.Date("2020-03-15") ~ "Pre-COVID",
      TRUE ~ "During-COVID"
    )
  )

#2. Create a reusable function to return an igraph object and its adjacency matrix
library(igraph)

create_network_sf <- function(df_subset) {
  df_edges <- df_subset %>%
    filter(toupper(trimws(transportflag)) == "Y") %>%
    filter(!is.na(hospitalname), !is.na(sceneaddress),
           !is.na(long_hosp), !is.na(lat_hosp),
           !is.na(long_racf), !is.na(lat_racf)) %>%
    mutate(
      from = toupper(trimws(sceneaddress)),
      to = toupper(trimws(hospitalname)),
      long_from = long_racf,
      lat_from = lat_racf,
      long_to = long_hosp,
      lat_to = lat_hosp
    ) %>%
    group_by(from, to, long_from, lat_from, long_to, lat_to) %>%
    summarise(weight = n(), .groups = "drop")
  
  # Create igraph object
  g <- graph_from_data_frame(df_edges %>% select(from, to, weight), directed = TRUE)
  
  # Create adjacency matrix (non-sparse for easy comparison)
  adj_matrix <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  
  # Create spatial features
  edges_sf <- df_edges %>%
    rowwise() %>%
    mutate(geometry = st_sfc(st_linestring(
      matrix(c(long_from, long_to, lat_from, lat_to), ncol = 2, byrow = FALSE)
    ), crs = 4326)) %>%
    st_as_sf()
  
  racf_nodes <- df_edges %>%
    select(name = from, long = long_from, lat = lat_from) %>%
    distinct() %>%
    mutate(type = "RACF")
  
  hospital_nodes <- df_edges %>%
    select(name = to, long = long_to, lat = lat_to) %>%
    distinct() %>%
    mutate(type = "Hospital")
  
  nodes_sf <- bind_rows(racf_nodes, hospital_nodes) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326)
  
  list(edges = edges_sf, nodes = nodes_sf, graph = g, adj = adj_matrix)
}


#3.Create 4 subsets using shared node list

net_weekdays <- create_network_sf(df %>% filter(daytype == "weekdays"))
net_weekends <- create_network_sf(df %>% filter(daytype == "weekends"))
net_precovid <- create_network_sf(df %>% filter(covid_period == "Pre-COVID"))
net_duringcovid <- create_network_sf(df %>% filter(covid_period == "During-COVID"))

#4.Base map and shared layout (bounding box)
library(rnaturalearth)
vic <- ne_states(country = "australia", returnclass = "sf") %>%
  filter(name == "Victoria")

bbox <- st_bbox(bind_rows(
  net_weekdays$nodes,
  net_weekends$nodes,
  net_precovid$nodes,
  net_duringcovid$nodes
))

#5.Define plotting function
library(ggplot2)

plot_network_map <- function(vic, edges_sf, nodes_sf, title, bbox) {
  ggplot() +
    geom_sf(data = vic, fill = "white", color = "black") +
    geom_sf(data = edges_sf, aes(size = weight), color = "#94B447", alpha = 0.3, show.legend = "line") +
    guides(size = guide_legend(override.aes = list(color = "#94B447", alpha = 1))) +
    geom_sf(data = nodes_sf, aes(color = type), size = 0.6, alpha = 0.5) +
    scale_color_manual(values = c("RACF" = "orange", "Hospital" = "steelblue")) +
    scale_size_continuous(range = c(0.2, 2)) +
    coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax), expand = FALSE) +
    theme_minimal() +
    labs(title = title, color = "Node Type", size = "Transfer Count")
}
#6.Generate the four maps
p1 <- plot_network_map(vic, net_weekdays$edges, net_weekdays$nodes, "Weekday Transfers", bbox)
p2 <- plot_network_map(vic, net_weekends$edges, net_weekends$nodes, "Weekend Transfers", bbox)
p3 <- plot_network_map(vic, net_precovid$edges, net_precovid$nodes, "Pre-COVID Transfers", bbox)
p4 <- plot_network_map(vic, net_duringcovid$edges, net_duringcovid$nodes, "During-COVID Transfers", bbox)


library(patchwork)

# Plot 1: Pre-COVID vs During-COVID
p_covid <- plot_network_map(vic, net_precovid$edges, net_precovid$nodes, "Pre-COVID Transfers", bbox) |
  plot_network_map(vic, net_duringcovid$edges, net_duringcovid$nodes, "During-COVID Transfers", bbox)

# Plot 2: Weekdays vs Weekends
p_week <- plot_network_map(vic, net_weekdays$edges, net_weekdays$nodes, "Weekday Transfers", bbox) |
  plot_network_map(vic, net_weekends$edges, net_weekends$nodes, "Weekend Transfers", bbox)

ggsave("comparison_covid.png", plot = p_covid, width = 12, height = 6, dpi = 300)
ggsave("comparison_week.png", plot = p_week, width = 12, height = 6, dpi = 300)

# Ensure both are same size or pad if needed

## Adjacency Matrices Comparison
# Dimensions (number of unique facilities)
dim(net_weekdays$adj)
dim(net_weekends$adj)
dim(net_precovid$adj)
dim(net_duringcovid$adj)

# Sparsity
mean(net_weekdays$adj > 0)
mean(net_weekends$adj > 0)
mean(net_precovid$adj > 0)
mean(net_duringcovid$adj > 0)

# Total transfers
sum(net_weekdays$adj)
sum(net_weekends$adj)
sum(net_precovid$adj)
sum(net_duringcovid$adj)

# Plot
library(tidyr)
library(dplyr)
library(Matrix)
library(ggplot2)

plot_diff_heatmap <- function(edges1, edges2, label1, label2, title_suffix = "", top_n = NULL) {
  # Combine all unique node names
  all_nodes <- union(
    unique(c(edges1$from, edges1$to)),
    unique(c(edges2$from, edges2$to))
  )
  
  # Create adjacency matrices
  adj1 <- create_adjacency_matrix(edges1, all_nodes)
  adj2 <- create_adjacency_matrix(edges2, all_nodes)
  
  # Compute difference
  adj_diff <-adj1
  
  # Convert to long format
  adj_diff_long <- as.data.frame(as.table(as.matrix(adj_diff)))
  colnames(adj_diff_long) <- c("from", "to", "difference")
  
  # Filter non-zero differences
  adj_diff_long <- adj_diff_long %>%
    filter(difference != 0)
  
  # If top_n is specified, apply filter
  if (!is.null(top_n)) {
    adj_diff_long <- adj_diff_long %>%
      mutate(abs_diff = abs(difference)) %>%
      arrange(desc(abs_diff)) %>%
      slice(1:top_n)
  }
  
  # Plot
  heatmap_plot <- ggplot(adj_diff_long, aes(x = from, y = to, fill = difference)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal(base_size = 10) +
    labs(
      title = paste("Differences in Transfers:", label2, "-", label1, title_suffix),
      x = "From (RACF)", y = "To (Hospital)", fill = "Δ Transfers"
    ) +
    theme(
      axis.text = element_blank()
    )
  
  # Return plot + matrices
  return(list(
    plot = heatmap_plot,
    adj_matrix_1 = adj1,
    adj_matrix_2 = adj2,
    adj_matrix_diff = adj_diff
  ))
}


#Result
result_day <- plot_diff_heatmap(net_weekdays$edges, net_weekends$edges, "Weekdays", "Weekends")
result_day$plot

result_cov <- plot_diff_heatmap(net_precovid$edges, net_duringcovid$edges, "Pre-COVID", "During-COVID")
result_cov$plot
result_cov$adj_matrix_1

#Adjust
library(dplyr)
library(ggplot2)

# Limit to top 200 absolute differences
adj_diff_top <- adj_diff_long %>%
  arrange(desc(abs(difference))) %>%
  slice(1:200)

# Plot
ggplot(adj_diff_top, aes(x = from, y = to, fill = difference)) +
  geom_tile(color = "grey90") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white", midpoint = 0,
    name = "Δ Transfers"
  ) +
  labs(
    title = "Top 200 Differences in Transfers (During-COVID - Pre-COVID)",
    x = "From (RACF)", y = "To (Hospital)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "right"
  )






# 1. Create difference matrices
adj_diff_covid <- net_duringcovid$adj - net_precovid$adj
adj_diff_daytype <- net_weekdays$adj - net_weekends$adj

# 2. Define a plotting function
plot_diff_matrix <- function(mat, title, file_name) {
  df_melt <- melt(mat)
  
  p <- ggplot(df_melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = title,
      x = "From (RACF)", y = "To (Hospital)",
      fill = "Δ Transfers"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  ggsave(file_name, plot = p, width = 8, height = 6, dpi = 300)
  return(p)
}

# 3. Generate and save both plots
p_covid <- plot_diff_matrix(adj_diff_covid, "Transfer Change: During-COVID vs Pre-COVID", "adj_diff_covid.png")
p_daytype <- plot_diff_matrix(adj_diff_daytype, "Transfer Change: Weekdays vs Weekends", "adj_diff_daytype.png")

##SBM
#1: Prepare adjacency matrices
library(igraph)

# Convert to graph objects
g_pre <- graph_from_data_frame(net_precovid$edges %>% select(from, to, weight), directed = TRUE)
g_covid <- graph_from_data_frame(net_duringcovid$edges %>% select(from, to, weight), directed = TRUE)

# Adjacency matrices
A_pre <- as_adjacency_matrix(g_pre, attr = "weight", sparse = TRUE)
A_covid <- as_adjacency_matrix(g_covid, attr = "weight", sparse = TRUE)

#2: Combine and embed
library(igraph)
# Create a unified set of nodes
# Step 1: Union of all node names
all_nodes <- union(V(g_pre)$name, V(g_covid)$name)

# Step 2: Reorder graphs to same node set
g_pre_full <- simplify(add_vertices(g_pre, nv = length(setdiff(all_nodes, V(g_pre)$name)),
                                    name = setdiff(all_nodes, V(g_pre)$name)))
g_covid_full <- simplify(add_vertices(g_covid, nv = length(setdiff(all_nodes, V(g_covid)$name)),
                                      name = setdiff(all_nodes, V(g_covid)$name)))

# Step 3: Sort vertices to match
g_pre_sorted <- permute(g_pre_full, match(all_nodes, V(g_pre_full)$name))
g_covid_sorted <- permute(g_covid_full, match(all_nodes, V(g_covid_full)$name))

# Step 4: Get adjacency matrices
A_pre <- as_adjacency_matrix(g_pre_sorted, attr = "weight", sparse = TRUE)
A_covid <- as_adjacency_matrix(g_covid_sorted, attr = "weight", sparse = TRUE)

# Now they have same dimensions
dim(A_pre)  # e.g., 1200 × 1200
dim(A_covid)  # also 1200 × 1200

# Step 5: Now you can concatenate
A_combined <- cbind(A_pre, A_covid)


# Combine side-by-side (if node IDs are aligned; otherwise use union and reorder)
A_combined <- cbind(as.matrix(A_pre), as.matrix(A_covid))

# Perform SVD embedding
#a real numeric sparse matrix
library(Matrix)

# Ensure both input matrices are sparse real
A_pre <- as(A_pre, "dgCMatrix")
A_covid <- as(A_covid, "dgCMatrix")

# Combine them
A_combined <- cbind(A_pre, A_covid)

# Confirm class
class(A_combined)
# Should return "dgCMatrix" (sparse, numeric)

library(sparsesvd)
svds <- sparsesvd(A_combined, k = 2)  # try d=5

# Embedding
yhat <- svds$v %*% diag(sqrt(svds$d))

#3: Fit GMM and test
library(mclust)

n1 <- ncol(A_pre)
n2 <- ncol(A_covid)

gmm_pre <- Mclust(yhat[1:n1, ])
gmm_covid <- Mclust(yhat[(n1 + 1):(n1 + n2), ])

# Compare
teststat <- -2 * (gmm_pre$loglik - gmm_covid$loglik)
df <- abs(gmm_pre$df - gmm_covid$df)
p_value <- pchisq(teststat, df = df, lower.tail = FALSE)

print(p_value)
#4：Visualize

# 1. Project embeddings using the SVD output
# Assume svds <- sparsesvd(A_combined, k = 5) already done
# Take the right singular vectors (V) and scale by sqrt of singular values
embedding <- svds$v %*% diag(sqrt(svds$d))

# 2. Add labels for which half is pre-COVID vs during-COVID
embedding_df <- as_tibble(embedding) %>%
  mutate(group = factor(rep(c("Pre-COVID", "During-COVID"),
                            each = ncol(A_combined) / 2)))

# 3. Plot in 2D space
ggplot(embedding_df, aes(x = V1, y = V2, color = group)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Spectral Embedding of RACF-Hospital Transfer Network",
    x = "Embedding Dimension 1",
    y = "Embedding Dimension 2",
    color = "Time Period"
  ) +
  theme_minimal()




#network structure metrics
library(igraph)

analyze_network <- function(edges_sf, nodes_sf) {
  # Convert to igraph
  edges_df <- st_drop_geometry(edges_sf)
  g <- graph_from_data_frame(d = edges_df %>% select(from, to, weight),
                             vertices = st_drop_geometry(nodes_sf), directed = TRUE)
  
  # Centrality measures
  deg_in <- degree(g, mode = "in")
  deg_out <- degree(g, mode = "out")
  btw <- betweenness(g, directed = TRUE)
  cls <- closeness(g, mode = "out")
  
  # Community detection
  g_undirected <- as.undirected(g, mode = "collapse")
  comm <- cluster_louvain(g_undirected)
  modularity_score <- modularity(comm)
  
  # Assign node metrics
  V(g)$in_degree <- deg_in
  V(g)$out_degree <- deg_out
  V(g)$betweenness <- btw
  V(g)$closeness <- cls
  V(g)$community <- membership(comm)
  
  # Create summary dataframe
  node_stats <- data.frame(
    name = V(g)$name,
    type = V(g)$type,
    in_degree = deg_in,
    out_degree = deg_out,
    betweenness = btw,
    closeness = cls,
    community = membership(comm)
  )
  
  list(stats = node_stats, modularity = modularity_score, communities = comm)
}
# Weekdays
result_weekdays <- analyze_network(net_weekdays$edges, net_weekdays$nodes)

# Weekends
result_weekends <- analyze_network(net_weekends$edges, net_weekends$nodes)

# Pre-COVID
result_precovid <- analyze_network(net_precovid$edges, net_precovid$nodes)

# During-COVID
result_duringcovid <- analyze_network(net_duringcovid$edges, net_duringcovid$nodes)

result_precovid$modularity
result_duringcovid$modularity
result_weekdays$modularity
result_weekends$modularity

head(result_precovid$stats[order(-result_precovid$stats$betweenness), ], 10)

write.csv(result_precovid$stats, "precovid_node_stats.csv", row.names = FALSE)
write.csv(result_duringcovid$stats, "duringcovid_node_stats.csv", row.names = FALSE)
write.csv(result_weekdays$stats, "weekdays_node_stats.csv", row.names = FALSE)
write.csv(result_weekends$stats, "weekends_node_stats.csv", row.names = FALSE)

# Load required libraries
library(tidyverse)
library(sf)

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Load transfer data
df <- read_csv("patient_clean.csv")

# Step 1: Filter valid transfers
transfers_all <- df %>%
  filter(toupper(trimws(transportflag)) == "Y") %>%
  filter(!is.na(long_racf), !is.na(lat_racf),
         !is.na(long_hosp), !is.na(lat_hosp))

# Step 2: Fisheye transform settings
ctr <- tibble(long = 145, lat = -37.8)
pwr <- 0.5

# Step 3: Transform coordinates
transfers_all_tf <- transfers_all %>%
  mutate(
    r_h = sqrt((long_hosp - ctr$long)^2 + (lat_hosp - ctr$lat)^2),
    th_h = atan2(lat_hosp - ctr$lat, long_hosp - ctr$long),
    lr_h = r_h^pwr,
    llong_hosp = lr_h * cos(th_h) + ctr$long,
    llat_hosp = lr_h * sin(th_h) + ctr$lat,
    
    r_a = sqrt((long_racf - ctr$long)^2 + (lat_racf - ctr$lat)^2),
    th_a = atan2(lat_racf - ctr$lat, long_racf - ctr$long),
    lr_a = r_a^pwr,
    llong_racf = lr_a * cos(th_a) + ctr$long,
    llat_racf = lr_a * sin(th_a) + ctr$lat
  )

# Step 4: Load Victoria boundary and convert to data frame
vic <- ne_states(country = "australia", returnclass = "sf") %>%
  filter(name == "Victoria") %>%
  st_transform(crs = 4326)

vic_coords <- st_coordinates(st_geometry(vic)[[1]])

vic_df <- as_tibble(vic_coords) %>%
  rename(long = X, lat = Y) %>%
  mutate(
    r = sqrt((long - ctr$long)^2 + (lat - ctr$lat)^2),
    th = atan2(lat - ctr$lat, long - ctr$long),
    lr = r^pwr,
    llong = lr * cos(th) + ctr$long,
    llat = lr * sin(th) + ctr$lat
  )

# Step 5: Plot map + fisheye network
library(tidyverse)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)

# Reuse previously defined variables: transfers_all_tf, vic_df, ctr, etc.

# Plot with enhancements
ggplot() +
  # Victoria boundary (fisheye)
  geom_polygon(data = vic_df,
               aes(x = llong, y = llat),
               fill = "grey95", color = "grey60", linewidth = 0.3) +
  
  # Transfer lines
  geom_segment(data = transfers_all_tf,
               aes(x = llong_racf, y = llat_racf,
                   xend = llong_hosp, yend = llat_hosp),
               color = "#94B447", alpha = 0.05) +
  
  # Hospital nodes
  geom_point(data = transfers_all_tf,
             aes(x = llong_hosp, y = llat_hosp, color = "Hospital"),
             size = 0.4, alpha = 0.6) +
  
  # RACF nodes
  geom_point(data = transfers_all_tf,
             aes(x = llong_racf, y = llat_racf, color = "RACF"),
             size = 0.4, alpha = 0.6) +
  
  # Legend
  scale_color_manual(name = "Node Type",
                     values = c("Hospital" = "darkgreen", "RACF" = "orange")) +
  
  # Scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  
  coord_fixed() +
  theme_minimal() +
  labs(title = "Fisheye Network of RACF to Hospital Transfers over Victoria")

+ coord_cartesian(
  xlim = c(144.7, 145.4),
  ylim = c(-38.2, -37.5)
)

ggsave("fisheye_network_victoria.png",
       plot = last_plot(), width = 12, height = 10, dpi = 300)

# Properties
library(dplyr)
library(igraph)

# Step 1: Convert to igraph object (again)
library(igraph)

# Create a vertex table from spatial nodes
vertices <- nodes_sf %>%
  st_drop_geometry() %>%
  distinct(name, .keep_all = TRUE)

# Create igraph object
g <- graph_from_data_frame(
  d = df_edges %>% select(from, to, weight),
  vertices = vertices,
  directed = TRUE
)
#Step 2: Degree Centrality
in_deg <- degree(g, mode = "in")
out_deg <- degree(g, mode = "out")
deg <- degree(g)  # total degree

hist(deg, breaks = 50, main = "Degree Distribution", col = "skyblue")
#Step 3: Centrality Measures
btw <- betweenness(g, directed = TRUE)
cls <- closeness(g, mode = "out")

# Top nodes
head(sort(btw, decreasing = TRUE), 10)
head(sort(cls, decreasing = TRUE), 10)

# Step 4: Clustering Coefficient

transitivity(g, type = "global")   # typically low in bipartite-style graphs
transitivity(g, type = "average")

#Step 5: Modularity & Communities
g_undirected <- as.undirected(g, mode = "collapse")
communities <- cluster_louvain(g_undirected)

# Modularity score
modularity(communities)

# Add to vertex attributes
V(g)$community <- membership(communities)

#Step 6: Combine Results with Location Info
V(g)$in_degree <- in_deg
V(g)$out_degree <- out_deg
V(g)$betweenness <- btw
V(g)$closeness <- cls

node_stats <- data.frame(
  name = V(g)$name,
  type = V(g)$type,
  in_degree = V(g)$in_degree,
  out_degree = V(g)$out_degree,
  betweenness = V(g)$betweenness,
  closeness = V(g)$closeness,
  community = V(g)$community
)

# Merge coordinates from nodes_sf
coords <- st_drop_geometry(nodes_sf)
node_stats <- left_join(node_stats, coords, by = c("name", "type"))

head(node_stats)

head(sort(degree(g, mode = "in"), decreasing = TRUE), 10)
