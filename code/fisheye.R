library(fisheye)
library(sf)
library(tidyverse)
load("data/lgas.rda")
# Really important to remove the empty polygons
lgas <- lgas %>%
  filter(!st_is_empty(.))
lga_nc <- st_transform(lgas, 3857) %>%
  na.omit()
lga_ncfe <- fisheye(lga_nc, centre = lga_nc[44, ], method = 'log', k = 4)
plot(st_geometry(lga_ncfe), col = "grey70", lwd = .2)

ggplot(lga_ncfe) + geom_sf()

lga_ncfe_df <- as.data.frame(st_coordinates(lga_ncfe))
ggplot(lga_ncfe_df, aes(x=X, y=Y, group=interaction(L1, L2, L3))) +
  geom_polygon(fill=NA, colour="grey80") +
  theme_map()

# Using the health regions
load("data/dhhs_sm.rda")

dhhs_nc <- st_transform(dhhs_sm, 3857)
# Missing CRS - could add but let's try to do the 
# transform manually so can do on all alements
# map, points, lines

dhhs_df <- as_Spatial(dhhs_sm) %>% fortify()
ggplot(dhhs_df, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill=NA, colour="grey80") +
  theme_map()

# Centre: 145, -37.8
ctr <- tibble(long=145, lat=-37.8)

pwr <- 0.5
dhhs_df_tf <- dhhs_df %>%
  mutate(r = sqrt((long-ctr$long)^2 + (lat-ctr$lat)^2),
         th = atan2(lat-ctr$lat, long-ctr$long)) %>%
  mutate(lr = r^pwr) %>%
  mutate(llong = lr * cos(th) + ctr$long,
         llat = lr * sin(th) + ctr$lat)
#ggplot(dhhs_df_log, aes(x=r, y=lr)) + geom_point()
#ggplot(dhhs_df_log, aes(x=th, y=sin(th))) + geom_point()

hospital_tf <- hospital %>% 
  rename(long = longitude,
         lat = latitude) %>%
  mutate(r = sqrt((long-ctr$long)^2 + (lat-ctr$lat)^2),
         th = atan2(lat-ctr$lat, long-ctr$long)) %>%
  mutate(lr = r^pwr) %>%
  mutate(llong = lr * cos(th) + ctr$long,
         llat = lr * sin(th) + ctr$lat)

transfers_all_tf <- transfers_all %>% 
  mutate(r_h = sqrt((long_hosp-ctr$long)^2 + (lat_hosp-ctr$lat)^2),
         th_h = atan2(lat_hosp-ctr$lat, long_hosp-ctr$long),
         r_a = sqrt((long_racf-ctr$long)^2 + (lat_racf-ctr$lat)^2),
         th_a = atan2(lat_racf-ctr$lat, long_racf-ctr$long)) %>%
  mutate(lr_h = r_h^pwr, lr_a = r_a^pwr) %>%
  mutate(llong_hosp = lr_h * cos(th_h) + ctr$long,
         llat_hosp = lr_h * sin(th_h) + ctr$lat, 
         llong_racf = lr_a * cos(th_a) + ctr$long,
         llat_racf = lr_a * sin(th_a) + ctr$lat)

# Check the movement of points
#ggplot(dhhs_df_tf) + 
#  geom_segment(aes(x=long, y=lat, 
#                   xend=llong, yend=llat), alpha=0.2) +
#  geom_point(aes(x=long, y=lat), colour="red", size=0.5) +
#  geom_point(aes(x=llong, y=llat), colour="orange", size=0.5)
  
ggplot() +
  geom_polygon(data=dhhs_df_tf, 
               aes(x=llong, y=llat, group=group), 
               fill=NA, colour="grey80") +
  geom_segment(data=transfers_all_tf,
                 aes(x=llong_racf, xend=llong_hosp,
                     y=llat_racf, yend=llat_hosp,
                     group=NA),
                 colour="#94B447", alpha=0.05) +
  geom_point(data=hospital_tf, aes(x=llong, y=llat), 
             colour="darkgreen", alpha=0.2) +
  theme_map() 
  