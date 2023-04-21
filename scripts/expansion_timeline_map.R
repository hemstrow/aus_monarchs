#Simple plot for Robinson projection of globe

library(rnaturalearth)
library(ggplot2)
library(sf)
library(ggspatial)


world <- ne_countries(scale = 'small', returnclass = 'sf')

polygon <- st_polygon(x = list(rbind(c(-0.0001, 90),
                                     c(0, 90),
                                     c(0, -90),
                                     c(-0.0001, -90),
                                     c(-0.0001, 90)))) |>
  st_sfc() |>
  st_set_crs(4326)

sf::sf_use_s2(FALSE)


world <- world |> st_difference(polygon)

world_robin <- st_transform(world, crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

sf::sf_use_s2(TRUE)

world <- ggplot() +
  geom_sf(data = world_robin, fill= "azure2")  + #background color for land
  coord_sf(xlim = c(-70e5, 100e5),
           ylim = c(-60e5, 58e5)) +
  theme_bw() +
  annotation_scale(text_cex = 1.5) + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24)) +
  annotation_north_arrow(style = north_arrow_fancy_orienteering(text_size = 20), location = "bl", 
                         pad_y = grid::unit(.05, "native"))

world
ggsave("results/map.jpg", world, "jpg", width = 15, height = 11) # note, need to annotate this with ranges and arrows to produce final image

