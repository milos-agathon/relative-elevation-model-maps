# 1. PACKAGES
#------------

install.packages("devtools")

devtools::install_github(
    "ropensci/terrainr"
)

pacman::p_load(
    terrainr, terra,
    sf, osmdata, classInt,
    ggplot2, rayshader
) # add gstat

# 2. DEFINE AREA
#---------------

#-110.621713,43.750439,-110.573660,43.781138

xmin <- -110.621713
ymin <- 43.750439
xmax <- -110.573660
ymax <- 43.781138

bb <- sf::st_sfc(
    sf::st_polygon(
        list(
            cbind(
                c(xmin, xmax, xmax, xmin, xmin),
                c(ymin, ymin, ymax, ymax, ymin)
            )
        )
    ),
    crs = 4326
)

# 3. GET LIDAR DATA
#------------------

dem <- terrainr::get_tiles(
    data = bb,
    output_prefix = "rem",
    side_length = 8e3,
    resolution = 1,
    services = "elevation",
    verbose = TRUE
)

dem_rast <- terra::rast(dem$elevation)
dem_rast <- terra::rast("rem_3DEPElevation_1_1.tif")

# 4. GET RIVER LINE
#------------------

river <- osmdata::opq(
    bbox = bb
) |>
    osmdata::add_osm_feature(
        key = "waterway",
        value = "river"
    ) |>
    osmdata::osmdata_sf()

river_sf <- river$osm_lines |>
    sf::st_intersection(
        bb
    ) |>
    sf::st_union() |>
    sf::st_cast(
        "LINESTRING"
    ) |>
    sf::st_as_sf()

terra::plot(dem_rast)
plot(
    sf::st_geometry(
        river_sf
    ),
    col = "white",
    add = TRUE
)

# 5. EXTRACT ELEVATION VALUES
#----------------------------

dem_rast_agg <- terra::aggregate(
    dem_rast,
    fact = 5
)

river_elev <- terra::extract(
    x = dem_rast_agg,
    y = terra::vect(river_sf),
    xy = TRUE,
    na.rm = TRUE
) |>
    na.omit()

names(river_elev)[2] <- "elevation"
nrow(river_elev)

# 6. DEFINE MODEL
#----------------

idw_model <- gstat::gstat(
    formula = elevation ~ 1,
    locations = ~ x + y,
    data = river_elev,
    nmax = nrow(river_elev)
)

# 7. PREDICT VALUES
#------------------

river_surface <- terra::interpolate(
    dem_rast_agg,
    idw_model,
    crs = terra::crs(dem_rast_agg)
)

# 8. REM
#--------

rem <- dem_rast_agg - river_surface

rem_final <- terra::resample(
    rem, dem_rast
)

# 9. FINAL PREPARATION
#---------------------

rem_df <- as.data.frame(
    rem_final,
    xy = TRUE
)

head(rem_df)

names(rem_df)[3] <- "elevation"

epsilon <- 1e-10

rem_df$elevation_log <- log1p(
    pmax(
        rem_df$elevation, epsilon
    )
)

breaks <- classInt::classIntervals(
    rem_df$elevation_log,
    n = 12,
    style = "fisher"
)$brks

cols <- hcl.colors(
    palette = "Mako",
    12, rev = TRUE
)

pie(rep(
    1, length(cols)
), col = cols)

pal <- cols[c(1, 8:12)]

pie(rep(
    1, length(pal)
), col = pal)

theme_for_the_win <- function() {
    theme_minimal() +
        theme(
            axis.line = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.background = element_rect(
                fill = "white", color = NA
            ),
            plot.margin = unit(
                c(
                    t = 0, r = 0,
                    l = 0, b = 0
                ), "cm"
            )
        )
}

rem_plot <- ggplot(
    rem_df, aes(
        x = x, y = y,
        fill = elevation_log
    )
) +
    geom_raster() +
    scale_fill_gradientn(
        breaks = breaks,
        colours = pal,
        name = ""
    ) +
    theme_for_the_win()

# 11. 3D PLOT
#------------

width <- ncol(rem_final) / 500
height <- nrow(rem_final) / 500

rayshader::plot_gg(
    ggobj = rem_plot,
    width = width,
    height = height,
    windowsize = c(
        width * 75,
        height * 75
    ),
    scale = 75,
    solid = FALSE,
    shadow = FALSE,
    shadow_intensity = 1,
    phi = 87,
    theta = 0,
    zoom = .64,
    multicore = TRUE
)

rayshader::render_camera(
    zoom = .545
)

u <- "https://dl.polyhaven.org/file/ph-assets/HDRIs/hdr/4k/rosendal_plains_2_4k.hdr"
hdri_file <- basename(u)

download.file(
    url = u,
    destfile = hdri_file,
    mode = "wb"
)

rayshader::render_highquality(
    filename = "3d-snake-river.png",
    preview = TRUE,
    light = FALSE,
    environment_light = hdri_file,
    intensity = .85,
    rotate_env = 90,
    parallel = TRUE,
    width = width * 500,
    height = height * 500,
    interactive = FALSE
)
