
# =============================================================================
# weather_station_summaries.R
#
# Purpose:
#   Summarize tree crown structure (LAI, tree height, crown base height,
#   crown effective radius) and land cover fractional composition within
#   circular buffers around each NYC weather station used to constrain the
#   neighborhood-scale urban climate model. Buffers are computed at 250 m
#   and 500 m radius around each station.
#
# Inputs:
#   - location_NYC_sites.csv
#       Weather station locations. Must contain columns "latitude_degrees_north" 
#       and "longitude_degrees_north" (WGS84, decimal degrees). Any other 
#       columns (e.g. a station name or ID) are carried through unchanged to the
#       output.
#   - tree_polygon_lai_estimates.gpkg
#       Vector polygon layer of individual tree crowns (produced by
#       get_lai_at_crowns.R), with attribute columns "LAI",
#       "max_tree_height_m" (total tree height), "crown_height_m" (height
#       of the crown overall), and "poly_id".
#   - landcover_nyc_2021_6in.tif
#       Single-band classified land cover raster (Zenodo record 14053441)
#       with integer classes:
#         1 = tree canopy, 2 = grass/shrub, 3 = bare ground, 4 = water,
#         5 = building, 6 = road, 7 = other impervious surface, 8 = railroad
#
# Output:
#   - weather_station_buffer_summaries.csv
#       One row per weather station x buffer radius (250 m, 500 m):
#         - mean_lai, mean_tree_height_m, mean_crown_base_height_m,
#           mean_crown_effective_radius_m: crown area-weighted means,
#           computed only over the tree canopy area within the buffer
#           (all units m or m^2)
#         - frac_cover_*: fractional cover (0-1) of each land cover class
#           within the buffer
#
# Notes / assumptions:
#   - Crown effective radius is computed as sqrt(area / pi) from each
#     crown polygon's actual area, not from the pre-existing
#     crown_diameter_m attribute.
#   - For crowns that straddle a buffer edge, only the portion of the crown
#     polygon that falls inside the buffer is used both to weight and to
#     compute LAI / height / effective radius, so partial crowns contribute
#     proportionally to their overlap with the buffer.
# =============================================================================

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(units)

# ---- File paths and parameters ----------------------------------------------

weather_station_path <- "F:/NYC/flux_sites/location_NYC_sites.csv"
tree_crown_polygon_path <- "tree_polygon_lai_estimates.gpkg"
landcover_raster_path <- "F:/NYC/Tree_Data/crowns/landcover_nyc_2021_6in.tif"
modeled_lai_monthly_path <- "modeled_lai_monthly.csv"
output_summary_path <- "weather_station_buffer_summaries.csv"

buffer_radii_m <- c(250, 500)

landcover_class_labels <- c(
  "1" = "tree_canopy",
  "2" = "grass_shrub",
  "3" = "bare_ground",
  "4" = "water",
  "5" = "building",
  "6" = "road",
  "7" = "other_impervious",
  "8" = "railroad"
)

# ---- Load input datasets -----------------------------------------------------

weather_stations <- read_csv(weather_station_path) |>
  janitor::clean_names() |> 
  st_as_sf(coords = c("longitude_degrees_east", "latitude_degrees_north"), 
           crs = 4326, remove = FALSE) |>
  mutate(station_row_id = row_number())

tree_crowns <- st_read(tree_crown_polygon_path, quiet = TRUE)

landcover <- rast(landcover_raster_path)

# Reproject weather stations into the crown polygon CRS so buffer distances
# can be specified in real-world meters below regardless of the underlying
# CRS units (e.g. if the crown data is in a state plane CRS using feet).
weather_stations <- st_transform(weather_stations, st_crs(tree_crowns))

# ---- Build circular buffers around each station, for each requested radius --

station_buffers <- buffer_radii_m |>
  map(function(radius_m) {
    weather_stations |>
      mutate(
        buffer_radius_m = radius_m,
        geometry = st_buffer(geometry, dist = set_units(radius_m, "m"))
      )
  }) |>
  bind_rows()

# ---- Crown area-weighted tree structure summaries ----------------------------

#' Compute crown area-weighted tree structure metrics within buffers.
#'
#' Intersects tree crown polygons with each buffer polygon, then computes
#' area-weighted mean LAI, tree height, crown base height, and crown
#' effective radius using only the portion of each crown that falls inside
#' the buffer as both the weight and the basis for the effective radius.
#'
#' @param crowns sf polygon layer of tree crowns with columns LAI,
#'   max_tree_height_m, crown_height_m.
#' @param buffers sf polygon layer of buffers, one row per station x
#'   radius, with columns station_row_id and buffer_radius_m.
#' @return Tibble with one row per station_row_id x buffer_radius_m, with
#'   crown area-weighted tree structure metrics.
summarize_tree_structure <- function(crowns, buffers) {
  crowns_near_buffers <- crowns |>
    select(LAI, max_tree_height_m, crown_height_m) |>
    st_filter(buffers, .predicate = st_intersects)

  crown_fragments <- st_intersection(
    crowns_near_buffers,
    buffers |> select(station_row_id, buffer_radius_m)
  ) |>
    st_collection_extract("POLYGON", warn = FALSE)

  crown_fragments$crown_fragment_area_m2 <- as.numeric(st_area(crown_fragments))

  crown_fragments |>
    st_drop_geometry() |>
    mutate(crown_effective_radius_m = sqrt(crown_fragment_area_m2 / pi)) |>
    group_by(station_row_id, buffer_radius_m) |>
    summarize(
      mean_lai = weighted.mean(LAI, crown_fragment_area_m2, na.rm = TRUE),
      mean_tree_height_m = weighted.mean(max_tree_height_m, crown_fragment_area_m2, na.rm = TRUE),
      mean_crown_base_height_m = weighted.mean(max_tree_height_m-crown_height_m, crown_fragment_area_m2, na.rm = TRUE),
      mean_crown_effective_radius_m = weighted.mean(crown_effective_radius_m, crown_fragment_area_m2, na.rm = TRUE),
      .groups = "drop"
    )
}

tree_structure_summary <- summarize_tree_structure(tree_crowns, station_buffers)

# ---- Land cover fractional cover summaries ------------------------------------

#' Compute fractional land cover cover within buffers.
#'
#' Uses exact (area-weighted) raster extraction so that land cover cells
#' only partially covered by a buffer contribute proportionally to the
#' fractional cover totals, rather than being fully included or excluded.
#'
#' @param landcover_raster Single-band SpatRaster of classified land cover
#'   (see class codes in the file header comment).
#' @param buffers sf polygon layer of buffers, one row per station x
#'   radius, with columns station_row_id and buffer_radius_m.
#' @param class_labels Named character vector mapping land cover integer
#'   codes (as strings) to human-readable class names.
#' @return Tibble with one row per station_row_id x buffer_radius_m, and
#'   one frac_cover_* column per land cover class giving fractional cover
#'   (0-1).
summarize_landcover_fractions <- function(landcover_raster, buffers, class_labels) {
  buffers_reproj <- st_transform(buffers, crs(landcover_raster))

  # Only carry the join keys (station_row_id, buffer_radius_m) forward from
  # buffers_reproj, rather than every station attribute column, so this
  # doesn't produce columns that collide (and get .x/.y suffixed) with the
  # same station attributes already present in station_lookup at the final
  # left_join.
  buffer_keys <- buffers_reproj |> dplyr::select(station_row_id, buffer_radius_m)

  # Extract land cover class values in each buffer radius
  # land cover class value and the fraction of that cell covered by the
  # polygon (columns: ID, <layer name>, fraction).
  # output is a list of dataframes, one for each buffer radius
  cover_values <- exactextractr::exact_extract(landcover_raster, buffers_reproj)
  buffer_ind <- 0
  cover_fractions <- lapply(cover_values,
         function(new_df){
           buffer_ind <<- buffer_ind + 1
           cat("\nWorking on data for buffer index", buffer_ind)
           summary_df <- new_df |> group_by(value) |>
                    summarize(total_cover = sum(coverage_fraction, na.rm=TRUE)) |>
                    ungroup() |>
                    mutate(fractional_cover = total_cover / sum(total_cover, na.rm=TRUE)) |>
                    dplyr::select(-total_cover) |>
                    pivot_wider(names_from=value, values_from=fractional_cover, names_prefix="LC_")
           if(nrow(summary_df) != 0)
             return(summary_df |>
                      cbind(buffer_keys[buffer_ind,]))
           return(summary_df)
         }) |>
    bind_rows()
  # Re-order names
  cover_fractions <- cover_fractions |>
    dplyr::select(names(buffer_keys), paste0("LC_", 1:8))
  names(cover_fractions) <- c(names(buffer_keys),
                              "tree_cover", "grass_cover", "bare_ground_cover",
                              "water_cover", "building_cover", "road_cover",
                              "other_impervious_cover", "railroad_cover")
  
  #         1 = tree canopy, 2 = grass/shrub, 3 = bare ground, 4 = water,
  #         5 = building, 6 = road, 7 = other impervious surface, 8 = railroad
  return(cover_fractions)
}

landcover_fraction_summary <- summarize_landcover_fractions(landcover, station_buffers, landcover_class_labels)

# ---- Add seasonal variability from MODIS (fixed, average across city) --------

#' Expand a per-site peak-season LAI value into 12 monthly LAI values.
#'
#' `modeled_lai_monthly.csv` gives a single city-wide-average seasonal curve
#' (deciduous canopy LAI runs from ~0 in winter dormancy up to a summer
#' peak), and each site's `mean_lai` computed above is treated as that
#' site's own peak-season value. Each month's site-level LAI is therefore
#' `mean_lai * (monthly_lai / max(monthly_lai))`, i.e. the fixed seasonal
#' curve is linearly rescaled from 0 to `mean_lai` at every site.
#'
#' @param peak_lai Numeric vector of per-site peak-season LAI values
#'   (the existing `mean_lai` column).
#' @param monthly_lai Numeric vector of length 12, city-average LAI for
#'   each month in calendar order (Jan-Dec).
#' @return Tibble with one row per element of `peak_lai` and columns
#'   `lai_01`..`lai_12`.
scale_lai_seasonally <- function(peak_lai, monthly_lai) {
  stopifnot(length(monthly_lai) == 12)
  seasonal_fraction <- monthly_lai / max(monthly_lai, na.rm = TRUE)
  monthly_lai_matrix <- outer(peak_lai, seasonal_fraction)
  colnames(monthly_lai_matrix) <- sprintf("lai_%02d", 1:12)
  as_tibble(monthly_lai_matrix)
}

modeled_lai_monthly <- read_csv(modeled_lai_monthly_path) |>
  arrange(doy)

stopifnot(nrow(modeled_lai_monthly) == 12)

tree_structure_summary <- tree_structure_summary |>
  bind_cols(scale_lai_seasonally(tree_structure_summary$mean_lai, modeled_lai_monthly$LAI)) |>
  select(-mean_lai)

# ---- Combine summaries and write output ---------------------------------------

station_lookup <- station_buffers |>
  st_drop_geometry() |>
  distinct(station_row_id, buffer_radius_m, .keep_all = TRUE)

weather_station_summary <- station_lookup |>
  left_join(tree_structure_summary, by = c("station_row_id", "buffer_radius_m")) |>
  left_join(landcover_fraction_summary, by = c("station_row_id", "buffer_radius_m"))

write_csv(weather_station_summary |> 
            dplyr::select(-geometry) |> 
            drop_na(lai_01) |> 
            mutate(across(everything(), ~ replace_na(.x, 0))), 
          output_summary_path)
