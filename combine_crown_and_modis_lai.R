
# Load individual tree LAI estimates and MODIS grid
# Generate a dataset for LAI at the 500 m MODIS scale
# Normalize it seasonally based on MODIS seasonality 

library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(phenofit)

rmse <- function(x1, x2){ return(sqrt(mean((x1 - x2)^2, na.rm=TRUE)))}

# Load the overall MODIS LAI dataset 
modis_dir <- "F:/NYC/MODIS/"
output_dir <- "F:/NYC/iTree_LAI_data/"
phenometrics <- terra::rast(paste0(modis_dir, "MODIS_phenology_mean.tif"))
crs(phenometrics) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
phenology <- terra::rast(paste0(modis_dir, "MODIS_LAI_NYC_2024_phenoseries.tif"))
crs(phenology) <- crs(phenometrics)

# Also load the urban MODIS LAI dataset (Dong et al 2025)
modis_urban_lai_dir <- "F:/NYC/MODIS_global_urban_dong_2025/"
# This is produced annually - get a long-term average: 
modis_urban_lai_files <- list.files(modis_urban_lai_dir, full.names=TRUE)
# Create a raster with the valid extent of the other MODIS datasets,
#   reprojected to match the Dong et al 2025 dataset's CRS
extent_raster <- terra::project(phenometrics[[1]], 
                                terra::rast(modis_urban_lai_files[[1]])) |> 
  terra::trim()
modis_urban_lai_all_years <- lapply(modis_urban_lai_files,
                                    function(filepath){
  return( terra::rast(filepath) |> 
            terra::crop(extent_raster) |> 
            terra::project(phenology[[1]]) )                                    
})
# Note - I was interested in tracking afforestation using the MODIS LAI product. Here's an example:
regress_lai_vs_year <- function(){
  # modis_urban_lai_all_years: list of n SpatRasters, each 24 bands (one raster per year)
  # years:   numeric vector length n, in the same order as `modis_urban_lai_all_years`
  years <- 2000:2022
  stopifnot(length(years) == length(modis_urban_lai_all_years))
  print("proceeding")
  
  # per-pixel regression across years, for one band
  regfun <- function(y, x) {
    if (sum(!is.na(y)) < 3) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
    fit <- tryCatch(lm(y ~ x), error = function(e) NULL)
    if (is.null(fit)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
    cf <- summary(fit)$coefficients
    c(r2        = summary(fit)$r.squared,   # R^2
      pvalue    = cf[2, 4],                 # p-value of slope
      slope     = cf[2, 1],                 # slope
      intercept = cf[1, 1])                 # intercept
  }
  
  nb  <- nlyr(modis_urban_lai_all_years[[1]])          # 24
  res <- vector("list", nb)
  
  print("building band stacks")
  for (b in seq_len(nb)) {
    band_stack <- rast(lapply(modis_urban_lai_all_years, function(r) r[[b]]))   # n layers (one per year)
    print("working on regression for band")
    print(b)
    res[[b]]   <- app(band_stack, fun = function(v) regfun(v, years))
  }
  
  # reassemble into four 24-band modis_urban_lai_all_years, one per statistic
  r2        <- rast(lapply(res, function(z) z[[1]]))
  pvalue    <- rast(lapply(res, function(z) z[[2]]))
  slope     <- rast(lapply(res, function(z) z[[3]]))
  intercept <- rast(lapply(res, function(z) z[[4]]))
  
  # preserve original band names
  bn <- names(modis_urban_lai_all_years[[1]])
  names(r2) <- names(pvalue) <- names(slope) <- names(intercept) <- bn
  
  return(list(r2, pvalue, slope, intercept))
}
modis_urban_lai_trends <- regress_lai_vs_year()
# Generate and save a plot showing actual LAI change over time
modis_urban_lai_slopes <- modis_urban_lai_trends[[3]][[1:12]]
names(modis_urban_lai_slopes) <- month.name
modis_urban_lai_plot <- ggplot() + 
  geom_spatraster(data=modis_urban_lai_slopes*22) + 
  scale_fill_gradient2(low = "#550055", mid = "white", high = "#00AA00", 
                       midpoint = 0, limits=c(-0.05, 0.05)*22) +
  facet_wrap(~lyr) + 
  theme_bw() + 
  scale_x_continuous(expand=c(0,0), breaks=c(-74.2, -73.8)) + 
  scale_y_continuous(expand=c(0,0), breaks=c(40.5, 40.7, 40.9)) +
  theme(strip.background = element_rect(fill="white")) + 
  labs(fill="LAI Change")
modis_urban_lai_plot
ggsave("output_plots/modis_urban_lai_plot.png", modis_urban_lai_plot, 
       height=4, width=6)
# Do the same showing R^2 between LAI and year 
modis_urban_lai_rsqd <- modis_urban_lai_trends[[1]][[1:12]]
names(modis_urban_lai_rsqd) <- month.name
modis_urban_lai_rsqd_plot <- ggplot() + 
  geom_spatraster(data=modis_urban_lai_rsqd) + 
  scale_fill_gradient2(low = "white", high = "#FF0000", 
                       limits=c(0, 1.0)) +
  facet_wrap(~lyr) + 
  theme_bw() + 
  scale_x_continuous(expand=c(0,0), breaks=c(-74.2, -73.8)) + 
  scale_y_continuous(expand=c(0,0), breaks=c(40.5, 40.7, 40.9)) +
  theme(strip.background = element_rect(fill="white")) + 
  labs(fill=bquote("LAI-Year " ~R^2* ""))
modis_urban_lai_rsqd_plot
ggsave("output_plots/modis_urban_lai_rsqd_plot.png", modis_urban_lai_rsqd_plot, 
       height=4, width=6)
# Last up, just get an averaged June-July-August value for all years
modis_urban_lai_summer_avg <- lapply(
  modis_urban_lai_all_years,
  function(new_rast){
    return(mean(new_rast[[6:8]],
                na.rm=TRUE))
  }
) |> 
  terra::rast() |> 
  mean(na.rm=TRUE)
names(modis_urban_lai_summer_avg) <- "urban_peak_lai"



# Load LAI data at tree crowns 
lai_points <- st_read("tree_lai_estimates.gpkg") |> 
  st_transform(crs(phenometrics))
# Rasterize LAI data
lai_rasterized <- rasterize(vect(lai_points |> drop_na(leaf_area_imputed_dbh)), 
                            phenometrics, 
                            field = "leaf_area_imputed_dbh", fun = "sum", background = 0)
lai_rasterized <- lai_rasterized / 463.31^2
names(lai_rasterized) <- "LAI"
lai_rasterized[lai_rasterized == 0] <- NA
terra::writeRaster(lai_rasterized, 
                   paste0(output_dir, "/lai_raster.tif"), 
                   overwrite=TRUE)
# Also get a rasterized estimate of total tree area in each grid cell 
total_crown_area <- rasterize(vect(lai_points |> 
                                     mutate(crown_area_m2 = pi*(crown_diameter_m)^2/4)), 
                              phenometrics, 
                              field = "crown_area_m2", fun = "sum", background = 0)
names(total_crown_area) <- "total_crown_area_m2"
total_crown_area[total_crown_area == 0] <- NA
terra::writeRaster(total_crown_area, 
                   paste0(output_dir, "/total_crown_area.tif"), 
                   overwrite=TRUE)
terra::writeRaster(total_crown_area / (463.31^2), 
                   paste0(output_dir, "/fractional_tree_cover.tif"), 
                   overwrite=TRUE)
mean_lai_of_crowns <- lai_rasterized * (463.31^2) / total_crown_area
names(mean_lai_of_crowns) <- "mean_crown_LAI"
terra::writeRaster(mean_lai_of_crowns, 
                   paste0(output_dir, "/mean_lai_of_crowns.tif"), 
                   overwrite=TRUE)

# Compare to the MODIS LAI product
modis_lai_estimate <- terra::rast(paste0(modis_dir, "MODIS_LAI_summer_NYC_mean.tif"))
crs(modis_lai_estimate) <- crs(phenometrics)
names(modis_lai_estimate) <- "MODIS_LAI"
evi_lai_estimate <- phenometrics[["EVI_Amplitude_1"]]
combined_lai <- c(modis_lai_estimate, evi_lai_estimate, lai_rasterized, 
                  modis_urban_lai_summer_avg, total_crown_area, mean_lai_of_crowns)
combined_lai_df <- combined_lai |> 
  as.data.frame()

# Count how many LAI cells were resolved with LiDAR+iTree but not MODIS
#  First, for the MODIS LAI product:
(combined_lai_df |> 
    drop_na(LAI) |> 
    summarize(num_invalid = sum(is.na(MODIS_LAI))))/nrow(combined_lai_df |> 
                                                           drop_na(LAI))
#  Next, for the phenometric EVI_Amplitude_1:
(combined_lai_df |> 
    drop_na(LAI) |> 
    summarize(num_invalid = sum(is.na(EVI_Amplitude_1))))/nrow(combined_lai_df |> 
                                                           drop_na(LAI))
#  Last, for the urban LAI dataset:
(combined_lai_df |> 
    drop_na(LAI) |> 
    summarize(num_invalid = sum(is.na(urban_peak_lai))))/nrow(combined_lai_df |> 
                                                                 drop_na(LAI))

modis_itree_lai_comparison_model <- lm(data=combined_lai_df |> 
                                         drop_na(LAI, MODIS_LAI) |> 
                                         filter(!is.infinite(LAI),
                                                !is.infinite(MODIS_LAI)), 
                                       LAI ~ MODIS_LAI)
summary(modis_itree_lai_comparison_model)
rmse <- function(x1, x2){ return(sqrt(mean((x1 - x2)^2, na.rm=TRUE)))}
rmse(combined_lai_df$LAI, combined_lai_df$MODIS_LAI)

MODIS_LAI_comparison_plot <- ggplot(combined_lai_df) + 
  geom_point(aes(x=MODIS_LAI, LAI)) + 
  geom_abline(intercept=0, slope=1, linetype="dashed", col="black") + 
  geom_abline(intercept=summary(modis_itree_lai_comparison_model)$coefficients[1,1], 
              slope=summary(modis_itree_lai_comparison_model)$coefficients[2,1], 
              linetype="solid", col="red") + 
  annotate("text", x=0.5, y=7.5, parse=TRUE, hjust=0,
           label = paste("R^2 == ", 
                         round(summary(modis_itree_lai_comparison_model)$adj.r.squared, 2))) + 
  annotate("text", x=0.5, y=7, hjust=0,
           label = paste("RMSE = ", 
                         round(rmse(combined_lai_df$LAI, combined_lai_df$MODIS_LAI), 2))) + 
  theme_bw() + 
  labs(x = "MODIS LAI Estimate", 
       y = "iTree-derived LAI Estimate") + 
  scale_x_continuous(limits=c(0,8), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0,8), expand=c(0,0))
MODIS_LAI_comparison_plot
ggsave("output_plots/MODIS_LAI_comparison_plot.png", MODIS_LAI_comparison_plot,
       height=4, width=4)


# Repeat the above, comparing vs. the urban LAI dataset
modis_urban_itree_lai_comparison_model <- lm(data=combined_lai_df |> 
                                               drop_na(LAI, urban_peak_lai) |> 
                                               filter(!is.infinite(LAI),
                                                      !is.infinite(urban_peak_lai)), 
                                             LAI ~ urban_peak_lai)
summary(modis_urban_itree_lai_comparison_model)
rmse <- function(x1, x2){ return(sqrt(mean((x1 - x2)^2, na.rm=TRUE)))}
rmse(combined_lai_df$LAI, combined_lai_df$urban_peak_lai)

MODIS_urban_LAI_comparison_plot <- ggplot(combined_lai_df) + 
  geom_point(aes(x=urban_peak_lai, mean_crown_LAI), alpha=0.05) + 
  geom_abline(intercept=0, slope=1, linetype="dashed", col="black") + 
  geom_abline(intercept=summary(modis_itree_lai_comparison_model)$coefficients[1,1], 
              slope=summary(modis_itree_lai_comparison_model)$coefficients[2,1], 
              linetype="solid", col="red") + 
  annotate("text", x=0.5, y=7.5, parse=TRUE, hjust=0,
           label = paste("R^2 == ", 
                         round(summary(modis_itree_lai_comparison_model)$adj.r.squared, 2))) + 
  annotate("text", x=0.5, y=7, hjust=0,
           label = paste("RMSE = ", 
                         round(rmse(combined_lai_df$mean_crown_LAI, combined_lai_df$urban_peak_lai), 2))) + 
  theme_bw() + 
  labs(x = "MODIS LAI Estimate", 
       y = "iTree-derived LAI Estimate") + 
  scale_x_continuous(limits=c(0,8), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0,8), expand=c(0,0))
MODIS_urban_LAI_comparison_plot
ggsave("output_plots/MODIS_urban_LAI_comparison_plot.png", MODIS_urban_LAI_comparison_plot,
       height=4, width=4)

modis_histogram_comparison_plot <- ggplot() + 
  geom_density(data=combined_lai_df |> drop_na(LAI) |> filter(!is.na(MODIS_LAI)), 
               aes(x=LAI), col="green3") + 
  geom_density(data=combined_lai_df |> drop_na(LAI) |> filter(!is.na(EVI_Amplitude_1)), 
               aes(x=LAI), col="orange") + 
  geom_density(data=combined_lai_df, 
               aes(x=LAI), col="black") + 
  annotate("text", x=2, y=0.7, hjust=0, label="MCD15A3H MODIS LAI Product", col="green3") + 
  annotate("text", x=2, y=0.65, hjust=0, label="MCD12Q2 MODIS EVI Phenology Amplitude", col="orange") + 
  annotate("text", x=2, y=0.6, hjust=0, label="LiDAR + iTree ITC LAI Product", col="black") + 
  theme_bw() + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab("LAI") + 
  ylab("Frequency")
modis_histogram_comparison_plot
ggsave("output_plots/modis_histogram_comparison_plot.png", modis_histogram_comparison_plot,
       height=4, width=4)


# Fit a phenology curve to the AVERAGE LAI values in each bimonthly period 
lai_phenology <- phenology |> 
  as.data.frame() |> 
  as.matrix() |> 
  colMeans(na.rm=TRUE) |> 
  t() |> 
  as.data.frame() |> 
  pivot_longer(1:25, 
               names_to="month", 
               values_to="LAI") |> 
  mutate(month = as.numeric(str_split_i(month, "_", 1)))
phenofit_results <- phenofit::curvefit(lai_phenology$LAI, (lai_phenology$month)*365/24+1)
applyPhenoModel <- function(t, pheno_model)
{
  pheno_model$par[1] + (pheno_model$par[2] - pheno_model$par[1]) * (1/(1 + exp(-pheno_model$par[4] * (t - pheno_model$par[3]))) + 1/(1 + exp(pheno_model$par[6] * (t - pheno_model$par[5]))) - 1)
}
modeled_lai_vals <- data.frame(LAI = applyPhenoModel(1:365, phenofit_results$model$Beck),
                               doy = 1:365)

# Compare to average dates in MODIS phenology product for high-LAI pixels:
phenometrics |> 
  c(modis_lai_estimate) |>
  as.data.frame() |> 
  drop_na(MODIS_LAI, Greenup_1) |> 
  filter(MODIS_LAI > 2) |> 
  as.matrix() |> 
  colMeans()

lai_seasonality_plot <- ggplot() + 
  geom_point(data=lai_phenology,
             aes(x=month/2, y=LAI)) + 
  geom_line(data=modeled_lai_vals,
            aes(x=doy/365*12, y=LAI)) + 
  geom_vline(xintercept=phenofit_results$model$Beck$par[[3]]/365*12, linetype="dashed") + 
  geom_vline(xintercept=phenofit_results$model$Beck$par[[5]]/365*12, linetype="dashed") + 
  theme_bw() + 
  xlab("Month") + 
  ylab("LAI") + 
  scale_x_continuous(limits=c(0,12), expand=c(0,0), breaks=(0:6)*2) +
  scale_y_continuous(limits=c(0,2.5), expand=c(0,0))
lai_seasonality_plot
ggsave("output_plots/lai_seasonality_plot.png", lai_seasonality_plot,
       height=4, width=6)


# How well does crown area alone explain LAI? 
lai_and_area_cell_df <- c(lai_rasterized, 
                          total_crown_area) |> 
  as.data.frame()
lai_crown_area_model <- lm(data=lai_and_area_cell_df,
                           LAI ~ total_crown_area_m2)
lai_and_area_cell_df <- lai_and_area_cell_df |> 
  mutate(LAI_predicted = predict(lai_crown_area_model, pick(LAI, total_crown_area_m2)))
LAI_from_area_plot <- ggplot(lai_and_area_cell_df) + 
  geom_point(aes(x=total_crown_area_m2/4e4, y=LAI), alpha=0.1) + 
  geom_abline(intercept=summary(lai_crown_area_model)$coefficients[1,1], 
              slope=summary(lai_crown_area_model)$coefficients[2,1]*1e6, 
              linetype="solid", col="red") + 
  annotate("text", x=2.5, y=5.5, parse=TRUE, hjust=0,
           label = paste("R^2 == ", 
                         round(summary(lai_crown_area_model)$adj.r.squared, 2))) + 
  annotate("text", x=2.5, y=5, hjust=0,
           label = paste("RMSE = ", 
                         round(rmse(lai_and_area_cell_df$LAI, lai_and_area_cell_df$LAI_predicted), 2))) + 
  theme_bw() + 
  scale_x_continuous(limits=c(0,5), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0,6), expand=c(0,0)) + 
  labs(
    x = expression("Grid-scale Total Crown Area (" * km^2 * ")"),
    y = expression("LAI")
  )
LAI_from_area_plot
ggsave("output_plots/LAI_from_area_plot.png", LAI_from_area_plot,
       height=4, width=4)

