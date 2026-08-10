
# =============================================================================
# summarize_modis_lai.R
#
# Purpose:
#   Exploratory visualization of MODIS-derived LAI and phenology (greenup /
#   greendown) data for NYC: spread of phenology dates, relationship between
#   peak LAI and phenology timing, and the seasonal LAI curve implied by the
#   2024 phenoseries raster.
#
# Inputs (all under data_dir, from MODIS products, see Zenodo/NASA sources):
#   - MODIS_LAI_summer_NYC_mean.tif: single-band summer-mean LAI raster
#     (band "Lai").
#   - MODIS_phenology_mean.tif: multi-band phenometrics raster (greenup,
#     mid-greenup, maturity, senescence, mid-greendown, dormancy dates, EVI
#     amplitude, etc).
#   - MODIS_LAI_NYC_2024_phenoseries.tif: multi-band raster with one LAI
#     estimate per ~15-day period across 2024.
#
# Output:
#   This script produces plots for interactive/exploratory review only; it
#   does not write any files to disk.
# =============================================================================

library(tidyverse)
library(terra)

data_dir <- "F:/NYC/MODIS/"

# Load LAI and seasonality data (SOS and EOS)
summer_lai <- terra::rast(paste0(data_dir, "MODIS_LAI_summer_NYC_mean.tif"))
phenometrics <- terra::rast(paste0(data_dir, "MODIS_phenology_mean.tif"))

# Load phenoseries
phenoseries_2024 <- terra::rast(paste0(data_dir, "MODIS_LAI_NYC_2024_phenoseries.tif"))
names(phenoseries_2024) <- paste0("LAI_2024_period_", sprintf("%02d", 1:(dim(phenoseries_2024)[[3]])))

# Combine data and get into data frame format
modis_df <- c(summer_lai,
              phenometrics,
              phenoseries_2024) %>%
  as.data.frame(xy=TRUE)
# Select phenoseries columns by name rather than a hardcoded column-index
# range, since that range depends on how many bands summer_lai/phenometrics
# contain and would silently mis-select columns if those inputs change.
modis_df_long <- modis_df %>%
  pivot_longer(starts_with("LAI_2024_period_"), names_prefix="LAI_2024_period_",
               names_to="period", values_to="period_Lai") %>%
  mutate(period = as.numeric(period))


# Visualize spread in Mid-greenup and Mid-greendown dates (inc. with high/low LAI sites).
# For each phenometric, the dashed density is restricted to pixels with a
# valid summer LAI value (Lai), while the solid density uses all pixels with
# a valid date for that phenometric, so the two lines show whether having
# tree cover (i.e. a resolvable LAI) shifts the phenology date distribution.
ggplot() +
  geom_density(data=modis_df %>% 
                 drop_na(Lai, Greenup_1), 
               aes(x=Greenup_1), linetype="dashed", col="lightgreen") + 
  geom_density(data=modis_df %>% 
                 drop_na(Greenup_1), 
               aes(x=Greenup_1), col="lightgreen") + 
  geom_density(data=modis_df %>% 
                 drop_na(Lai, MidGreenup_1), 
               aes(x=MidGreenup_1), linetype="dashed", col="green1") + 
  geom_density(data=modis_df %>% 
                 drop_na(MidGreenup_1), 
               aes(x=MidGreenup_1), col="green1") + 
  geom_density(data=modis_df %>% 
                 drop_na(Lai, Maturity_1), 
               aes(x=Maturity_1), linetype="dashed", col="green3") + 
  geom_density(data=modis_df %>% 
                 drop_na(Maturity_1), 
               aes(x=Maturity_1), col="green3") + 
  geom_density(data=modis_df %>% 
                 drop_na(Lai, Senescence_1), 
               aes(x=Senescence_1), linetype="dashed", col="orange1") + 
  geom_density(data=modis_df %>% 
                 drop_na(Senescence_1), 
               aes(x=Senescence_1), col="orange1") +    
  geom_density(data=modis_df %>% 
                 drop_na(Lai, Senescence_1), 
               aes(x=MidGreendown_1), linetype="dashed", col="red") + 
  geom_density(data=modis_df %>% 
                 drop_na(MidGreendown_1), 
               aes(x=MidGreendown_1), col="red") +      
  geom_density(data=modis_df %>% 
                 drop_na(Lai, Dormancy_1), 
               aes(x=Dormancy_1), linetype="dashed", col="maroon") + 
  geom_density(data=modis_df %>% 
                 drop_na(Dormancy_1), 
               aes(x=Dormancy_1), col="maroon") +         
  geom_smooth(data=modis_df_long %>% 
                drop_na(period_Lai, period) %>%
                filter(Lai > 2), 
              aes(x=(period-1)*365/24, period_Lai/50), col="black") + 
  theme_bw()

# Visualize relationship between LAI and greenup / greendown 
ggplot(data=modis_df %>%
         drop_na(Lai, MidGreenup_1)) + 
  geom_density_2d_filled(aes(x=Lai, y=MidGreenup_1))
ggplot(data=modis_df %>%
         drop_na(Lai, MidGreendown_1)) + 
  geom_density_2d_filled(aes(x=Lai, y=MidGreendown_1))




