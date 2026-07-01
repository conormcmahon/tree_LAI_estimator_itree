
library(tidyverse)
library(terra)

data_dir <- "F:/NYC/MODIS/"

# Load LAI and seasonality data (SOS and EOS)
summer_LAI <- terra::rast(paste0(data_dir, "MODIS_LAI_summer_NYC_mean.tif"))
phenometrics <- terra::rast(paste0(data_dir, "MODIS_phenology_mean.tif"))

# Load phenoseries
phenoseries_2024 <- terra::rast(paste0(data_dir, "MODIS_LAI_NYC_2024_phenoseries.tif"))
names(phenoseries_2024) <- paste0("LAI_2024_period_", sprintf("%02d", 1:(dim(phenoseries_2024)[[3]])))

# Combine data and get into data frame format
MODIS_df <- c(summer_LAI, 
              phenometrics, 
              phenoseries_2024) %>%
  as.data.frame(xy=TRUE)
MODIS_df_long <- MODIS_df %>%
  pivot_longer(29:53, names_prefix="LAI_2024_period_", 
               names_to="period", values_to="period_Lai") %>%
  mutate(period = as.numeric(period))


# Visualize spread in Mid-greenup and Mid-greendown dates (inc. with high/low LAI sites)
# Visualize spread in Mid-greenup and Mid-greendown dates (inc. with high/low LAI sites)
ggplot() + 
  geom_density(data=MODIS_df %>% 
                 drop_na(Lai, Greenup_1), 
               aes(x=Greenup_1), linetype="dashed", col="lightgreen") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(Greenup_1), 
               aes(x=Greenup_1), col="lightgreen") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(Lai, MidGreenup_1), 
               aes(x=MidGreenup_1), linetype="dashed", col="green1") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(MidGreenup_1), 
               aes(x=MidGreenup_1), col="green1") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(Lai, Maturity_1), 
               aes(x=Maturity_1), linetype="dashed", col="green3") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(Maturity_1), 
               aes(x=Maturity_1), col="green3") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(Lai, Senescence_1), 
               aes(x=Senescence_1), linetype="dashed", col="orange1") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(Senescence_1), 
               aes(x=Senescence_1), col="orange1") +    
  geom_density(data=MODIS_df %>% 
                 drop_na(Lai, Senescence_1), 
               aes(x=MidGreendown_1), linetype="dashed", col="red") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(MidGreendown_1), 
               aes(x=MidGreendown_1), col="red") +      
  geom_density(data=MODIS_df %>% 
                 drop_na(Lai, Dormancy_1), 
               aes(x=Dormancy_1), linetype="dashed", col="maroon") + 
  geom_density(data=MODIS_df %>% 
                 drop_na(Dormancy_1), 
               aes(x=Dormancy_1), col="maroon") +         
  geom_smooth(data=MODIS_df_long %>% 
                drop_na(period_Lai, period) %>%
                filter(Lai > 2), 
              aes(x=(period-1)*365/24, period_Lai/50), col="black") + 
  theme_bw()

# Visualize relationship between LAI and greenup / greendown 
ggplot(data=MODIS_df %>%
         drop_na(Lai, MidGreenup_1)) + 
  geom_density_2d_filled(aes(x=Lai, y=MidGreenup_1))
ggplot(data=MODIS_df %>%
         drop_na(Lai, MidGreendown_1)) + 
  geom_density_2d_filled(aes(x=Lai, y=MidGreendown_1))




