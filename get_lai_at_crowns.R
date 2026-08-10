
# =============================================================================
# get_lai_at_crowns.R
#
# Purpose:
#   Estimate per-tree Leaf Area Index (LAI) at individual LiDAR-delineated
#   crowns, using allometric equations from Nowak (1996) (and a comparison
#   equation from Timilsina et al.) applied to crown height, crown diameter,
#   DBH (measured or genus-imputed), and a genus-specific shading
#   coefficient. Also fits/validates the supporting DBH-imputation model and
#   produces diagnostic plots comparing LAI under several modeling choices
#   (measured vs. imputed DBH, with vs. without species information).
#
# Inputs:
#   - tree_ndvi_median_jja_jf_2021polygons.rds: per-crown species ID /
#     genus classification and NDVI summaries.
#   - treeobjects_2021_nyc.csv: per-crown LiDAR morphometrics (height,
#     crown radius, DBH where measured).
#   - crown_height_model.rds: fitted model predicting crown base height
#     from total height, crown width, and genus.
#   - species_shading_coefficients.csv: genus-level shading coefficients.
#   - top_25_genera_nyc_ufia_2024.csv: UFIA urban forestry survey estimates,
#     used only for an external validation comparison.
#   - nyc_class_tree_genus_polygons_v2.gpkg: crown polygon geometries.
#
# Outputs:
#   - tree_lai_estimates.csv: one row per crown, with LAI estimates and all
#     intermediate allometric variables (no geometry).
#   - tree_lai_estimates.gpkg / tree_polygon_lai_estimates.gpkg: the same
#     data joined to crown centroid / crown polygon geometry, respectively.
#   - output_plots/*.png: diagnostic plots (DBH model fit, LAI sensitivity
#     to imputed DBH and lumped species, LAI vs. DBH/crown diameter/height
#     by genus, comparison to UFIA genus-level totals).
#
# To try tomorrow...
# Clean up all code
# Document code
# Functionalize things more (esp different test cases)
# Final validation vs. actual iTree?
# More validation vs. DBH estimate for Nowak, both for Timilsina, and DBH for UTD
# Simple comparisons of how things would change if all trees shifted genera?
# (reach, later) Simple comparisons of how things would change if all trees grew a bit?
# Write up some kind of small report about allometry for large trees / concern with Nowak approach
# =============================================================================

library(tidyverse)
library(sf)
library(lme4)
library(terra)
library(janitor)

# Root mean square error function
rmse <- function(x1, x2){ return(sqrt(mean((x1 - x2)^2, na.rm=TRUE)))}

# Load data files with tree species and with crown geometry (height, radius)
tree_species_info <- read_rds("F:/NYC/Tree_Data/crowns/tree_ndvi_median_jja_jf_2021polygons.rds") |>
  janitor::clean_names() |> 
  group_by(poly_id) |> 
  slice_head(n=1)
crown_morphometrics <- read_csv("F:/NYC/Tree_Data/crowns/treeobjects_2021_nyc.csv") |> 
  janitor::clean_names() |> 
  mutate(poly_id = 1:n())
# Combine the species and morphometrics datasets
tree_data <- tree_species_info |>  
  merge(crown_morphometrics, by="poly_id")
# Clean up some confusing names
tree_data <- tree_data |> 
  mutate(dbh_cm = dbh*2.54,
         crown_diameter_m = radius*2*0.3048,
         max_tree_height_m = height*0.3048,
         ground_area_projected_m_sq = pi*crown_diameter_m^2/4) |> 
  dplyr::select(-c(dbh, tnc_shape_area, height_max_ft, radius, height, shape_length, shape_area))
# Correct implausible tree heights and DBH
tree_data <- tree_data |> 
  mutate(max_tree_height_m = case_when(
    max_tree_height_m < 1 ~ 1,
    max_tree_height_m > 40 ~ 40,
    TRUE ~ max_tree_height_m
  )) |> 
  mutate(dbh_cm = case_when(
    dbh_cm < 1 ~ 1,
    dbh_cm > 200 ~ 200,
    TRUE ~ dbh_cm
  )) 
rm(crown_morphometrics, tree_species_info)

# Estimate DBH for trees where it wasn't measured directly 
# Some DBH values are implausible - remove these
dbh_model <- lmer(
  dbh_cm ~ max_tree_height_m + crown_diameter_m + 0 + (max_tree_height_m | genus_merged),
  data = tree_data |> drop_na(dbh_cm, max_tree_height_m, crown_diameter_m, genus_predicted),
  REML = TRUE
)
tree_data <- tree_data |> 
  mutate(dbh_pred = predict(dbh_model, 
                            newdata=tree_data,
                            allow.new.levels=TRUE))
tree_data <- tree_data |> 
  mutate(dbh_pred = ifelse(dbh_pred < 1, 1, dbh_pred),
         dbh_pred = ifelse(dbh_pred > 200, 200, dbh_pred))
dbh_model_assessment <- lm(data=tree_data, dbh_cm ~ dbh_pred)
summary(dbh_model_assessment)
dbh_rmse <- rmse(tree_data$dbh_cm, tree_data$dbh_pred)
dbh_rmse
dbh_estimate_plot <- ggplot()  + 
  theme_minimal() + 
  geom_density_2d_filled(data=tree_data, 
                         aes(x=dbh_pred, y=dbh_cm), 
                         bins=10) + 
  geom_point(data=tree_data |> 
               filter(genus_merged %in% c("Quercus", "Acer", "Prunus", "Platanus", "Pyrus")) |> 
               group_by(genus_merged) |> 
               sample_n(size=100),
             aes(x=dbh_pred, y=dbh_cm, col=genus_merged)) + 
  geom_abline(intercept=0, slope=1, linetype="dashed") + 
  geom_abline(intercept=summary(dbh_model_assessment)$coefficients[1,1], 
              slope=summary(dbh_model_assessment)$coefficients[2,1], 
              linetype="solid", col="red") + 
  annotate("text", x=5, y=95, parse=TRUE, hjust=0,
           label = paste("R^2 == ", 
                         round(summary(dbh_model_assessment)$adj.r.squared, 2))) + 
  annotate("text", x=5, y=85, hjust=0,
           label = paste("RMSE = ", 
                         round(dbh_rmse, 2))) + 
  scale_fill_manual(values = colorRampPalette(c("#FFFFFF33", "#00008BFF"), alpha = TRUE)(10)) +
  scale_x_continuous(limits=c(0, 100), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0, 100), expand=c(0,0)) + 
  guides(fill = "none") + 
  theme(legend.position = "inside",
        legend.position.inside = c(1.0, 0.0),
        legend.justification = c(1, 0)) + 
  labs(x = "DBH Predicted",
       y = "DBH Measured", 
       color = "Genus")
dbh_estimate_plot
ggsave("output_plots/dbh_estimate_plot.png", dbh_estimate_plot, height=4, width=4)

# Load model to estimate crown base height from total tree height, crown diameter, and genus
crown_height_model <- readRDS("crown_height_model.rds")
# Apply model to get crown base height
tree_data <- tree_data |> 
  mutate(crown_height_ft = predict(crown_height_model,
                                   newdata=tree_data |> 
                                     mutate(max_height = max_tree_height_m/0.3048,
                                            avg_crown_width = crown_diameter_m/0.3048,
                                            genus = genus_merged),
                                   allow.new.levels=TRUE)) |> 
  # Remove implausible crown heights (clamp to between 1 m and the tree's
  # own total height, converting the model's feet output to meters)
  mutate(crown_height_m = case_when(
    crown_height_ft < 1/0.3048 ~ 1,
    crown_height_ft > max_tree_height_m/0.3048 ~ max_tree_height_m,
    TRUE ~ crown_height_ft*0.3048
  )) |>
  dplyr::select(-crown_height_ft)
# Also make a version which has no species information: 
tree_data <- tree_data |> 
  mutate(crown_height_lumped_spp_ft = predict(crown_height_model,
                                           newdata=tree_data |> 
                                             mutate(max_height = max_tree_height_m/0.3048,
                                                    avg_crown_width = crown_diameter_m/0.3048,
                                                    genus = ""),
                                           allow.new.levels=TRUE)) |> 
  # Remove implausible crown heights
  # NOTE: the lower-bound threshold here (1) is compared directly against
  # crown_height_lumped_spp_ft, i.e. it clamps at 1 FOOT. The equivalent
  # branch above for crown_height_m clamps at 1/0.3048 ft, i.e. 1 METER.
  # These two blocks look like they're meant to apply the same clamp policy
  # - flagging the mismatch rather than silently changing which one is
  # "correct", since either could be intentional and this affects the
  # already-published tree_lai_estimates.csv.
  mutate(crown_height_lumped_spp = case_when(
    crown_height_lumped_spp_ft < 1 ~ 1,
    crown_height_lumped_spp_ft > max_tree_height_m/0.3048 ~ max_tree_height_m,
    TRUE ~ crown_height_lumped_spp_ft*0.3048
  )) |>
  dplyr::select(-crown_height_lumped_spp_ft)

# Sanity check for how many cases fall outside of expected geometry bounds
# We expect crown_height_m / crown_diameter_m to range from about 0.5 to 2.0
tree_data |> 
  mutate(total_count = n(),
         height_width_ratio = crown_height_m / crown_diameter_m,
         crown_ratio_bin = case_when(
           height_width_ratio < 0.5 ~ "      HWR < 0.5",
           height_width_ratio < 1.0 ~ "0.5 < HWR < 1.0",
           height_width_ratio < 2.0 ~ "1.0 < HWR < 2.0",
           TRUE ~ "2.0 < HWR"
         )) |> 
  group_by(factor(crown_ratio_bin, levels=c("      HWR < 0.5", "0.5 < HWR < 1.0", "1.0 < HWR < 2.0", "2.0 < HWR"))) |> 
  summarize(count=n(),
            total_count=first(total_count)) |> 
  mutate(fraction=count/total_count)
ggplot(tree_data) + 
  geom_histogram(aes(x=crown_height_m / crown_diameter_m)) + 
  geom_vline(xintercept=c(0.5, 1.0, 2.0), col="red") + 
  scale_x_continuous(limits=c(0, 5), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme_bw() 

# Load species-specific data on shading coefficients 
shading_coeff <- read_csv("species_shading_coefficients.csv")
# NOTE - Sophora japonica was renamed to Styphnolobium japonica:
shading_coeff <- shading_coeff |> 
  mutate(genus = ifelse(genus=="Sophora", "Styphnolobium", genus))
shading_coeff_genus_average <- shading_coeff |> 
  group_by(genus) |> 
  summarize(shading_coeff_mean = mean(shading_coefficient, na.rm=TRUE),
            shading_coeff_sd = sd(shading_coefficient, na.rm=TRUE)) |> 
  mutate(genus_merged = genus) |> 
  dplyr::select(-genus)
# Add shade information to tree dataframe
tree_data <- tree_data |> 
  left_join(shading_coeff_genus_average)
# For genera without shading coefficient species-specific offsets, just set it to zero for now
tree_data <- tree_data |> 
  replace_na(list("shading_coeff_mean" = 0))

# Get list of genera included in NYC classifier
classified_genera <- unique((tree_data |> drop_na(genus_predicted))$genus_predicted) |> 
  sort()
labelled_genera <- unique((tree_data |> drop_na(genus_ref))$genus_ref) |> 
  sort()
fraction_trees_classified <- nrow(tree_data |> drop_na(genus_merged))/nrow(tree_data)

# Estimate Leaf Area and Leaf Area Index for each tree
tree_data <- tree_data |> 
  mutate(height_width_ratio = crown_height_m / crown_diameter_m,
         height_width_ratio_lumped_spp = crown_height_lumped_spp / crown_diameter_m,
         # Correction for cases with crowns with very wide or very skinny crowns
         height_for_lai = case_when(
           height_width_ratio > 2   ~ crown_diameter_m * 2,
           TRUE                     ~ crown_height_m
         ),
         crown_diameter_for_lai = case_when(
           height_width_ratio < 0.5 ~ crown_height_m * 2,
           TRUE                     ~ crown_diameter_m
         ),
         # Correction for exceptionally large crowns
         height_for_lai = case_when(
           height_for_lai > 12 ~ 12,
           height_for_lai < 1  ~ 1,
           TRUE                ~ height_for_lai
         ),
         crown_diameter_for_lai = case_when(
           crown_diameter_for_lai > 14 ~ 14,
           crown_diameter_for_lai < 1  ~ 1,
           TRUE                        ~ crown_diameter_for_lai
         ),
         # Now again for lumped spp
         height_for_lai_lumped_spp = case_when(
           height_width_ratio_lumped_spp > 2   ~ crown_diameter_m * 2,
           TRUE                     ~ crown_height_lumped_spp
         ),
         crown_diameter_for_lai_lumped_spp = case_when(
           height_width_ratio_lumped_spp < 0.5 ~ crown_height_lumped_spp * 2,
           TRUE                     ~ crown_diameter_m
         ),
         # Correction for exceptionally large crowns
         height_for_lai_lumped_spp = case_when(
           height_for_lai_lumped_spp > 12 ~ 12,
           height_for_lai_lumped_spp < 1  ~ 1,
           TRUE                ~ height_for_lai_lumped_spp
         ),
         crown_diameter_for_lai_lumped_spp = case_when(
           crown_diameter_for_lai_lumped_spp > 14 ~ 14,
           crown_diameter_for_lai_lumped_spp < 1  ~ 1,
           TRUE                        ~ crown_diameter_for_lai_lumped_spp
         ),
         # Ongoing 
         crown_area_total = pi*crown_diameter_for_lai*(height_for_lai + crown_diameter_for_lai)/2,
         crown_area_total_lumped_spp = pi*crown_diameter_for_lai_lumped_spp*(height_for_lai_lumped_spp + crown_diameter_for_lai_lumped_spp)/2,
         shading_coeff = 0.0617 * log(dbh_cm) + 0.615 + shading_coeff_mean,
         shading_coeff_imputed_dbh = 0.0617 * log(dbh_pred) + 0.615 + shading_coeff_mean,
         shading_coeff_lumped_spp = 0.0617 * log(dbh_pred) + 0.615) |> 
  mutate(shading_coeff = pmin(pmax(shading_coeff, 0), 1),
         shading_coeff_imputed_dbh = pmin(pmax(shading_coeff_imputed_dbh, 0), 1),
         shading_coeff_lumped_spp = pmin(pmax(shading_coeff_lumped_spp, 0), 1)) |> 
  mutate(crown_diameter_eff = ifelse(height_width_ratio < 0.5, crown_height_m*2, crown_diameter_m),
         log_leaf_area = -4.3309 + 0.2942*height_for_lai + 0.7312*crown_diameter_for_lai + 5.7217*shading_coeff - 0.0148*crown_area_total, 
         leaf_area = exp(log_leaf_area)*ifelse(height_width_ratio>2, height_width_ratio/2, 1)*ifelse(height_width_ratio<0.5, 0.5/height_width_ratio, 1)*(crown_diameter_eff / crown_diameter_for_lai)^2,
         LAI = leaf_area/ground_area_projected_m_sq,
         log_leaf_area_imputed_dbh = -4.3309 + 0.2942*height_for_lai + 0.7312*crown_diameter_for_lai + 5.7217*shading_coeff_imputed_dbh - 0.0148*crown_area_total, 
         leaf_area_imputed_dbh = exp(log_leaf_area_imputed_dbh)*ifelse(height_width_ratio>2, height_width_ratio/2, 1)*ifelse(height_width_ratio<0.5, 0.5/height_width_ratio, 1)*(crown_diameter_eff/crown_diameter_for_lai)^2,
         LAI_imputed_dbh = leaf_area_imputed_dbh/ground_area_projected_m_sq,
         log_leaf_area_lumped_spp = -4.3309 + 0.2942*height_for_lai_lumped_spp + 0.7312*crown_diameter_for_lai_lumped_spp + 5.7217*shading_coeff_lumped_spp - 0.0148*crown_area_total_lumped_spp, 
         leaf_area_lumped_spp = exp(log_leaf_area_lumped_spp)*ifelse(height_width_ratio_lumped_spp>2, height_width_ratio_lumped_spp/2, 1)*ifelse(height_width_ratio_lumped_spp<0.5, 0.5/height_width_ratio_lumped_spp, 1)*(crown_diameter_eff/crown_diameter_for_lai_lumped_spp)^2,
         LAI_lumped_spp = leaf_area_lumped_spp/ground_area_projected_m_sq,
         log_leaf_area_timilsina = -3.21 + 0.16*height_for_lai + 0.43*crown_diameter_for_lai + 4.80*shading_coeff - 0.004*crown_area_total, 
         leaf_area_timilsina = exp(log_leaf_area_timilsina)*ifelse(height_width_ratio>2, height_width_ratio/2, 1)*ifelse(height_width_ratio<0.5, 0.5/height_width_ratio, 1)*(crown_diameter_eff / crown_diameter_for_lai)^2,
         LAI_timilsina = leaf_area_timilsina/ground_area_projected_m_sq)

# Cut off unrealistic LAI - force to range from 0 to 15
tree_data <- tree_data |> 
  mutate(
    LAI = case_when(
      is.na(LAI) ~ NA,
      LAI < 0    ~ 0,
      LAI > 15   ~ 15,
      TRUE       ~ LAI
    ),
    LAI_imputed_dbh = case_when(
      is.na(LAI_imputed_dbh) ~ NA,
      LAI_imputed_dbh < 0    ~ 0,
      LAI_imputed_dbh > 15   ~ 15,
      TRUE                   ~ LAI_imputed_dbh
    ),
    LAI_lumped_spp = case_when(
      is.na(LAI_lumped_spp) ~ NA,
      LAI_lumped_spp < 0    ~ 0,
      LAI_lumped_spp > 15   ~ 15,
      TRUE                   ~ LAI_lumped_spp
  ))

# Add in the LAI for trees in a dense forest setting, and for trees in medium density
tree_data <- tree_data |> 
  mutate(LAI_forest = -log((1-shading_coeff)/0.65),
         LAI_medium_density = (LAI_forest + LAI)/2,
         LAI_forest_imputed_dbh = -log((1-shading_coeff_imputed_dbh)/0.65),
         LAI_medium_density_imputed_dbh = (LAI_forest_imputed_dbh + LAI_imputed_dbh)/2)

# Compare LAI estimates with imputed vs. measured DBH
imputed_dbh_impact_model <- lm(data=tree_data, 
                               LAI ~ LAI_imputed_dbh)
summary(imputed_dbh_impact_model)
rmse(tree_data$LAI, tree_data$LAI_imputed_dbh)
imputed_dbh_impact_plot <- ggplot() + 
  geom_density_2d_filled(data=tree_data,
                         aes(x=LAI, y=LAI_imputed_dbh), bins=12) +
  geom_abline(intercept=0, slope=1, linetype="dashed", col="black") + 
  geom_abline(intercept=summary(imputed_dbh_impact_model)$coefficients[1,1], 
              slope=summary(imputed_dbh_impact_model)$coefficients[2,1], 
              linetype="solid", col="red") + 
  annotate("text", x=0.5, y=6, parse=TRUE, hjust=0,
           label = paste("R^2 == ", 
                         round(summary(imputed_dbh_impact_model)$adj.r.squared, 2))) + 
  annotate("text", x=0.5, y=5.5, hjust=0,
           label = paste("RMSE = ", 
                         round(rmse(tree_data$LAI, tree_data$LAI_imputed_dbh), 2))) + 
  scale_x_continuous(limits=c(0, 7), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) + 
  scale_fill_manual(values=colorRampPalette(c("white", "blue"))(12)) + 
  theme_bw() +
  xlab("LAI (from Measured DBH)") +
  ylab("LAI (from Modeled DBH)") + 
  theme(legend.position = "none")
imputed_dbh_impact_plot
ggsave("output_plots/imputed_dbh_impact_plot.png", imputed_dbh_impact_plot, 
       height=4, width=4)


# Compare LAI estimates with vs. without species information on shading 
species_impact_model <- lm(data=tree_data, 
                               LAI ~ LAI_lumped_spp)
summary(species_impact_model)
rmse(tree_data$LAI, tree_data$LAI_lumped_spp)
species_impact_plot <- ggplot() + 
  geom_density_2d_filled(data=tree_data,
                         aes(x=LAI, y=LAI_lumped_spp), bins=12) +
  geom_abline(intercept=0, slope=1, linetype="dashed", col="black") + 
  geom_abline(intercept=summary(species_impact_model)$coefficients[1,1], 
              slope=summary(species_impact_model)$coefficients[2,1], 
              linetype="solid", col="red") + 
  annotate("text", x=0.5, y=6, parse=TRUE, hjust=0,
           label = paste("R^2 == ", 
                         round(summary(species_impact_model)$adj.r.squared, 2))) + 
  annotate("text", x=0.5, y=5.5, hjust=0,
           label = paste("RMSE = ", 
                         round(rmse(tree_data$LAI, tree_data$LAI_lumped_spp), 2))) + 
  scale_x_continuous(limits=c(0, 7), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) + 
  scale_fill_manual(values=colorRampPalette(c("white", "blue"))(12)) + 
  theme_bw() +
  xlab("LAI (With Species Information)") +
  ylab("LAI (No Species-specific Coefficients)") + 
  theme(legend.position = "none")
species_impact_plot
ggsave("output_plots/species_impact_plot.png", species_impact_plot, 
       height=4, width=4)


# Visualize some of the data spread for the larger LAI dataset with estimated DBH
LAI_dbh_plot <- ggplot(tree_data) + 
  geom_density_2d_filled(aes(x=dbh_cm, y=LAI), bins=12) + 
  geom_smooth(aes(x=dbh_cm, y=LAI), col="red",
              method = "gam", formula = y ~ s(x, k = 20)) + 
  scale_x_continuous(limits=c(0, 100), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) + 
  scale_fill_manual(values=colorRampPalette(c("white", "blue"))(12)) + 
  theme_bw() + 
  xlab("Diameter at Breast Height (cm)") + 
  ylab("LAI") +
  theme(legend.position = "none")
LAI_dbh_plot
ggsave("output_plots/LAI_dbh_plot.png", LAI_dbh_plot, 
       height=4, width=4)

LAI_dbh_genera_plot <- ggplot(tree_data |> 
                                filter(genus_merged %in% classified_genera) |> 
                                group_by(genus_merged) |> 
                                mutate(dbh_q05 = quantile(dbh_cm, 0.05, na.rm=TRUE),
                                       dbh_q95 = quantile(dbh_cm, 0.95, na.rm=TRUE)) |> 
                                filter(dbh_cm > dbh_q05, dbh_cm < dbh_q95)) + 
  geom_density_2d_filled(aes(x=dbh_cm, y=LAI), bins=12) + 
  geom_smooth(aes(x=dbh_cm, y=LAI, group=genus_merged, col=genus_merged),
              method = "gam", formula = y ~ s(x, k = 8), se=FALSE) + 
  scale_x_continuous(limits=c(0, 100), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) + 
  scale_fill_manual(values=colorRampPalette(c("white", "black"))(12)) + 
  theme_bw() + 
  xlab("Diameter at Breast Height (cm)") + 
  ylab("LAI") + 
  guides(fill = "none", 
         col = guide_legend(ncol=2)) + 
  #theme(legend.position = "inside",
  #      legend.position.inside = c(1, 0),
  #      legend.justification = c(1, 0)) +
  labs(col="Genus")
LAI_dbh_genera_plot
ggsave("output_plots/LAI_dbh_genera_plot.png", LAI_dbh_genera_plot, 
       height=5.5, width=8)



LAI_crown_diameter_plot <- ggplot(tree_data) + 
  geom_density_2d_filled(aes(x=crown_diameter_m, y=LAI), bins=12) + 
  geom_smooth(aes(x=crown_diameter_m, y=LAI), col="red",
              method = "gam", formula = y ~ s(x, k = 20)) + 
  scale_x_continuous(limits=c(0, 20), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) + 
  scale_fill_manual(values=colorRampPalette(c("white", "blue"))(12)) + 
  theme_bw() + 
  xlab("Crown Diameter (m)") + 
  ylab("LAI") +
  theme(legend.position = "none")
LAI_crown_diameter_plot
ggsave("output_plots/LAI_crown_diameter_plot.png", LAI_crown_diameter_plot, 
       height=4, width=4)

LAI_crown_diameter_genera_plot <- ggplot(tree_data |> 
                                           filter(genus_merged %in% classified_genera) |> 
                                           group_by(genus_merged) |> 
                                           mutate(crown_diameter_q05 = quantile(crown_diameter_m, 0.05, na.rm=TRUE),
                                                  crown_diameter_q95 = quantile(crown_diameter_m, 0.95, na.rm=TRUE)) |> 
                                           filter(crown_diameter_m > crown_diameter_q05, 
                                                  crown_diameter_m < crown_diameter_q95)) + 
  geom_density_2d_filled(aes(x=crown_diameter_m, y=LAI), bins=12) + 
  geom_smooth(aes(x=crown_diameter_m, y=LAI, group=genus_merged, col=genus_merged),
              method = "gam", formula = y ~ s(x, k = 2), se=FALSE) + 
  scale_x_continuous(limits=c(0, 40), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) + 
  scale_fill_manual(values=colorRampPalette(c("white", "black"))(12)) + 
  theme_bw() + 
  xlab("Crown Diameter (m)") + 
  ylab("LAI") + 
  guides(fill = "none", 
         col = guide_legend(ncol=2)) + 
  #theme(legend.position = "inside",
  #      legend.position.inside = c(1, 0),
  #      legend.justification = c(1, 0)) +
  labs(col="Genus")
LAI_crown_diameter_genera_plot
ggsave("output_plots/LAI_crown_diameter_genera_plot.png", LAI_crown_diameter_genera_plot, 
       height=5.5, width=8)


LAI_height_plot <- ggplot(tree_data) + 
  geom_density_2d_filled(aes(x=max_tree_height_m, y=LAI), bins=12) + 
  geom_smooth(aes(x=max_tree_height_m, y=LAI), col="red",
              method = "gam", formula = y ~ s(x, k = 20)) + 
  scale_x_continuous(limits=c(0, 20), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) + 
  scale_fill_manual(values=colorRampPalette(c("white", "blue"))(12)) + 
  theme_bw() + 
  xlab("Tree Height (m)") + 
  ylab("LAI") +
  theme(legend.position = "none")
LAI_height_plot
ggsave("output_plots/LAI_height_plot.png", LAI_height_plot, 
       height=4, width=4)

LAI_height_genera_plot <- ggplot(tree_data |> 
                                   filter(genus_merged %in% classified_genera) |> 
                                   group_by(genus_merged) |> 
                                   mutate(max_tree_height_q05 = quantile(max_tree_height_m, 0.05, na.rm=TRUE),
                                          max_tree_height_q95 = quantile(max_tree_height_m, 0.95, na.rm=TRUE)) |> 
                                   filter(max_tree_height_m > max_tree_height_q05, 
                                          max_tree_height_m < max_tree_height_q95)) + 
  geom_density_2d_filled(aes(x=max_tree_height_m, y=LAI), bins=12) + 
  geom_smooth(aes(x=max_tree_height_m, y=LAI, group=genus_merged, col=genus_merged),
              method = "gam", formula = y ~ s(x, k = 2), se=FALSE) + 
  scale_x_continuous(limits=c(0, 40), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) + 
  scale_fill_manual(values=colorRampPalette(c("white", "black"))(12)) + 
  theme_bw() + 
  xlab("Tree Height (m)") + 
  ylab("LAI") + 
  guides(fill = "none", 
         col = guide_legend(ncol=2)) + 
  #theme(legend.position = "inside",
  #      legend.position.inside = c(1, 0),
  #     legend.justification = c(1, 0)) +
  labs(col="Genus")
LAI_height_genera_plot
ggsave("output_plots/LAI_height_genera_plot.png", LAI_height_genera_plot, 
       height=5.5, width=8)


LAI_dbh_errorbar_plots <- ggplot(tree_data |> 
                                   filter(genus_merged %in% classified_genera) |> 
                                   group_by(genus_merged) |> 
                                   summarize(dbh_q25 = quantile(dbh_pred, 0.25, na.rm=TRUE),
                                             dbh_q50 = quantile(dbh_pred, 0.50, na.rm=TRUE),
                                             dbh_q75 = quantile(dbh_pred, 0.75, na.rm=TRUE),
                                             LAI_q25 = quantile(LAI_imputed_dbh, 0.25, na.rm=TRUE),
                                             LAI_q50 = quantile(LAI_imputed_dbh, 0.50, na.rm=TRUE),
                                             LAI_q75 = quantile(LAI_imputed_dbh, 0.75, na.rm=TRUE))) + 
  geom_point(aes(x=dbh_q50, y=LAI_q50)) + 
  geom_errorbar(aes(x=dbh_q50, y=LAI_q50, xmin=dbh_q25, xmax=dbh_q75, col=genus_merged)) + 
  geom_errorbar(aes(x=dbh_q50, y=LAI_q50, ymin=LAI_q25, ymax=LAI_q75, col=genus_merged)) + 
  theme_bw() + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.99, 0.01),
        legend.justification = c(1, 0)) +
  labs(x="Diameter at Breast Height (cm)",
       y="Leaf Area Index", 
       col="Genus") + 
  scale_x_continuous(limits=c(0,100), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0,7), expand=c(0,0)) +
  guides(col = guide_legend(ncol=2))
LAI_dbh_errorbar_plots
ggsave("output_plots/LAI_dbh_errorbar_plots.png", LAI_dbh_errorbar_plots, 
       height=5.5, width=5.5)

LAI_genus_plot <- ggplot(tree_data |> 
                           drop_na(genus_merged) |> 
                           group_by(genus_merged) |> 
                           mutate(count = n()) |> 
                           filter(count > 5000,
                                  ! (genus_merged %in% c("Ailanthus", "Styphnolobium"))) |> 
                           arrange(genus_merged)) + 
  geom_boxplot(aes(x=genus_merged, y=LAI, group=genus_merged), alpha=0.01) + 
  scale_y_continuous(limits=c(0, 10)) + 
  theme_bw() + 
  xlab("Genus") + 
  ylab("LAI") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
LAI_genus_plot
ggsave("output_plots/LAI_genus_plot.png", LAI_genus_plot, 
       height=4, width=8)


write_csv(tree_data, "tree_lai_estimates.csv")

# Add in the spatial information 
geometry_data <- st_read(
  dsn = "F:/NYC/Tree_Data/crowns/miller_scientific_data/nyc_class_tree_genus_polygons_v2.gpkg", 
  query = "SELECT Poly_ID, geom FROM nyc_class_tree_genus_polygons_v2"
)
tree_points <- st_centroid(geometry_data)
tree_points <- tree_points |> 
  mutate(poly_id = Poly_ID) |> 
  dplyr::select(-Poly_ID) |> 
  left_join(tree_data)
st_write(tree_points %>% st_as_sf(),
         "tree_lai_estimates.gpkg", append=FALSE)

tree_polygons <- geometry_data |> 
  mutate(poly_id = Poly_ID) |> 
  dplyr::select(-Poly_ID) |> 
  left_join(tree_data)
st_write(tree_polygons %>% st_as_sf(),
         "tree_polygon_lai_estimates.gpkg", append=FALSE)



# Also try comparing vs. the UFIA urban forestry estimates:
ufia_top_genera <- read_csv("top_25_genera_nyc_ufia_2024.csv") |> 
  janitor::clean_names()
names(ufia_top_genera) <- c("genus_merged", paste0("ufia_", names(ufia_top_genera)))

genus_summary_vs_ufia <- tree_data |> 
  filter(genus_merged %in% classified_genera) |> 
  group_by(genus_merged) |> 
  #drop_na(dbh_cm) |> 
  summarize(dbh_q25 = quantile(dbh_pred, 0.25, na.rm=TRUE),
            dbh_q50 = quantile(dbh_pred, 0.50, na.rm=TRUE),
            dbh_q75 = quantile(dbh_pred, 0.75, na.rm=TRUE),
            LAI_q25 = quantile(LAI_imputed_dbh, 0.25, na.rm=TRUE),
            LAI_q50 = quantile(LAI_imputed_dbh, 0.50, na.rm=TRUE),
            LAI_q75 = quantile(LAI_imputed_dbh, 0.75, na.rm=TRUE),
            leaf_area_total = sum(leaf_area, na.rm=TRUE),
            num_trees = n()) |> 
  left_join(ufia_top_genera)

ggplot(genus_summary_vs_ufia) + 
  geom_point(aes(x=num_trees, ufia_number_of_trees))
summary(lm(data=genus_summary_vs_ufia, 
           ufia_number_of_trees ~ num_trees))
ggplot(genus_summary_vs_ufia) + 
  geom_point(aes(x=leaf_area_total, ufia_leaf_area))
summary(lm(data=genus_summary_vs_ufia, 
           leaf_area_total ~ ufia_leaf_area))


