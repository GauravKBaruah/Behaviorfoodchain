# Display session information 
# sessionInfo()
# R version 4.4.2 (2024-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_Germany.utf8  LC_CTYPE=English_Germany.utf8   
# [3] LC_MONETARY=English_Germany.utf8 LC_NUMERIC=C                    
# [5] LC_TIME=English_Germany.utf8    
# 
# time zone: Europe/Berlin
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.4   ggplot2_3.5.1
# 
# loaded via a namespace (and not attached):
#   [1] labeling_0.4.3    R6_2.5.1          tidyselect_1.2.1  farver_2.1.2     
# [5] magrittr_2.0.3    gtable_0.3.6      glue_1.8.0        tibble_3.2.1     
# [9] pkgconfig_2.0.3   generics_0.1.3    lifecycle_1.0.4   cli_3.6.3        
# [13] viridisLite_0.4.2 scales_1.3.0      grid_4.4.2        vctrs_0.6.5      
# [17] withr_3.0.2       compiler_4.4.2    rstudioapi_0.17.1 tools_4.4.2      
# [21] munsell_0.5.1     pillar_1.10.1     colorspace_2.1-1  crayon_1.5.3     
# [25] rlang_1.1.4      
# 


#######################################################################
# tri-trophic food-chain persistence
#### Figure 3

#
# 1. Setup
#
rm(list=ls())
sessionInfo()

library(ggplot2)
library(survival)
library(coxme)
library(dplyr)
library(AICcmodavg)
library(ggsurvfit)

#
# 2. Load data
#
c1 <- read.csv("Data_aphid_ladybird.csv")
str(c1)

# Subset and calculate survival flags
cox_dat <- subset(c1, Sampling_day > -1) %>%
  mutate(
    aphid_survival = ifelse(Aphids > 0, 1, 0),
    ladybug_survival = ifelse(Total_ladybirds > 0, 1, 0),
    Plant_survival = Plant_survival, #  Plant_survival is coded 1=alive, 0=dead
    Extinction = ifelse(aphid_survival == 1 & ladybug_survival == 1 & Plant_survival == 1, 0, 1)
  )

# Prepare dataset
cox_df <- cox_dat %>%
  select(Cage_ID, Sampling_day, Position_of_cage, Aphid_clone, Climate_chamber, 
         Treatment, No_of_leaves, Plant_Age, Extinction)

# Convert to factors
cox_df <- cox_df %>%
  mutate(across(c(Climate_chamber, Aphid_clone, Cage_ID, Position_of_cage, Treatment), as.factor))

#
# 3. Full food chain persistence model
#
# Full mixed model
mixed_full_model <- coxme(Surv(Sampling_day, Extinction) ~ Treatment + Plant_Age + No_of_leaves +
                            (1|Climate_chamber/Position_of_cage) + (1|Aphid_clone/Cage_ID),
                          data = cox_df)
summary(mixed_full_model)

# Reduced models
reduced1 <- coxme(Surv(Sampling_day, Extinction) ~ Treatment + Plant_Age +
                    (1|Climate_chamber/Position_of_cage) + (1|Aphid_clone/Cage_ID),
                  data = cox_df)

reduced2 <- coxme(Surv(Sampling_day, Extinction) ~ Treatment +
                    (1|Climate_chamber/Position_of_cage) + (1|Aphid_clone/Cage_ID),
                  data = cox_df)

reduced3 <- coxme(Surv(Sampling_day, Extinction) ~ Treatment + (1|Aphid_clone/Cage_ID), 
                  data = cox_df)

# Model selection
aictab(list(mixed_full_model, reduced1, reduced2, reduced3), 
       c("full", "reduced1", "reduced2", "reduced3"))

# Let's assume reduced3 is best
best_mixed_model <- reduced3
summary(best_mixed_model)

# Extract fixed effects
fixed_effects_df <- as.data.frame(summary(best_mixed_model)$coefficients)
colnames(fixed_effects_df) <- c("estimate", "exp_estimate", "se", "z_value", "p_value")
rownames(fixed_effects_df) <- gsub("Treatment", "", rownames(fixed_effects_df))

#
# 4. Plot full food chain persistence
#

treatment_colors <- c("Dropper" = "#B2182B", "Mix" = "#F4A582", "Non-dropper" = "#2166AC")

(p_persistence <- survfit2(Surv(Sampling_day, Extinction) ~ Treatment, data = cox_df) %>%
    ggsurvfit(linetype_aes = TRUE) +
    labs(x = "Days", y = "Persistence probability of full food chain") +
    add_confidence_interval() +
    scale_color_manual(values = treatment_colors) +
    scale_fill_manual(values = treatment_colors) +
    scale_linetype_manual(values = c("1342", "1342", "1342")) +
    theme_bw()+
    theme(legend.position  = c(0.2,0.1),
          legend.background = element_rect(#fill="lightblue",
            size=0.5, linetype="solid", 
            colour ="darkblue")))

# Effect size plot
(p_effects <- fixed_effects_df %>% 
    mutate(Treatment = rownames(.)) %>%
    ggplot(aes(x=Treatment, y = exp_estimate, col=Treatment)) +
    geom_point(size=4) +
    geom_linerange(aes(ymin = exp_estimate - se, ymax = exp_estimate + se), size=1) +
    labs(x="Treatment", y = "Extinction rate estimate (hazard ratio)") +
    theme_bw() +
    geom_hline(yintercept = 1, linetype="dashed", col="#B2182B") +
    scale_color_manual(values = treatment_colors) +
    theme(legend.position = "none"))


# Create species-specific extinction datasets
cox_aphid <- cox_dat %>% 
  mutate(Species = "Aphid", Extinction = ifelse(aphid_survival == 1, 0, 1))

cox_ladybird <- cox_dat %>% 
  mutate(Species = "Ladybird", Extinction = ifelse(ladybug_survival == 1, 0, 1))

cox_plant <- cox_dat %>% 
  mutate(Species = "Plant", Extinction = ifelse(Plant_survival == 1, 0, 1))


# Combine species-specific datasets
cox_df_species <- bind_rows(
  cox_aphid,
  cox_ladybird,
  cox_plant
)

# Make sure Species is factor
cox_df_species$Species <- factor(cox_df_species$Species, levels = c("Plant", "Aphid", "Ladybird"))

surv_species <- survfit2(Surv(Sampling_day, Extinction) ~ Treatment + Species, 
                         data = cox_df_species)

# Define linetypes for species
species_linetypes <- c("Plant" = "twodash", 
                       "Aphid" = "solid", 
                       "Ladybird" = "dashed")

# Plot
(p_species_combined <- surv_species %>%
    ggsurvfit(linetype_aes = TRUE) +
    labs(x = "Days", y = "Persistence probability (species-specific)") +
    add_confidence_interval() +
    scale_color_manual(values = c("#B2182B","#B2182B","#B2182B",
                                  "#F4A582","#F4A582","#F4A582",
                                  "#2166AC","#2166AC","#2166AC"))+
    scale_fill_manual(values = c("#B2182B","#B2182B","#B2182B",
                                 "#F4A582","#F4A582","#F4A582",
                                 "#2166AC","#2166AC","#2166AC"))+
    scale_linetype_manual(values = c("dotted", "longdash", "solid", 
                                     "dotted","longdash", "solid", 
                                     "dotted","longdash","solid" ))+
    theme_bw() +
    theme(legend.position  = c(0.25,0.25),
          legend.background = element_rect(#fill="lightblue",
            size=0.5, linetype="solid", 
            colour ="darkblue")))

# Show plot
p_species_combined


#
# 6. save final combined figures
#

# png(file = "Figure 3.png",width = 12,height = 6,units = 'in', res = 600)
ggpubr::ggarrange(p_persistence , p_effects,
                  p_species_combined, nrow = 1,ncol=3, labels=c("A","B","C"),
                  widths = c(1, 0.8, 1.3)  # Adjust these to  desired relative proportions
)
# dev.off()
