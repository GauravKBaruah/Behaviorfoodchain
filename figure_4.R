
# > sessionInfo()
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
#   [1] ggsurvfit_1.1.0  AICcmodavg_2.3-4 dplyr_1.1.4      coxme_2.2-22     bdsmatrix_1.3-7 
# [6] survival_3.8-3   ggplot2_3.5.1   
# 
# loaded via a namespace (and not attached):
#   [1] generics_0.1.3    tidyr_1.3.1       rstatix_0.7.2     lattice_0.22-6   
# [5] magrittr_2.0.3    grid_4.4.2        Matrix_1.7-1      backports_1.5.0  
# [9] Formula_1.2-5     purrr_1.0.2       scales_1.3.0      abind_1.4-8      
# [13] cli_3.6.3         rlang_1.1.4       crayon_1.5.3      cowplot_1.1.3    
# [17] munsell_0.5.1     splines_4.4.2     withr_3.0.2       tools_4.4.2      
# [21] parallel_4.4.2    ggsignif_0.6.4    colorspace_2.1-1  ggpubr_0.6.0     
# [25] VGAM_1.1-13       unmarked_1.5.0    broom_1.0.7       vctrs_0.6.5      
# [29] R6_2.5.1          stats4_4.4.2      lifecycle_1.0.4   car_3.1-3        
# [33] MASS_7.3-64       pkgconfig_2.0.3   pillar_1.10.1     gtable_0.3.6     
# [37] glue_1.8.0        Rcpp_1.0.14       tibble_3.2.1      tidyselect_1.2.1 
# [41] rstudioapi_0.17.1 farver_2.1.2      xtable_1.8-4      nlme_3.1-166     
# [45] carData_3.0-5     labeling_0.4.3    compiler_4.4.2   
# > 


#######################################################################
#  state transitions over time in tri-trophic food-chain
#### Figure 4

rm(list=ls()) 
sessionInfo()

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(lubridate)
library(ggsurvfit)
library(survival)
library(coxme)
library(dplyr)
library(plyr)
library(tidyverse)
library(msm)
library(survminer)
library(broom)

state_dat<-read.csv(file="Data_aphid_ladybird.csv")
str(state_dat)

state_dat<-subset(state_dat,state_dat$Sampling_day>-1)

# Define survival flags
state_dat <- state_dat %>%
  mutate(
    aphid_survival = ifelse(Aphids > 0, 1, 0),
    ladybug_survival = ifelse((Total_ladybirds) > 0, 1, 0)
  )

# Summarize to one row per timepoint per cage
state_df <- state_dat %>% 
  group_by(Cage_ID, Sampling_day, Treatment, No_of_leaves, Plant_Age, Position_of_cage) %>%
  dplyr::summarise(across(c(aphid_survival, ladybug_survival, Plant_survival), max)) %>%
  ungroup() %>%
  mutate(
    Aphids = ifelse(aphid_survival > 0, "Aphid", "0"),
    LB = ifelse(ladybug_survival > 0, "LB", "0"),
    Plant = ifelse(Plant_survival > 0, "Plant", "0"),
    state = paste(Aphids, LB, Plant, sep = "-"),
    Cage_ID = as.character(Cage_ID)
  ) %>%
  filter(state != "NA-NA-NA", Sampling_day > 1) %>%
  mutate(Sampling_since = Sampling_day - 1)

# Create a unique state-to-number 
unique_states <- sort(unique(state_df$state))
unique_states
state_map <- setNames(seq_along(unique_states), unique_states)
state_map

# > state_map
# 0-LB-0     0-LB-Plant     Aphid-LB-0 Aphid-LB-Plant 
# 1              2              3              4 

# Map states to numeric codes
state_df <- state_df %>%
  mutate(state_num = state_map[state])

# Create numeric transition table
transition_table <- statetable.msm(state = state_num, subject = Cage_ID, data = state_df)
transition_table 


#ALLOWED TRANSITIONS

# > state_map
# 0-LB-0     0-LB-Plant     Aphid-LB-0 Aphid-LB-Plant 
# 1              2              3              4 


#ALLOWED TRANSITIONS
Q <- rbind(
  c( 0,  0,  0,  0),  # 1: 0-LB-0 doesnt go to anything absorbing state
  c( 1,  0,  0,  0),  # 2: 0-LB-Plant → can go to 0-LB-0
  c( 1,  0,  0,  0),  # 3: Aphid-LB-0 → 0-LB-0
  c( 1,  1,  1,  0)  # 4: Aphid-LB-Plant → can go to all 
)

rownames(Q) <- colnames(Q) <- c("0-LB-0", "0-LB-Plant", "Aphid-LB-0", "Aphid-LB-Plant")
Q



Q_init <- crudeinits.msm(state_num ~ Sampling_since, subject=Cage_ID,
                         data=state_df, qmatrix=Q)
Q_init

state_df$Treatment <- as.factor(state_df$Treatment)

state_df <- state_df %>% arrange(Cage_ID, Sampling_since)
state_df$next_state_num <- ave(state_df$state_num, state_df$Cage_ID, FUN = function(x) c(x[-1], NA))

# Find all transitions
transition_counts <- table(paste(state_df$state_num, state_df$next_state_num))
transition_counts
observed_transitions <- names(transition_counts[transition_counts > 0])
observed_transitions
observed_transitions <- na.omit(observed_transitions)
observed_transitions

# Fit model MARKOV MODEL 
mod <- msm(
  state_num ~ Sampling_since,
  subject = Cage_ID,
  data = state_df,
  qmatrix = Q,
  covariates = list(
    "4-3" = ~ Treatment, #CHECK THE SIGNIFICANCE OF TRANSITION FROM STATE 4 TO 3
    "4-2" = ~ Treatment#,   #CHECK THE SIGNIFICANCE OF TRANSITION FROM STATE 4 TO 2
    # "4-1" = ~ Treatment   #CHECK THE SIGNIFICANCE OF TRANSITION FROM STATE 4 TO 1
  ),
  control = list(trace = TRUE, REPORT = 1),
  gen.inits = TRUE  # let msm estimate good initial values
)

summary(mod)

#Hazard rates of the transitions tested
hzmsm<-hazard.msm(mod)
hzmsm

###plot
# Remap states into biological labels
state_labels <- c(
  "0-LB-0" = "Ladybird only",
  "0-LB-Plant" = "Plant–Ladybird",
  "Aphid-LB-0" = "Aphid–Ladybird",
  "Aphid-LB-Plant" = "Plant-Aphid-Ladybird"
)

state_df$state_label <- state_labels[state_df$state]

# Summarise empirical prevalence using tally()
empirical_prevalence <- state_df %>%
  group_by(Sampling_since, Treatment, state_label) %>%
  tally(name = "count") %>%
  group_by(Sampling_since, Treatment) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Plot
ggplot(empirical_prevalence, aes(x = Sampling_since, y = proportion, color = state_label)) +
  geom_line(size = 1.2) +
  facet_wrap(~ Treatment, nrow = 1) +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette = "Dark2", name = "Community state") +
  labs(x = "Sampling day", y = "Observed state proportion") +
  theme(legend.position = "bottom")

# geom_area plot below

library(tidyr)

# First, create a full grid of all combinations
all_combinations <- expand_grid(
  Sampling_since = unique(state_df$Sampling_since),
  Treatment = unique(state_df$Treatment),
  state_label = unique(state_labels)
)

# Merge this full grid with  empirical data
empirical_complete <- all_combinations %>%
  left_join(empirical_prevalence, by = c("Sampling_since", "Treatment", "state_label")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    proportion = ifelse(is.na(proportion), 0, proportion)
  )

empirical_complete$state_label <-factor(empirical_complete$state_label ,levels= 
                                          c("Ladybird only" , "Aphid–Ladybird",
                                            "Plant–Ladybird" , "Plant-Aphid-Ladybird" ))

# Now plot again using geom_area
# png(file = "Figure 4.png",width = 9 ,height = 4,units = 'in', res = 600)
ggplot(empirical_complete, aes(x = Sampling_since, y = proportion, fill = state_label)) +
  geom_area(alpha = 0.9, position = "stack") +
  facet_wrap(~ Treatment, nrow = 1) +
  theme_classic(base_size = 14) +
  guides(fill = guide_legend(reverse = TRUE))+
  scale_fill_brewer(palette = "BrBG", name = "Food-chain state") +
  labs(x = "Sampling day", y = "Observed state proportion") +
  theme(legend.position = "bottom")
# dev.off()
