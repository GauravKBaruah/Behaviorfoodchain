# README for the paper titled "Behavioral variation affects persistence of an experimental food-chain."

*Authors: Pragya Singh, Gaurav Baruah, Caroline Müller*

*This R scripts and data follows GNU General public license v3.0.
# R scripts for reproducing the figures

*1.* `figure_2.R` : R script to reproduce figure 2.
*2.* `figure_3.R` : R script to reproduce figure 3.
*3.* `figure_4.R` : R script to reproduce figure 4.
*4.* `figure_5.R` : R script to reproduce figure 5. This is Empirical Dynamical Analysis (EDM). We use the for this
script Redm package of which the R package version has to be 0.6.9.


# Rawdata for "Data_aphid_ladybird.csv" file

Treatment : Treatment of replicate cage
Cage_ID : Cage Id of replicate cage
Sampling_day : Day on which sampling was done
Date_setup_cage : date on which cage was setup
Aphid_clone : aphid clone ID
plant_height : Height of plant (cm)
Plant_Age : in weeks
Planted_on : date of planting
No_of_leaves : No. of leaves on the plants
Position_of_cage : Position of cage in climate chamber
Climate_chamber : Climate chamber used
Aphids : Aphids abundance (adults + nymphs)
Ladybird_adults : Adult ladybird abudnance 
Ladybird_larvae : larval ladybird abundance
Total_ladybirds : Total ladybird abundance (adult + larvae)
Plant_survival : Plant survived (1) or dead (0)
Aphids_orig : Aphid abudnance without imputation
Total_ladybirds_orig: Total ladybird abundance (adult + larvae) without imputation
Plant_survival_orig: Plant survived (1) or dead (0) without imputation
Aphids_flag : Aphid abudnance Whether observed or imputed (interpolated/forward filled)
Ladybirds_flag : Ladybird abundance Whether observed or imputed (interpolated/forward filled)
Survival_flag : Plant abudance Whether observed or imputed (interpolated/forward filled)