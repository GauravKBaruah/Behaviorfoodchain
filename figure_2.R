# Behavioral variation affects persistence of an experimental food-chain
# Authors: Pragya Singh, Gaurav Baruah,Caroline Müller

# Clear the R environment
rm(list = ls())  

#######################################################################
# Time-series of predator, prey abundance and plant survival
### Figure 2  

# Read the CSV data file into the object 'c'
c <- read.csv("Data_aphid_ladybird.csv")

# Display column names and structure of the dataset
colnames(c)
str(c)

# Check how many records exist for each cage
table(c$Cage_ID)

# Quick  plot to visualize data collection timeline per cage
ggplot(c, aes(x = Cage_ID, y = Sampling_day)) + 
  geom_point()

# data before predator added
sub.data <- subset(c, c$Sampling_day > -1)
table(sub.data$Cage_ID)
# Plot sampling days per cage 
ggplot(sub.data, aes(x = Cage_ID, y = Sampling_day)) + 
  geom_point()
# Subset data from pre-predator phase 
sub.data1 <- subset(c, c$Sampling_day < 0)
# Check data of pre-treatment 
table(sub.data1$Sampling_day)
table(sub.data1$Sampling_day, sub.data1$Cage_ID)
table(sub.data1$Sampling_day, sub.data1$Treatment)

library(dplyr)
# Load the ggplot2 library for plotting
library(ggplot2)

# Extract extinction days per cage
plant_extinction <- c %>%
  filter(Plant_survival == 0) %>%
  group_by(Cage_ID, Treatment) %>%
  summarise(Extinction_day = min(Sampling_day), .groups = "drop")

# Join extinction info back with full data to get aphid abundance at extinction day
extinction_points <- plant_extinction %>%
  left_join(c, by = c("Cage_ID", "Treatment")) %>%
  filter(Sampling_day == Extinction_day) %>%
  select(Cage_ID, Treatment, Extinction_day, Aphids)

# png(file = "Figure 2.png",width = 12 ,height = 4,units = 'in', res = 600)
ggplot() +
  # Ladybird abundance lines (solid)
  geom_line(data = c,
            aes(x = Sampling_day,
                y = log(Total_ladybirds + 1),
                group = Cage_ID,
                color = as.factor(Cage_ID),
                linetype = "Ladybirds"),
            size = 1.75) +
  
  # Aphid abundance lines (dashed)
  geom_line(data = c,
            aes(x = Sampling_day,
                y = log(Aphids + 1),
                group = Cage_ID,
                color = as.factor(Cage_ID),
                linetype = "Aphids"),
            size = 1) +
  
  # Predator introduction vertical line
  geom_vline(xintercept = 0, color = "firebrick", linetype = "dashed", size = 1.25) +
  
  # Plant extinction markers: star at aphid abundance, same cage color
  geom_text(data = extinction_points,
            aes(x = Extinction_day, y = log(Aphids + 1), fill = "black"), 
            # as.factor(Cage_ID)),
            label = "\u2605",  # Unicode for solid star ★
            size = 8) +
  # Axis label
  ylab("Abundance of aphids and ladybirds\n(log-transformed)") +
  xlab("Sampling day") +
  
  # Facet by Treatment
  facet_wrap(~Treatment) +
  
  # Legend for species (only for linetypes)
  scale_linetype_manual(name = "Species",
                        values = c("Ladybirds" = "solid", "Aphids" = "longdash")) +
  
  # Keep same color scale for cage IDs (reused automatically)
  scale_color_viridis_d(guide = "none") +
  
  # General theme
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))
