library(sswr)

library(sampling)
library(dplyr)

library(stratification)

library(ggplot2)
library(viridis)

Sys.setenv(JAVA_HOME = "C:/Program Files/Eclipse Adoptium/jdk-17.0.16.8-hotspot")
library(spcosa)

library(sp)

library(survey)

library(surveyplanning)

library(patchwork)

source("final_project_functions.R")
source("final_plot_functions.R")

BCI <- read.csv(
  "soildata.csv",
  sep = ";",           # columns separated by semicolon
  dec = ",",           # decimals written with commas
  stringsAsFactors = FALSE
)

n <- 60 # the total sample size 
H <- 5 # the total number of strata

# First, we should choose the type of stratification to apply to the entire population

# The whole population may be unobservable for the study variable of interest. However, we may observe other characteristics of the population,
# such as the values of related covariates or the coordinates of population units.

###### Cum-root-f stratification based on the covariate Fe

# Note that cum-root-f stratification is a deterministic algorithm
crfFe_strata <- strata.cumrootf(
  x = BCI$Fe, n = n, Ls = H, nclass = 500) # nclass is the number of bins of the histogram
bh_Fe <- crfFe_strata$bh # stratum boundaries

crfFe_strata
head(BCI)

# visualization of the resulting cum-root-f stratification based on the Fe

colors <- viridis(H, option = "D")

crfFe_labels <- c(
  paste0("<", round(bh_Fe[1],3)),
  paste0(round(bh_Fe[1:(length(bh_Fe)-1)],3), "–", round(bh_Fe[2:length(bh_Fe)],3)),
  paste0(">", round(bh_Fe[length(bh_Fe)],3))
)

crfFe_plot <- ggplot(BCI, aes(x = x, y = y, fill = as.factor(crfFe_strata$stratumID))) +
  geom_tile(color = "grey80", linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = colors,
                    name = "Fe",
                    labels = crfFe_labels) +
  coord_fixed(xlim = range(BCI$x),
              ylim = range(BCI$y),
              expand = FALSE) +
  labs(
    title = "Cum-root-f stratification of BCI, using Fe as a stratification variable.",
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()

crfFe_plot

#ggsave("crf_stratification_Fe.png", plot = crfFe_plot, width = 16, height = 7, dpi = 300, units ="cm")


###### Geographical stratification

BCI_geo <- BCI
gridded(BCI_geo) <- ~ x + y # here we are defining the grid structure based on the coordinates x and y

set.seed(271199) # We should set a seed because each time the k-means algorithm runs, it starts with new random starting points
geo_strata <- stratify(
  object = BCI_geo, nStrata = H, nTry = 100, equalArea = FALSE)
summary(as.factor(geo_strata@stratumId))

# visualization of the resulting geographical stratification

geo_plot <- ggplot(BCI, aes(x = x, y = y, fill = as.factor(geo_strata@stratumId))) +
  geom_tile(color = "grey80", linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = colors, name = "Fe") +
  coord_fixed(xlim = range(BCI$x),
              ylim = range(BCI$y),
              expand = FALSE) +
  labs(
    title = "Geographical stratification.",
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()

geo_plot

#ggsave("geo_stratification.png", plot = geo_plot, width = 16, height = 7, dpi = 300, units ="cm")

combined_plot <- crfFe_plot / geo_plot +
  plot_layout(ncol = 1, heights = c(1, 1))  # tři řádky

ggsave("combined_stratifications.png", combined_plot, width = 8, height = 10, dpi = 300)

######

head(BCI)
SIM_SEQ <- seq(1,1000)
SEQ_n <- seq(30, 250, by = 10)

###### Sample size vs Standard error

###

nlist_H3 <- RUN_PARALLEL_SIMULATIONS(BCI, SEQ_n, 3, SIM_SEQ)
plot_SE_H3 <- plot_SE_boxplots(nlist_H3$srswor, nlist_H3$geo, nlist_H3$crfFe, SEQ_n, 3)

###

nlist_H5 <- RUN_PARALLEL_SIMULATIONS(BCI, SEQ_n, 5, SIM_SEQ)
plot_SE_H5 <- plot_SE_boxplots(nlist_H5$srswor, nlist_H5$geo, nlist_H5$crfFe, SEQ_n, 5)

plot_SE_H_3_5 <- plot_SE_H3 / plot_SE_H5

# Save the combined plot
ggsave("SE_H_3_5.pdf", plot_SE_H_3_5, width = 12, height = 15)


###

nlist_H7 <- RUN_PARALLEL_SIMULATIONS(BCI, SEQ_n, 7, SIM_SEQ)
plot_SE_H7 <- plot_SE_boxplots(nlist_H7$srswor, nlist_H7$geo, nlist_H7$crfFe, SEQ_n, 7)

###

nlist_H10 <- RUN_PARALLEL_SIMULATIONS(BCI, SEQ_n, 10, SIM_SEQ)
plot_SE_H10 <- plot_SE_boxplots(nlist_H10$srswor, nlist_H10$geo, nlist_H10$crfFe, SEQ_n, 10)

plot_SE_H_7_10 <- plot_SE_H7 / plot_SE_H10

# Save the combined plot
ggsave("SE_H_7_10.pdf", plot_SE_H_7_10, width = 12, height = 15)


###### Sample size vs Design effect

plot_desef_H3 <- plot_desef_boxplots(nlist_H3$geo, nlist_H3$crfFe, SEQ_n, 3)
plot_desef_H5 <- plot_desef_boxplots(nlist_H5$geo, nlist_H5$crfFe, SEQ_n, 5)

plot_desef_H_3_5 <- plot_desef_H3 / plot_desef_H5

# Save the combined plot
ggsave("desef_H_3_5.pdf", plot_desef_H_3_5, width = 12, height = 15)



plot_desef_H7 <- plot_desef_boxplots(nlist_H7$geo, nlist_H7$crfFe, SEQ_n, 7)
plot_desef_H10 <- plot_desef_boxplots(nlist_H10$geo, nlist_H10$crfFe, SEQ_n, 10)

plot_desef_H_7_10 <- plot_desef_H7 / plot_desef_H10

# Save the combined plot
ggsave("desef_H_7_10.pdf", plot_desef_H_7_10, width = 12, height = 15)


###### Sample size vs mz
mz_true <- mean(BCI$Ca)

plot_mz_H3 <- plot_mz_boxplots(nlist_H7$srswor, nlist_H3$geo, nlist_H3$crfFe, SEQ_n, 3, mz_true)
plot_mz_H5 <- plot_mz_boxplots(nlist_H7$srswor, nlist_H5$geo, nlist_H5$crfFe, SEQ_n, 5, mz_true)

plot_mz_H_3_5 <- plot_mz_H3 / plot_mz_H5

# Save the combined plot
ggsave("mz_H_3_5.pdf", plot_mz_H_3_5, width = 12, height = 15)



plot_mz_H7 <- plot_mz_boxplots(nlist_H7$srswor, nlist_H7$geo, nlist_H7$crfFe, SEQ_n, 7, mz_true)
plot_mz_H10 <- plot_mz_boxplots(nlist_H7$srswor, nlist_H10$geo, nlist_H10$crfFe, SEQ_n, 10, mz_true)

plot_mz_H_7_10 <- plot_mz_H7 / plot_mz_H10

# Save the combined plot
ggsave("mz_H_7_10.pdf", plot_mz_H_7_10, width = 12, height = 15)
