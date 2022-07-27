?mite
# I recall the vegan package, since it is the one required by mite dataset
library(vegan)

### Import the datasets
# Since mite_env is a .txt, I use the read.table() function.
  # With sep = " ", I treat the missing values as empty cells, without the value
  # With header = T, I take the first row and put it as columns header
  # With na.strings = c("NA", ""), I make R recognize as NA all the NA and the missing values
mite_env <- read.table("data/mite_env.txt",
                       sep = " ", 
                       header = T, 
                       na.strings = c("NA", ""))
mite_env
# I then check the number of rows and columns of the dataset 
nrow(mite_env)
ncol(mite_env)
# mite_env has 70 observations (sampling sites, rows) and 5 environmental variables (columns)

# Since mite is a .csv, I use the read.csv() function.
mite <- read.csv("data/mite.csv")
mite
# I then check the number of rows and columns of the dataset 
nrow(mite)
ncol(mite)
# mite has 70 observations (sampling sites, rows) and 35 variables (species, columns)



### Check the datasets, their structure
?mite

str(mite_env)
# mite_env is a dataframe which contains information about the environmental variables measured for each species. There are: Substrate density in g/L (SubsDens), Water content of the substrate in g/L (WatrCont), Substrate type (Substrate), Shrub density (Shrub) and Microtopography (Topo)
# the variables SubsDens and WatrCont are numeric data types
# the variables Substrate, Shrub and Topo are factors, character data types

str(mite)
# mite is a dataframe composed by integer data type. It contains information on the absolute frequency of species (columns) for each sampling site (rows)

# I then have to transform the variables "Substrate", "Shrub" and "Topo" of mite_env into factor data type.
# "Shrub" is an ordered factor data type

unique(mite_env$Substrate)
mite_env$Substrate <- factor(mite_env$Substrate,
                             levels = c("Sphagn1", "Sphagn2", "Sphagn3", "Sphagn4", "Litter", "Barepeat", "Interface"))

unique(mite_env$Shrub)
mite_env$Shrub <- factor(mite_env$Shrub,
                         levels = c("None","Few","Many", NA),
                         ordered = T)

unique(mite_env$Topo)
mite_env$Topo <- factor(mite_env$Topo,
                        levels = c("Blanket", "Hummock"))



### Summary statistics of mite_env and removal of NAs
summary(mite_env)

# I notice NAs in "SubsDens" and in "Shrub". 
# I double check with the is.na() function
is.na(mite_env)
# I find the coordinates of the NA values
index.na <- which(is.na(mite_env), arr.ind = T)
index.na
# The NAs are in the row 2, column 1 ("SubsDens"), and row 48, column 4 ("Shrub")
# I assing to index.na the number of the rows where the NA values are.
index.na <- index.na[, 1] 

# I remove the rows containing them:
mite_env_wna <- na.omit(mite_env) # I assign to mite_env_wna the same dataframe of mite_env, but without the rows that contained NAs
mite_env_wna
# Now, mite_env_wna has 68 observations.

# I remove the same rows in the "mite" dataframe
mite_wna <- mite[-index.na, ]
mite_wna
str(mite_wna)
# Now, mite_wna has 68 observations.



### Graphics and univariate distribution of environmental variables
# For numerical data (SubsDens, WatrCont) I use histograms, while for character data (substrate, Shrub, Topo) I use barplots
# With "main" I put the title of the graph
# With "xlab" and "ylab" I put the title of the x axis and y axis, respectively
# with "xlim" and "ylim" I set the limits of the x and y axis

hist(mite_env_wna$SubsDens,
     main = "Substrate Density",
     xlab = "g/L",
     ylim = c(0, 30))

hist(mite_env_wna$WatrCont,
     main = "Water Content of the Substrate", 
     xlab = "g/L",
     xlim = c(0, 1000),
     ylim = c(0, 25))

barplot(table(mite_env_wna$Substrate),
        main = "Substrate type",
        xlab = "Type",
        ylab = "Frequency",
        ylim = c(0, 30))

barplot(table(mite_env_wna$Shrub),
        main = "Shrub density",
        xlab = "Levels",
        ylab = "Frequency")

barplot(table(mite_env_wna$Topo),
        main = "Microtopography",
        xlab = "Levels",
        ylab = "Frequency",
        ylim = c(0, 50))


# I export two graphs: the histogram for Water Content, and the barplot for the Substrate
# With the png() function I set the characteristics of the image I'll obtain. 
  # With "filename" I set the path and the name of the image.
  # With "width" and "height" I set the dimension of the image
  # With "res" I set the resolution of the image
# With the hist() function I obtain the histogram, as we did before (the same will be with the barplot() function)
# With the dev.off() function I close the exporting process

png(filename = "outputs/Histogram_WaterContent.png",
    width = 2000,
    height = 2000,
    res = 400)
hist(mite_env_wna$WatrCont,
     main = "Water Content of the Substrate", 
     xlab = "g/L",
     xlim = c(0, 1000),
     ylim = c(0, 25))
dev.off()

png(filename = "outputs/Barplot_Substrate.png",
    width = 3000,
    height = 3000,
    res = 400)
barplot(table(mite_env_wna$Substrate),
        main = "Substrate type",
        xlab = "Type",
        ylab = "Frequency",
        ylim = c(0, 30))
dev.off()



### Convert community matrix (mite_wna) into a presence/absence matrix

mite_pa <- decostand(mite_wna, method = "pa")
mite_pa



### Calculate species richness for each site

Spec_rich <- specnumber(mite_wna)
Spec_rich

# Add Species Richness into environmental variables dataframe (mite_env_wna) as a new column
mite_env_wna$SpecRich <- Spec_rich
mite_env_wna



### Plot species richness with respect to the numeric environmental variables (Substrate Density and Water Content)
# In the argument of the plot function, I put firstly the two variables I want to plot together. 
  # In "data", I put the dataframe where the environmental variables are  
  # "cex" is used to set the dimension of the circles that represent each point
plot(SpecRich ~ SubsDens,
     data = mite_env_wna,
     main = "Species Richness distribution with respect to Substrate Density",
     xlab = "Substrate Density",
     ylab = "Species Richness",
     ylim = c(0,30),
     cex = 1)

plot(SpecRich ~ WatrCont,
     data = mite_env_wna,
     main = "Species Richness distribution with respect to Water Content",
     xlab = "Water Content of the substrate",
     ylab = "Species Richness",
     ylim = c(0,30),
     cex = 1)


### Test the linear correlation between Species Richness and numerical environmental variables Substrate Density and Water Content

# I first test the correlation between Species Richness and Substrate Density. I use the correlation test.
# Null hypothesis (H0) =  correlation = 0 -> there is no correlation between Species Richness and Substrate Density
# Alternative hypothesis = correlation is not equal to 0 (two-sided, or positive or negative correlation)
cor.test(mite_env_wna$SpecRich, 
         mite_env_wna$SubsDens) 
# p-value = 0.6034 -> p-value is greater than 0.1, so it is not significant -> there is no evidence against H0 -> there is no correlation between the two variables
# cor = -0.06411822 -> there is no correlation between Substrate Density and Species Richness

# I then test the correlation between Species Richness and Water Content. I use the correlation test.
# Null hypothesis (H0) =  correlation = 0 -> there is no correlation between Species Richness and Water Content
# Alternative hypothesis = correlation is not equal to 0 (two-sided, or positive or negative correlation)
cor.test(mite_env_wna$SpecRich, 
         mite_env_wna$WatrCont)
# p-value = 2.157e-11 -> p-value is much smaller than 0.01, so it is very significant -> there is strong evidence against H0 -> there is correlation between the two variables
# cor = -0.7038521 -> there is significant negative correlation between Water Content and Species Richness



### Regression model for the correlated numeric environmental variable (Water Content)

# Since the correlation is linear, I use the linear regression model
lrm_WatrCont <- lm(SpecRich ~ WatrCont, 
                   data = mite_env_wna)
lrm_WatrCont
# Intercept = 24.23089: for water content = 0, there are 24 species (estimated).
# Slope = -0.02247: for each increase of one unit of Water Content, the number of species decreases by 0.02247

# To see graphically the negative correlation, I can obtain the regression line on the scatterplot with the abline() function
plot(SpecRich ~ WatrCont,
     data = mite_env_wna,
     main = "Species Richness distribution with respect to Water Content",
     xlab = "Water Content of the Substrate",
     ylab = "Species Richness",
     ylim = c(0,30),
     cex = 1)

abline(lrm_WatrCont$coefficients[1], 
       lrm_WatrCont$coefficients[2], 
       col = "blue")



### Inspection of the summary of the linear regression model

summary(lrm_WatrCont)
# Residuals (distance between the regression line and the actual observed value): they range from -7.9707 to 6.8692 
# Pr(>|t|): p-values for both estimated coefficients are very low, << 0.01. For the intercept is < 2e-16; for WatrCont is 2.16e-11 -> the model is very significant
  # this is shown also by the three stars (***)
# Adjusted R-squared = 0.4878 -> the regression line model explains 48.78% of the variability of our data -> it is moderate 


