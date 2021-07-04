library(SRS)
library(RColorBrewer)
options(scipen=999) #prevent scientific notation at plot ticks

#load dataset
tmo <- read.csv("tmo.tsv",sep = "\t", row.names = 1)

#subsample 100 samples
tmo_subset <- tmo[sample(1:ncol(tmo),100)]

###FIGURE 1
#launch shiny app
SRS.shiny.app(tmo_subset)
#screenshots taken from the browser app

###FIGURE 2
#sample 100 R colors (excluding shades of gray)
colors_100 = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],100)

#draw SRS curves
#figure 2a
SRScurve(tmo_subset, step = 10^3, max.sample.size = 10^5,
         xlab = "sequencing depth", ylab = "observed ASVs", 
         col = colors_100)
par(las = 1)
abline(v=10^4)
axis(1, at = seq(0, 10^5, 10^4), labels = seq(0, 10^5, 10^4))
axis(2, at = seq(0, 250, 25), labels = seq(0, 250, 25))

#figure 2b
SRScurve(tmo_subset, step = 10^3, max.sample.size = 10^5, 
         rarefy.comparison = T,rarefy.comparison.legend = T,
         xlab = "sequencing depth", ylab = "observed ASVs", 
         col = colors_100, lty = c("solid","dotdash"))
par(las = 1)
abline(v=10^4)
axis(1, at = seq(0, 10^5, 10^4), labels = seq(0, 10^5, 10^4))
axis(2, at = seq(0, 250, 25), labels = seq(0, 250, 25))