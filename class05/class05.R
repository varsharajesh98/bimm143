#' ---
#' title: "Class 05 - Data Visualization and Graphs in R"
#' author: "Varsha Rajesh"
#' date: "Jan 23, 2020"
#' ---

plot(1:10, col = "blue", type = "o")

#Import/read input data file
weight_read <- read.table("bimm143_05_rstats/weight_chart.txt")
weight_read <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

#Plot weight_chart as scatterplot
weight_plot <- plot(weight_read$Age, weight_read$Weight)
#Plot weight_chart as line plot
weigh_plot <- plot(weight_read$Age, weight_read$Weight, typ = "o", pch = 15, cex = 1.5, lwd = 2,
              ylim = c(2,10), xlim = c(0,9), xlab = "Age(months)", ylab = "Weight(kg)", main = "Age Vs Weight of New-born Babies")

#Barplot for mouse genome features
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", sep = "\t", header = TRUE)
par(mar=c(5,12,9,4))
barplot(mouse$Count, horiz = TRUE, names.arg = mouse$Feature, las = 1, main="Mouse Genome Details", cex.names = 0.5)



#Barplot of Male-Female Counts
mfc <- read.delim("bimm143_05_rstats/male_female_counts.txt")
colors = c(rainbow(nrow(mfc)))
par(mar = c(5,5,5,5))
barplot(mfc$Count, names.arg = mfc$Sample, col = colors, las = 3, cex.names = 0.8, main = "Male and Female Counts", ylab = "Counts")
barplot(mfc$Count, names.arg = mfc$Sample, col = c("red", "blue"), las = 3, cex.names = 0.8, main = "Male and Female Counts", ylab = "Counts")


