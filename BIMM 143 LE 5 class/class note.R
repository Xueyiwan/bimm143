#' ---
#' title: "Class5 Data exploration and visualization in R"
#' author: "Xueyi Wan"
#' output: github_document
#' ---

# Class5 Data visualization
x <- rnorm(1000)
# gives 1000 random normal distribution data points

# some summary stats
mean(x)
sd(x)

summary(x)
boxplot(x)

hist(x)
rug(x)


# read.table instructions
#   read.csv and read.csv2 are identical to read.table except for the 
# defaults. They are intended for reading ‘comma separated value’ 
# files (‘.csv’) or (read.csv2) the variant used in countries that 
# use a comma as decimal point and a semicolon as field separator. 
# Similarly, read.delim and read.delim2 are for reading delimited 
# files, defaulting to the TAB character for the delimiter. 

# Section 2 scatter plots
# lets read our input file

weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
plot(weight$Age, weight$Weight, type = "o", pch = 15, cex = 1.5, lwd = 2, 
     ylim = c(2,10), xlab = "Age(months)", ylab = "Weight(kg)", 
     main = "Baby weight with age", col = "blue")
plot(1:5, pch=1:5, cex=1:5)

# Section 2 Barplot

mouse <- read.delim("bimm143_05_rstats/feature_counts.txt")
par(mar = c(3.1, 11.1, 4.1, 2))
barplot(mouse$Count, horiz = TRUE, names.arg = mouse$Feature, 
        main = "Number of features in the mouse GRCm38 genome", 
        las = 1, xlim = c(0,80000))

# Section 2 Histograms

x <- c(rnorm(10000), rnorm(10000)+4)
hist(x, breaks = 80)

# Section 3 Providing color vectors

people <- read.delim("bimm143_05_rstats/male_female_counts.txt", 
                     header = TRUE)
par(mar = c(5,5,2,2))
# make all bars different colors using the rainbow() function
# note the nrow() function
barplot(people$Count, names.arg = people$Sample, las = 2, 
        ylab = "Counts", cex.names = 0.7, 
        col = rainbow(nrow(people)))
# alternate male & female samples
# pass a 2 color vector to eh col parameter
barplot(people$Count, names.arg = people$Sample, las = 2, 
        ylab = "Counts", cex.names = 0.7, 
        col = c("blue2", "red2"))

# Section 3 coloring by value

genes <- read.delim("bimm143_05_rstats/up_down_expression.txt")
# How many genes are detailed in this file (how many rows)
nrow(genes)
# To determine how many genes are up, down and unchanging
table(genes$State)

plot(genes$Condition1, genes$Condition2, col = genes$State, 
    xlab = "Expression condition 1", ylab = "Expression condition 2")
# to know the order of states levels
levels(genes$State)
# using the correct colors in the correct order
palette(c("blue", "grey", "red"))

# Section 3 Dynamic use of color

meth <- read.delim("bimm143_05_rstats/expression_methylation.txt")
# How many genes are in this dataset
nrow(meth)

# use densCols() function to make a new color vector 
dcols <- densCols(meth$gene.meth, meth$expression)
# plot changing the plot character ('pch') to a solid circle
plot(meth$gene.meth, meth$expression, col = dcols, pch = 20)

# retrict to the genes that have more than zero expression values
# Find the indices of genes with above 0 expression
inds <- meth$expression > 0
# Make a density color vector for these genes
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds])
# plot just these genes
plot(meth$gene.meth[inds], meth$expression[inds], col = dcols, 
     pch = 20)

# change the colramp used by the densCols() function to go between
# blue, green, red and yellow with the colorRampPalette() function
dcols.custom <- densCols(meth$gene.meth[inds],
                         meth$expression[inds], 
                         colramp = colorRampPalette(c("blue2", 
                                                    "green2", 
                                                    "red2", 
                                                    "yellow")))
plot(meth$gene.meth[inds], meth$expression[inds], 
     col = dcols.custom, pch = 20)


# About this document
sessionInfo()
