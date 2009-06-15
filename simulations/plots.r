
require (graphics)
## require (CarbonEL)
require (stats)

toPDF = TRUE

alignColor = "red"
alignTKFColor = "green"
alignLongColor = "blue"
alignStemlocColor = "yellow"
alignStemlocAmaColor = "purple"

structureColor = "red"

alignShape = 16
alignTKFShape = 15
alignLongShape = 1
alignStemlocShape = 2
alignStemlocAmaShape = 20

structureShape = 16

## default plot filename
plotfile = "ama.pdf"
plotdim = c (6, 4)

## open PDF
if (toPDF) {
  pdf (file = plotfile, width = plotdim[1], height = plotdim[2])
} else {
  quartz (width = plotdim[1], height = plotdim[2])
}

## set up plot
par (family = "serif")
par (cex = 1)

par (oma = c (0, 0, 0, 0))

## read data
data = read.table ("inferred/indiegram.analysis", header = TRUE, sep = "\t")

## average datapoints for each outgroup length

align = array()
alignTKF = array()
alignLong = array()
alignStemloc = array()
alignStemlocAma = array()

alignTCS = array()
alignTCSTKF = array()
alignTCSLong = array()
alignTCSStemloc = array()
alignTCSStemlocAma = array()

structure = array()

bins = unique (data$z)
for (bin in bins) {

  align = rbind (
    align,
    c (bin, mean (subset (data$AlignAcc,
                          data$z == bin))))
  alignTKF = rbind (
    alignTKF,
    c (bin, mean (subset (data$AlignAccTKF,
                          data$z == bin))))
  alignLong = rbind (
    alignTKF,
    c (bin, mean (subset (data$AlignAccLong,
                          data$z == bin))))
  alignStemloc = rbind (
    alignStemloc,
    c (bin, mean (subset (data$AlignAccStemloc,
                          data$z == bin))))
  alignStemlocAma = rbind (
    alignStemlocAma,
    c (bin, mean (subset (data$AlignAccStemlocAma,
                          data$z == bin))))

  alignTCS = rbind (
    alignTCS,
    c (bin, mean (subset (data$AlignTCS,
                          data$z == bin))))
  alignTCSTKF = rbind (
    alignTCSTKF,
    c (bin, mean (subset (data$AlignTCSTKF,
                          data$z == bin))))
  alignTCSLong = rbind (
    alignTCSTKF,
    c (bin, mean (subset (data$AlignTCSLong,
                          data$z == bin))))
  alignTCSStemloc = rbind (
    alignTCSStemloc,
    c (bin, mean (subset (data$AlignTCSStemloc,
                          data$z == bin))))
  alignTCSStemlocAma = rbind (
    alignTCSStemlocAma,
    c (bin, mean (subset (data$AlignTCSStemlocAma,
                          data$z == bin))))

  structure = rbind (
    structure,
    c (bin, mean (subset (data$StructureReconstruction,
                          data$z == bin))))
}

## remove NA rows
align = align[-1,]
alignTKF = alignTKF[-1,]
alignLong = alignLong[-1,]
alignStemloc = alignStemloc[-1,]
alignStemlocAma = alignStemlocAma[-1,]

alignTCS = alignTCS[-1,]
alignTCSTKF = alignTCSTKF[-1,]
alignTCSLong = alignTCSLong[-1,]
alignTCSStemloc = alignTCSStemloc[-1,]
alignTCSStemlocAma = alignTCSStemlocAma[-1,]

structure = structure[-1,]

plot.new()
plot.window (xlim = c (0, 2.5),
             ylim = c (0, 1))
axis (1,
      c (0, 0.5, 1, 1.5, 2, 2.5))
axis (2)
title (xlab = "Outgroup branch length",
       ylab = "Alignment metric accuracy (AMA)",
       cex.lab = 1.2)

points (align,
        col = alignColor, pch = alignShape, cex = 1)
lines (lowess (align), col = alignColor, lwd = 1.5)


points (alignTKF,
        col = alignTKFColor, pch = alignTKFShape, cex = 1)
lines (lowess (alignTKF), col = alignTKFColor, lwd = 1.5)


points (alignLong,
        col = alignLongColor, pch = alignLongShape, cex = 1)
lines (lowess (alignLong), col = alignLongColor, lwd = 1.5)


points (alignStemloc,
        col = alignStemlocColor, pch = alignStemlocShape, cex = 1)
lines (lowess (alignStemloc), col = alignStemlocColor, lwd = 1.5)


points (alignStemlocAma,
        col = alignStemlocAmaColor, pch = alignStemlocAmaShape, cex = 1)
lines (lowess (alignStemlocAma), col = alignStemlocAmaColor, lwd = 1.5)



legend (x = "bottomleft",
        legend = c ("Indiegram", "Stemloc", "Stemloc-AMA",  "TKF91", "Long Indel"),
        col = c (alignColor, alignStemlocColor, alignStemlocAmaColor, alignTKFColor, alignLongColor),
        pch = c (alignShape, alignStemlocShape, alignStemlocAmaShape, alignTKFShape, alignLongShape),
        xjust = 1,
        inset = 0.02)

if (toPDF) {
  dev.off()
}


# TCS
## default plot filename
plotfile = "tcs.pdf"
plotdim = c (6, 4)

## open PDF
if (toPDF) {
  pdf (file = plotfile, width = plotdim[1], height = plotdim[2])
} else {
  quartz (width = plotdim[1], height = plotdim[2])
}

## set up plot
par (family = "serif")
par (cex = 1)

par (oma = c (0, 0, 0, 0))

plot.new()
plot.window (xlim = c (0, 2.5),
             ylim = c (0, 1))
axis (1,
      c (0, 0.5, 1, 1.5, 2, 2.5))
axis (2)
title (xlab = "Outgroup branch length",
       ylab = "Total column score (TCS)",
       cex.lab = 1.2)

points (alignTCS,
        col = alignColor, pch = alignShape, cex = 1)
lines (lowess (alignTCS), col = alignColor, lwd = 1.5)


points (alignTCSTKF,
        col = alignTKFColor, pch = alignTKFShape, cex = 1)
lines (lowess (alignTCSTKF), col = alignTKFColor, lwd = 1.5)


points (alignTCSLong,
        col = alignLongColor, pch = alignLongShape, cex = 1)
lines (lowess (alignTCSLong), col = alignLongColor, lwd = 1.5)


points (alignTCSStemloc,
        col = alignStemlocColor, pch = alignStemlocShape, cex = 1)
lines (lowess (alignTCSStemloc), col = alignStemlocColor, lwd = 1.5)


points (alignTCSStemlocAma,
        col = alignStemlocAmaColor, pch = alignStemlocAmaShape, cex = 1)
lines (lowess (alignTCSStemlocAma), col = alignStemlocAmaColor, lwd = 1.5)



legend (x = "bottomleft",
        legend = c ("Indiegram", "Stemloc", "Stemloc-AMA",  "TKF91", "Long Indel"),
        col = c (alignColor, alignStemlocColor, alignStemlocAmaColor, alignTKFColor, alignLongColor),
        pch = c (alignShape, alignStemlocShape, alignStemlocAmaShape, alignTKFShape, alignLongShape),
        xjust = 1,
        inset = 0.02)

if (toPDF) {
  dev.off()
}



# structure accuracy
plotfile = "structure.pdf"
plotdim = c (6, 4)

## open PDF
if (toPDF) {
  pdf (file = plotfile, width = plotdim[1], height = plotdim[2])
} else {
  quartz (width = plotdim[1], height = plotdim[2])
}

## set up plot
par (family = "serif")
par (cex = 1)

par (oma = c (0, 0, 0, 0))

## read data
data = read.table ("inferred/indiegram.analysis", header = TRUE, sep = "\t")


plot.new()
plot.window (xlim = c (0, 2.5),
             ylim = c (0, 1))
axis (1,
      c (0, 0.5, 1, 1.5, 2, 2.5))
axis (2)
title (xlab = "Outgroup branch length",
       ylab = "Accuracy of ancestral structure",
       cex.lab = 1.2)

points (structure,
        col = structureColor, pch = structureShape, cex = 1)
lines (lowess (structure), col = structureColor, lwd = 1.5)


if (toPDF) {
  dev.off()
}




# count exactly correct alignments; bin every 5 z-coords
delta = bins[2] - bins[1]
chunks = 5
binwidth = delta * chunks
bigbins = unique (round ((bins / delta) / chunks)) * binwidth + binwidth / 2

# manually remove last bin. Number of bins is hardcoded here.... UGH
bigbins = bigbins[-6]

# get bin names
#binnamefun <- function(x) cat (x - binwidth/2, " to ", x + binwidth/2)
#binname = apply(bigbins, 1, binnamefun)
binname = c("0 to 0.5", "0.5 to 1", "1 to 1.5", "1.5 to 2", "2 to 2.5");

# do binning
exact = array()
for (bin in bigbins) {

  bmin = bin - binwidth / 2
  bmax = bin + binwidth / 2

  exact = cbind (
    exact,
    c (sum (subset (data$AlignAcc,
                          data$z >= bmin & data$z < bmax & data$AlignAcc == 1)),
 
       sum (subset (data$AlignAccStemloc,
                          data$z >= bmin & data$z < bmax & data$AlignAccStemloc == 1)),

       sum (subset (data$AlignAccStemlocAma,
                          data$z >= bmin & data$z < bmax & data$AlignAccStemlocAma == 1)),	

       sum (subset (data$AlignAccTKF,
                          data$z >= bmin & data$z < bmax & data$AlignAccTKF == 1)),

       sum (subset (data$AlignAccLong,
                          data$z >= bmin & data$z < bmax & data$AlignAccLong == 1))))
}

exact = exact[,-1]


# histogram
plotfile = "perfect.pdf"
plotdim = c (6, 4)

## open PDF
if (toPDF) {
  pdf (file = plotfile, width = plotdim[1], height = plotdim[2])
} else {
  quartz (width = plotdim[1], height = plotdim[2])
}

## set up plot
par (family = "serif")
par (cex = 1)

par (oma = c (0, 0, 0, 0))

barplot(exact,
	names.arg = binname,
	beside = TRUE,
	xlab = "Outgroup branch length",
	ylab = "Perfect alignments",
        col = c (alignColor, alignStemlocColor, alignStemlocAmaColor, alignTKFColor, alignLongColor))

legend (x = "topright",
        legend = c ("Indiegram", "Stemloc", "Stemloc-AMA",  "TKF91", "Long Indel"),
        col = c (alignColor, alignStemlocColor, alignStemlocAmaColor, alignTKFColor, alignLongColor),
        pch = c (15, 15, 15, 15, 15),
        xjust = 1,
        inset = 0.02)


if (toPDF) {
  dev.off()
}
