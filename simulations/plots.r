
require (graphics)
require (CarbonEL)
require (stats)

toPDF = TRUE

alignColor = "blue"
structureColor = "red"

alignShape = 16
structureShape = 17

## default plot filename
plotfile = "simulations.pdf"
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
structure = array()
bins = unique (data$z)
for (bin in bins) {
  align = rbind (
    align,
    c (bin, mean (subset (data$AlignAcc,
                          data$z == bin))))
  structure = rbind (
    structure,
    c (bin, mean (subset (data$StructureReconstruction,
                          data$z == bin))))
}

## remove NA rows
align = align[-1,]
structure = structure[-1,]

plot.new()
plot.window (xlim = c (0, 2.5),
             ylim = c (0.4, 1))
axis (1,
      c (0, 0.5, 1, 1.5, 2, 2.5))
axis (2)
title (xlab = "Outgroup branch length",
       ylab = "Inference accuracy",
       cex.lab = 1.2)

points (align,
        col = alignColor, pch = alignShape, cex = 1)
lines (lowess (align), col = alignColor, lwd = 1.5)
points (structure,
        col = structureColor, pch = structureShape, cex = 1)
lines (lowess (structure), col = structureColor, lwd = 1.5)

legend ("bottomleft", legend = c ("Structural reconstruction accuracy", "Alignment accuracy"),
        col = c (structureColor, alignColor),
        pch = c (structureShape, alignShape),
        xjust = 1,
        inset = 0.02)

if (toPDF) {
  dev.off()
}
