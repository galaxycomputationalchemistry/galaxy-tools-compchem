#!/usr/bin/env Rscript

options(stringAsfactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

library(bio3d)

dcdfile <- args[1]
pdbfile <- args[2]

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

method <- args[3]
selection <- args[4]
domain <- args[5]

output <- args[6]
pca_plot <- args[7]
pca_cluster  <- args[8]
pc1_rmsf <- args[9]


if (selection == "string") {
    inds <- atom.select(pdb, string = domain)
} 
if (selection == "elety") {
    inds <- atom.select(pdb, elety = domain)
}
if (selection == "resid") {
    inds <- atom.select(pdb, resid = domain)
}
if (selection == "segid") {
    inds <- atom.select(pdb, segid = domain)
}
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

if (method == "FALSE") {
    pc <- pca.xyz(xyz[,inds$xyz], use.svd=FALSE)
}
if (method == "TRUE") {
    pc <- pca.xyz(xyz[,inds$xyz], use.svd=TRUE)
}

write.table(pc$au[,1:2:3], file = output, row.names = TRUE, col.names = FALSE, quote =FALSE, sep="\t")

png(pca_plot)
plot(pc, col=bwr.colors(nrow(xyz)) )
dev.off()

png(pca_cluster)
hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
plot(pc, col=grps)
dev.off()

png(pc1_rmsf)
plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc$au[,2], typ="l", col="blue")
dev.off()

