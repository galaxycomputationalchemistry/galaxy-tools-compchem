#!/usr/bin/env Rscript

options(stringAsfactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

library(bio3d)

dcdfile <- args[1]
pdbfile <- args[2]

output <- args[3]
dccm_plot <- args[4]

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)

cij<-dccm(xyz[,ca.inds$xyz])

write.table(cij, file = output, row.names = TRUE, col.names = FALSE, quote =FALSE, sep="\t")

png(dccm_plot)
plot(cij)
dev.off()


