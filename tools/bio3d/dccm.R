#!/usr/bin/env Rscript

options(stringAsfactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

library(bio3d)
require(lattice)

dcdfile <- args[1]
pdbfile <- args[2]

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

selection <- args[3]
domain <- args[4]

output <- args[5]
dccm_plot <- args[6]

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

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
cij<-dccm(xyz[,inds$xyz])

write.table(cij, file = output, row.names = TRUE, col.names = FALSE, quote =FALSE, sep="\t")

png(dccm_plot)
plot(cij)
dev.off()


