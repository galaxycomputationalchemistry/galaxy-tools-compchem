#!/usr/bin/env Rscript

options(stringAsfactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

library(bio3d)

dcdfile <- args[1]
pdbfile <- args[2]
selection <- args[3]

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)


if (selection == "string") {
    domain <- args[4]
    output <- args[5]
    rmsf_plot <- args[6]
    inds <- atom.select(pdb, string = domain)
} 
if (selection == "resno") {
    res1 <- args[4]
    res2 <- args[5]
    output <- args[6]
    rmsf_plot <- args[7]
    inds <- atom.select(pdb, resno=res1:res2)
} 
if (selection == "elety") {
    domain <- args[4]
    output <- args[5]
    rmsf_plot <- args[6]
    inds <- atom.select(pdb, elety = domain)
}
if (selection == "resid") {
    domain <- args[4]
    output <- args[5]
    rmsf_plot <- args[6]
    inds <- atom.select(pdb, resid = domain)
}
if (selection == "segid") {
    domain <- args[4]
    output <- args[5]
    rmsf_plot <- args[6]
    inds <- atom.select(pdb, segid = domain)
}

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

rf <- rmsf(xyz[,inds$xyz])

write.table(rf, file = output, row.names = TRUE, col.names = FALSE, quote =FALSE, sep="\t")

png(rmsf_plot)
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")
dev.off()

