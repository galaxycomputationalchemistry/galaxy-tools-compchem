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
    rmsd_plot <- args[6]
    rmsd_hist <- args[7]
    inds <- atom.select(pdb, string = domain)
} 
if (selection == "resno") {
    res1 <- args[4]
    res2 <- args[5]
    output <- args[6]
    rmsd_plot <- args[7]
    rmsd_hist <- args[8]
    inds <- atom.select(pdb, resno=res1:res2)
} 
if (selection == "elety") {
    domain <- args[4]
    output <- args[5]
    rmsd_plot <- args[6]
    rmsd_hist <- args[7]
    inds <- atom.select(pdb, elety = domain)
}
if (selection == "resid") {
    domain <- args[4]
    output <- args[5]
    rmsd_plot <- args[6]
    rmsd_hist <- args[7]
    inds <- atom.select(pdb, resid = domain)
}
if (selection == "segid") {
    domain <- args[4]
    output <- args[5]
    rmsd_plot <- args[6]
    rmsd_hist <- args[7]
    inds <- atom.select(pdb, segid = domain)
}

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

rd <- rmsd(xyz[1,inds$xyz], xyz[,inds$xyz])

write.table(rd, file = output, row.names = TRUE, col.names = FALSE, quote =FALSE, sep="\t")

png(rmsd_plot)
plot(rd, typ="l", ylab="RMSD (Å)", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)
dev.off()

png(rmsd_hist)
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
lines(density(rd), typ="l", col="red", lty=2, lwd=2)
dev.off()
