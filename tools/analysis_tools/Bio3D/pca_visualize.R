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

pdb1 <- args[6]
nc_pc1 <- args[7]
pdb2 <- args[8]
nc_pc2 <- args[9]
pdb3 <- args[10]
nc_pc3 <- args[11]

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

p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file=pdb1)
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file=pdb2)
p3 <- mktrj.pca(pc, pc=3,b=pc$au[,3], file=pdb3)

write.ncdf(p1, nc_pc1)
write.ncdf(p2, nc_pc2)
write.ncdf(p3, nc_pc3)
