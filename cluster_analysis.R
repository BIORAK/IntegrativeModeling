# Load the following libraries
# Load the "fpc" library for cluster statistics
library(fpc)

# Load the "bio3d" library for all vs all pairwise rmsd calculation
library(bio3d)

# Load the "A2R" library for coloring the dendrogram
#source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
#source("/new/lbiorak/R/x86_64-pc-linuxm2[lower.tri(rmsdValues)] <- NA-gnu-library/2.15/A2R.R")

# Read the PDB Filenames and chains from the input file
pdbfiles = read.delim("pdbid_chains.txt",header=F)

# Store the PDB filenames in filenames array
filenames = array()
for(i in 0:length(pdbfiles$V1))
{
  filenames[i] = paste("top_models/", pdbfiles$V1[i], ".pdb", sep = "")
}

pdbs <- pdbaln(filenames)
core <- core.find(pdbs)

gaps.pos <- gap.inspect(pdbs$xyz)

inds <- print(core, vol = 1)
write.pdb(xyz = pdbs$xyz[1, inds$xyz], resno = pdbs$resno[1,inds$atom], file = "quick_core.pdb")
xyz <- fit.xyz(fixed = pdbs$xyz[1, ], mobile = pdbs, fixed.inds = inds$xyz, mobile.inds = inds$xyz, pdbext = "", outpath = "core_fitlsq/", full.pdbs = TRUE, het2atom = TRUE)

# All vs all pairwise rmsd matrix
rmsdValues <- rmsd(xyz[, gaps.pos$f.inds])
rownames(rmsdValues) <- pdbfiles$V1
colnames(rmsdValues) <- pdbfiles$V1

hc <- hclust(as.dist(rmsdValues), method="complete")

# Get the number of clusters in the tree by cutting at a height
clustStats <- cluster.stats(as.dist(rmsdValues), cutree(hc, h = 5), G2=FALSE, G3=FALSE,silhouette=FALSE)

# Print the cluster statistics
clusStat <- c("Cluster Statistics")
underLine <- c("#################")
clusNum <- c("Number of clusters ->", clustStats$cluster.number)
clusAvgBet <- c("Average between clusters ->", clustStats$average.between)

# Color the 
#A2Rplot(hc, k = clustStats$cluster.number,col.up = "gray", col.down = rainbow(clustStats))
