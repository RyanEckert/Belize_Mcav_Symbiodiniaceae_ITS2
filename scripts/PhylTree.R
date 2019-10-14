setwd("~/path/to/working/directory")
set.seed(10)

library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")

cbcGen = read.csv("CBCgen3.csv")
CBCgen = df2genind(cbcGen[2:10], ind.names = cbcGen$Sample, sep = ",")
CBCgen

CBCdist = provesti.dist(CBCgen)
CBCdist

theTree = CBCdist %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
tree = plot(theTree)
add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.


