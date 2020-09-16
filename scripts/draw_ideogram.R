library("RIdeogram")


dual_comp <- read.table(snakemake@input[[1]], sep = "\t", header = T, stringsAsFactors = F)
synteny_comp <- read.table(snakemake@input[[2]], sep="\t", header=T, stringsAsFactors = F)

k <- strsplit(snakemake@input[[3]], "/")
z <- k[[1]][1:length(k[[1]])-1]
p <- paste(z, collapse="/")
setwd(p)

ideogram(karyotype=dual_comp,synteny=synteny_comp)
convertSVG("chromosome.svg", device = "png")
