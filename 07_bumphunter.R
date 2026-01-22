########################## file information #############################
# Author: Bochao Zhang
# Date: 2024-9-12
# Description:
# Use bumphunter to find regions with allele imbalance
#########################################################################

options(scipen=999)

library(optparse)
library(stringr)

## Parse options
option_list <- list(
	make_option(c("-i", "--input"), type = "character", default = NULL,
    	help = "Path to input file"),
	make_option("--cluster_hap", type = "integer", default = 3000,
		help = "Maximum gap allowed when clustering"),
	make_option("--normal_span", type = "integer", default = 10000,
		help = "Normal copy number smooth span"),
	make_option("--tumor_span", type = "integer", default = 200,
		help = "Tumor copy number smooth span"),
	make_option(c("-c", "--cutoff"), type = "numeric", default = 0.7,
		help = "Copy number diffence cutoff")
	make_option(c("-o", "--output"), type = "character", default = NULL,
		help = "Output directory")
)

opt_parser = OptionParser(option_list=option_list, usage="Use bumphunter to find regions with allele imbalance")
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Specify input file using -i or --input.n", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("Specify output directory using -o or --output.n", call.=FALSE)
}

library(bumphunter, quietly=TRUE)

## Read input table
infile <- sprintf(opt$input)
data <- read.csv(infile)

## Sort data based on position
pos <- list(pos1=strtoi(data$pos))
chr <- rep(paste0("chr",c(6)), times = sapply(pos,length))
pos <- unlist(pos, use.names=FALSE)
idx <- order(pos)
pos <- pos[idx]

## Create output data
outdata <- data.frame(pos=pos)
outdata$tumor1 <- data$tumor1[idx]
outdata$tumor2 <- data$tumor2[idx]
outdata$normal1 <- data$normal1[idx]
outdata$normal2 <- data$normal2[idx]

## Clustering
cl <- clusterMaker(chr, pos, maxGap = opt$cluster_hap)

## Smooth normal copy number
n1 <- data$normal1
n1 <- n1[idx]
smoothed <- loessByCluster(n1, pos, cluster=cl, bpSpan=opt$normal_span, minNum=7, minInSpan=5)
for (i in 1:length(smoothed$smoothed)) {
    if (smoothed$smoothed[i][1]) {
        n1[i] <- smoothed$fitted[i]
    }
}

n2 <- data$normal2
n2 <- n2[idx]
smoothed <- loessByCluster(n2, pos, cluster=cl, bpSpan=opt$normal_span, minNum=7, minInSpan=5)
for (i in 1:length(smoothed$smoothed)) {
    if (smoothed$smoothed[i][1]) {
        n2[i] <- smoothed$fitted[i]
    }
}

outdata$smoothed_normal1 <- n1
outdata$smoothed_normal2 <- n2

## Segment hap1 tumor copy number
t1 <- data$tumor1
t1 <- t1[idx]

segs <- getSegments(t1, cl, cutoff=c(0.7,1.3))

up_idx <- c()
dn_idx <- c()
zero_idx <- c()
for (ind in segs$upIndex){
    up_idx <- append(up_idx, ind)
}
for (ind in segs$dnIndex){
    dn_idx <- append(dn_idx, ind)
}
for (ind in segs$zeroIndex){
    zero_idx <- append(zero_idx, ind)
}

## Smooth hap1 tumor copy number
smoothed <- loessByCluster(t1[up_idx], pos[up_idx], cluster=cl[up_idx], bpSpan=opt$tumor_span, minNum=2, minInSpan=2)
for (i in 1:length(smoothed$smoothed)) {
    if (smoothed$smoothed[i][1]) {
        t1[up_idx[i]] <- smoothed$fitted[i]
    }
}

smoothed <- loessByCluster(t1[dn_idx], pos[dn_idx], cluster=cl[dn_idx], bpSpan=opt$tumor_span, minNum=2, minInSpan=2)
for (i in 1:length(smoothed$smoothed)) {
    if (smoothed$smoothed[i][1]) {
        t1[dn_idx[i]] <- smoothed$fitted[i]
    }
}

smoothed <- loessByCluster(t1[zero_idx], pos[zero_idx], cluster=cl[zero_idx], bpSpan=opt$tumor_span, minNum=2, minInSpan=2)
for (i in 1:length(smoothed$smoothed)) {
    if (smoothed$smoothed[i][1]) {
        t1[zero_idx[i]] <- smoothed$fitted[i]
    }
}

outdata$smoothed_tumor1 <- t1

## Segment hap2 tumor copy number
t2 <- data$copynumber2
t2 <- t2[idx]

segs <- getSegments(t2, cl, cutoff=c(0.7,1.3))

up_idx <- c()
dn_idx <- c()
zero_idx <- c()
for (ind in segs$upIndex){
    up_idx <- append(up_idx, ind)
}
for (ind in segs$dnIndex){
    dn_idx <- append(dn_idx, ind)
}
for (ind in segs$zeroIndex){
    zero_idx <- append(zero_idx, ind)
}

## Smooth hap2 tumor copy number
smoothed <- loessByCluster(t2[up_idx], pos[up_idx], cluster=cl[up_idx], bpSpan=opt$tumor_span, minNum=2, minInSpan=2)
for (i in 1:length(smoothed$smoothed)) {
    if (smoothed$smoothed[i][1]) {
        t2[up_idx[i]] <- smoothed$fitted[i]
    }
}

smoothed <- loessByCluster(t2[dn_idx], pos[dn_idx], cluster=cl[dn_idx], bpSpan=opt$tumor_span, minNum=2, minInSpan=2)
for (i in 1:length(smoothed$smoothed)) {
    if (smoothed$smoothed[i][1]) {
        t2[dn_idx[i]] <- smoothed$fitted[i]
    }
}

smoothed <- loessByCluster(t2[zero_idx], pos[zero_idx], cluster=cl[zero_idx], bpSpan=opt$tumor_span, minNum=2, minInSpan=2)
for (i in 1:length(smoothed$smoothed)) {
    if (smoothed$smoothed[i][1]) {
        t2[zero_idx[i]] <- smoothed$fitted[i]
    }
}

outdata$smoothed_tumor2 <- t1

## Normalize allele copy number difference
normalized_diff <- (t1-t2) / (t1+t2) * 2
outdata$normalized_diff <- diff

## Output copy numbers
outfile <- sprintf("%s/smoothed_cn.csv", opt$output)
write.csv(outdata, outfile, row.names=FALSE, quote=FALSE)

## Find copy number differences bumps
mat <- matrix(c(rep(0L, length(normalized_diff)), normalized_diff), nrow=length(normalized_diff), ncol=2)
design <- cbind(rep(1,2), rep(c(0,1), each=1))
bumps <- bumphunter(mat, design, chr, pos, cl, coef=2, cutoff=opt.cutoff, bpSpan=opt$cluster_hap, nullMethod="permutation", smooth=FALSE, B=250, verbose=FALSE)

## Output bumps
outfile <- sprintf("%s/bumps.csv", opt$output)
write.csv(bumps$table, outfile, row.names=FALSE, quote=FALSE)
