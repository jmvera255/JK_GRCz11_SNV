
# read in command line arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2){
  write(paste("Insufficient arguments provided!\n",
              "Usage: Rscript plot-stats.R  <samples_name> <query_stats.txt>\n", 
              sep = "\n"), 
        stderr())
  quit(status=1, save="no");
}
library(ggplot2)

# define sample name and input file name from command line
sample <- as.character(args[1])

# build df of data
data <- read.delim(args[2], sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(data) <- c("Qual", "RMS_Qual", "Allele_Depth", "Total_Depth")

# define quantile targets for ablines
quants <- c(25, NA, 75, 75)
lims <- list(c(0,1000), c(10,50), c(0, 40),c(0,50))
ticks <- list(seq(100,1000, by=100), seq(15,50, by=5),
              seq(5,40, by=5), seq(5, 50, by=5))

# plot
fid <- paste0(sample, ".stats_summary.png")
png(fid, width=800, height=600)
par(mfrow=c(2,2), mar=c(4,5,3,3), oma=c(2,2,2,2))

for(m in 1:ncol(data)){
  metric <- colnames(data)[m]
  d <- density(data[,m])
  plot(d, main=paste(sample, metric), xlim=lims[[m]], xaxt="n", ylim=c(0, max(d[['y']])*1.2), cex = 1.15)
  polygon(d, col="grey65", border="black")
  
  if(!is.na(quants[m])){
    q <- quantile(data[,m], quants[m]/100)
    abline(v=q, col="blue2", lwd=2)

    legend("topright", col="blue2", lwd=1, bty="n", 
           legend=paste0(" ", quants[m], "th \nquant \n (", q, ") "))
  }
  axis(1, at = ticks[[m]], labels=TRUE, pos=0)
}
dev.off()

# melt data and save for later plotting w/ggplot and all samples
data$Sample <- sample
melt.data <- reshape2::melt(data, id.vars="Sample")
sample <- stringr::str_replace_all(sample, "-", ".")
assign(paste0(sample, ".melt"), melt.data)
rm(melt.data)
save(list = ls(pattern="melt"), file = paste0(sample, ".melt.RData"))

