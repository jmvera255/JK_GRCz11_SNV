
# read in command line arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1){
  write(paste("No samples list provided!\n",
              "Usage: Rscript calc_homozygous_score.R <samples.txt>\n", 
              "Where samples.txt is a list of header names for merged vcf files, one per line\n",
              sep = "\n"), 
        stderr())
  quit(status=1, save="no");
}

# Load libraries silently
library(stringr)
library(dplyr)
library(ggplot2)

# Process Count Results
samples <- scan(args[1], character(), sep="\n");

hits.list <- list()
for(k in c("10", "25", "50")){
  hits.list[[k]] <- as.data.frame(matrix(NA, ncol=10, nrow=length(samples)))
  for(s in 1:length(samples)){
    sample <- samples[s]
  
    fid <- paste0(k, "kb/", sample, ".", k, "kb.counts.txt")
    counts <- read.delim(fid, header=FALSE, sep="\t", stringsAsFactors = FALSE)
    colnames(counts) <- c("Region", "WT.het", "Mut.het", "Mut.hom")
    
    counts <- mutate(counts, Het.ratio = WT.het/Mut.het)
    counts$Het.ratio[is.infinite(counts$Het.ratio)] <- 1
    counts$Het.ratio[is.na(counts$Het.ratio)] <- 1
    counts <- mutate(counts, Score = Het.ratio * Mut.hom)
    
    # prep df for plotting
    counts <- mutate(counts, Chr = as.numeric(str_split(Region, ":", simplify = TRUE)[,1]))
    temp <- str_split(counts$Region, ":", simplify = TRUE)[,2]
    counts$Start <- as.numeric(str_split(temp, "-", simplify = TRUE)[,1])
    counts <- arrange(counts, Chr, Start)
    counts$Pos <- seq(1, nrow(counts), by = 1)
    h <- which.max(counts$Score)
    hits.list[[k]][s,] <- c(sample, counts[h,])
    
    # get chrom positions
    chrom.pos <- counts %>% group_by(Chr) %>% 
      summarise(Start = min(Pos), Stop = max(Pos), Mid = median(Pos)) %>% as.data.frame
    #x.ticks <- sort(c(chrom.pos$Mid, chrom.pos$Start[-1]))
    x.ticks <- chrom.pos$Start[-1]
    x.labels <- rep("", times=length(x.ticks))
    label.positions <- seq(1, length(x.ticks), by = 2)
    for(x in 1:length(label.positions)){
      px <- label.positions[x]
      x.labels[px] <- chrom.pos$Chr[x] 
    }
    
    # make plot
    fid <- paste0(sample, ".", k, "kb.genome-wide_score.png")
    png(fid, height=450, width = 1000)
    p <- ggplot(counts, aes(x=Pos, y=Score)) + 
      geom_line(color="dodgerblue3") +
      geom_point(color="dodgerblue3", show.legend = FALSE) + 
      labs(title=paste0(sample, " ", k, "kb Windows"), y="Homozygosity Score", x="\nChromosome") +
      theme_classic() + scale_y_continuous(expand = c(0.007,0.007)) +
      scale_x_continuous(expand = c(0.005,0.005), breaks = x.ticks, labels = NULL) +
      theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"), plot.margin = unit(c(1,1,3,1), "lines"),
            axis.line = element_line(size = 1, linetype = "solid"), 
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10), 
            axis.title.x = element_text(vjust=-0.5), 
            axis.ticks.length=unit(.25, "cm")) +
      annotate(geom="text", x = counts$Pos[h], y = counts$Score[h]*1.04, label=counts$Region[h],
               color="black")
    
    # add x-axis chrom labels btwn tick marks
    for(i in 1:nrow(chrom.pos)){
      p <- p + annotate(geom="text", x=chrom.pos$Mid[i], y=-1 * counts$Score[h]/30, label=as.character(chrom.pos$Chr[i]))
    }
    p <- p + coord_cartesian(xlim=c(1,max(counts$Pos)*1.05), ylim = c(0, max(counts$Score)*1.07),clip = "off")
    
    print(p)
    invisible(dev.off())
  
  }
}

# process hits.list
for(k in c("10", "25", "50")){
  df <- hits.list[[k]]
  colnames(df) <- c("Sample", "Region", "WT.het", "Mut.het", "Mut.hom", 
                    "Het.ratio", "Score", "Chr", "Start", "Pos")
  df <- select(df, -Pos)
  # set start/stop using BED file formats and write to file
  df$Start <- df$Start - 1
  df$Stop <- df$Start + as.numeric(k)*1000
  df$kb <- as.numeric(k)
  fid <- paste0(k, "kb.top_score_regions.txt")
  write.table(df, file=fid, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  hits.list[[k]] <- df
}

# print top scoring regions for each window size
for(k in c("10", "25", "50")){
  df <- hits.list[[k]][,1:7]
  print(df)
}
