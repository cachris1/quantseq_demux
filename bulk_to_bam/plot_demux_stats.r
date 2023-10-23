library(dplyr)
library(ggplot2)


do_something <- function(summary_tsv, out_path) {
    demux_stats <- read.table()
    p <- ggplot(demux_stats, aes(y = written_reads)) + geom_histogram()
    ggsave(p, out_path)
}

do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads)