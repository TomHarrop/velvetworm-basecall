# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

summary_file <- snakemake@input[[1]]
reads_to_keep <- snakemake@output[[1]]

full_summary <- fread(summary_file)
fwrite(full_summary[passes_filtering == TRUE, .(read_id)],
       file = reads_to_keep,
       col.names = FALSE)

sessionInfo()
