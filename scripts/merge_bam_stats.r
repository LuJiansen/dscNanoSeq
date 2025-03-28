library(dplyr)
library(data.table)

args <- commandArgs(TRUE)
indir <- args[1]
pattern <- args[2]
output <- args[3]

stats_files <- list.files(path = indir, pattern = pattern, full.names = T)

stats <- lapply(stats_files, function(x) {
    df <- fread(x)
    colnames(df) <- c("barcode", "nFrags", "datasize", "min_len", "max_len", "mean_len", "median_len", "MAQ20", "MAQ30")
    return(df)
}) %>% do.call(rbind, .)

head(stats)

stats_m <- stats %>%
    dplyr::group_by(barcode) %>%
    dplyr::summarise(
        min_len = min(min_len),
        max_len = max(max_len),
        median_len = mean(median_len),
        mean_len = sum(mean_len * nFrags) / sum(nFrags),
        MAQ20 = sum(MAQ20 * nFrags) / sum(nFrags),
        MAQ30 = sum(MAQ30 * nFrags) / sum(nFrags),
        nFrags = sum(nFrags),
        datasize = sum(datasize)
    )

write.table(stats_m, output, row.names = F, col.names = T, quote = F, sep = '\t')

