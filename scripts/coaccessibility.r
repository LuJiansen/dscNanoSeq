library(parallel)
library(ArchR)
library(gUtils)
library(pbmcapply)
library(dplyr)

args <- commandArgs(TRUE)

#setwd("./analysis/GM12878_merged_noerror")
# frag_file <- "./GM12878_merged_fragment.bed.gz"
# peak_file <- "./dscNanoATAC_GM12878_noerror_ArchR_peaks_local500k_q0.01_rep2.txt"

frag_file <- args[1]
peak_file <- args[2]
ncores <- args[3]
output <- args[4]

# Load fragment and peak data
print("Loading fragment and peak data...")
frag <- read.table(frag_file, comment.char = "")
colnames(frag) <- c("chrom", "start", "end", "id", "mapq", "strand")

frag <- frag %>% arrange(chrom, start)
frag$length <- frag$end - frag$start
frag$id <- make.unique(frag$id)

# Convert fragment to flanking sites
frag_flank <- rbind(
    frag[, c("chrom", "start", "id", "length")] %>% `colnames<-`(c("chrom", "start", "id", "length")) %>% mutate(side = "+"),
    frag[, c("chrom", "end", "id", "length")] %>% `colnames<-`(c("chrom", "start", "id", "length")) %>% mutate(side = "-")
)
frag_flank_gr <- dt2gr(frag_flank)

# + means the left side of fragment is slected
# - means the right side of fragment is selected
# when considering the relative position of peaks,
# + mean the fragment located in right side of peak,
# - mean the fragment located in left side of peak.

peak <- fread(peak_file)
colnames(peak)[1:3] <- c("chrom", "start", "end")
# peak <- peak[peak$chrom == chrom, ]
peak$idx <- 1:nrow(peak)
peak_gr <- dt2gr(peak)

target_gr <- resize(peak_gr, width = 1e3, fix = "center")
control_gr <- resize(peak_gr, width = 1e5, fix = "center")

run_ks_test <- function(frag, target, control, side, cores = 30){
    pbmclapply(unique(target$idx), function(idx) {
        target_ov <- frag[frag$side == side] %*% target[target$idx == idx]
        if (length(target_ov) > 0) {
            control_ov <- frag[frag$side == side, ] %*% control[control$idx == idx, ]
            tg_len <- target_ov$length
            ct_len <- control_ov$length
            pval <- ks.test(x = tg_len, y = ct_len)$p.value

            return(c(idx, length(tg_len), length(ct_len), median(tg_len), median(ct_len), pval))
        } else {
            return(c(idx, NA, NA, NA, NA, NA))
        }
    }, mc.cores = cores, ignore.interactive = T) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>% 
    `colnames<-`(c("idx","npeaks","ncontrol","peak_len","control_len","pval"))
}

# Run KS test for each chromosome and side of the peak
print("Running KS test...")
right <- run_ks_test(frag_flank_gr, target_gr, control_gr, "+", ncores)
left <- run_ks_test(frag_flank_gr, target_gr, control_gr, "-", ncores)

ks_res <- full_join(right, left, by = "idx", suffix = c(".right",".left"))
ks_res$idx <- as.integer(ks_res$idx)
ks_res <- left_join(ks_res, peak[,c("chrom","start","end","idx")])

write.table(ks_res, output, row.names = F, col.names = F, quote = F, sep = '\t')

