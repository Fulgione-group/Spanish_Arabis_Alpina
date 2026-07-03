library(ggplot2)
library(data.table)


###Load results
files <- sprintf("WCAN2_vs_CECAN.chr%s.xpehh.out", 1:8)

xpehh_all <- rbindlist(lapply(seq_along(files), function(i) {
  dt <- fread(files[i])
  dt[, chr := paste0("chr", i)]
  dt
}), use.names = TRUE, fill = TRUE)

# genome-wide histogram
hist(xpehh_all$xpehh, breaks = 100,
     main = "XP-EHH genome-wide",
     xlab = "XP-EHH")

# genome-wide quantiles
quantile(xpehh_all$xpehh,
         probs = c(0.99, 0.995, 0.999),
         na.rm = TRUE)

# upper 0.5% tail
cutoff <- quantile(xpehh_all$xpehh, probs = 0.995, na.rm = TRUE)
cutoff2 <- quantile(xpehh_all$xpehh, probs = 0.005, na.rm = TRUE)

tail <- xpehh_all[xpehh >= cutoff]
tail <- tail[order(chr, pos)]

# merge nearby tail SNPs within each chromosome
tail[, region_id := cumsum(c(TRUE, chr[-1] != chr[-.N] | diff(pos) > 20000))]

regions <- tail[, .(
  start = min(pos),
  end = max(pos),
  n_snps = .N,
  max_xpehh = max(xpehh),
  mean_xpehh = mean(xpehh),
  top_snp = pos[which.max(xpehh)]
), by = .(chr, region_id)]

regions[, width := end - start + 1]
regions <- regions[order(-n_snps, -max_xpehh)]

regions

head(regions)

head( xpehh_all)
xpehh_all_wcan2_vs_cecan <- xpehh_all
# chromosome offsets for Manhattan-like plot
chr_lengths <- xpehh_all[, .(chr_len = max(pos, na.rm = TRUE)), by = chr]
chr_lengths[, offset := cumsum(shift(chr_len, fill = 0))]

xpehh_all <- merge(xpehh_all, chr_lengths[, .(chr, offset)], by = "chr")
xpehh_all[, genome_pos := pos + offset]

# per-chromosome plot
ggplot(xpehh_all[seq(1, nrow(xpehh_all), by=100)], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  geom_hline(yintercept = cutoff2, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome"
  )





###Check the highest 3P-CLR peaks Table S9.
#WCAN2 peaks 

#SC5D
max(xpehh_all$xpehh[xpehh_all$chr=="chr6" & xpehh_all$pos>=28457738 & xpehh_all$pos<=28459744])
ggplot(xpehh_all[xpehh_all$chr=="chr6" & xpehh_all$pos> 28000000 & xpehh_all$pos< 29000000], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  geom_hline(yintercept = cutoff2, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome"
  )


#twd40, brrc2, vip4 
max(xpehh_all$xpehh[xpehh_all$chr=="chr8" & xpehh_all$pos>=44666284 & xpehh_all$pos<=44678737])
max(xpehh_all$xpehh[xpehh_all$chr=="chr8" & xpehh_all$pos>=46043497 & xpehh_all$pos<=46056637])
max(xpehh_all$xpehh[xpehh_all$chr=="chr8" & xpehh_all$pos>=46056760 & xpehh_all$pos<=46060745])

ggplot(xpehh_all[xpehh_all$chr=="chr8" & xpehh_all$pos> 44043497 & xpehh_all$pos< 47043497], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome")


#missing data for AA7G12830, SNDN2, DET1 p1, DET1p2, APX2, AA7G12800, LAX2, CSP4: missing data
ggplot(xpehh_all[xpehh_all$chr=="chr7" & xpehh_all$pos> 9021563 & xpehh_all$pos< 11425637], aes(x = pos, y = xpehh) ) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome")+
  geom_vline(xintercept = 10730418)+
  geom_vline(xintercept = 10408014)













  
  
  
############ WCAN1 vs WCAN2
  
files <- sprintf("WCAN1_vs_WCAN2.chr%s.xpehh.out", 1:8)

xpehh_all <- rbindlist(lapply(seq_along(files), function(i) {
  dt <- fread(files[i])
  dt[, chr := paste0("chr", i)]
  dt
}), use.names = TRUE, fill = TRUE)

# genome-wide histogram
hist(xpehh_all$xpehh, breaks = 100,
     main = "XP-EHH genome-wide",
     xlab = "XP-EHH")

# genome-wide quantiles
quantile(xpehh_all$xpehh,
         probs = c(0.99, 0.995, 0.999),
         na.rm = TRUE)

# upper 0.5% tail
cutoff <- quantile(xpehh_all$xpehh, probs = 0.995, na.rm = TRUE)
cutoff2 <- quantile(xpehh_all$xpehh, probs = 0.005, na.rm = TRUE)

tail <- xpehh_all[xpehh >= cutoff]
tail <- tail[order(chr, pos)]

# merge nearby tail SNPs within each chromosome
tail[, region_id := cumsum(c(TRUE, chr[-1] != chr[-.N] | diff(pos) > 20000))]

regions <- tail[, .(
  start = min(pos),
  end = max(pos),
  n_snps = .N,
  max_xpehh = max(xpehh),
  mean_xpehh = mean(xpehh),
  top_snp = pos[which.max(xpehh)]
), by = .(chr, region_id)]

regions[, width := end - start + 1]
regions <- regions[order(-n_snps, -max_xpehh)]

regions

head(regions)



  



###Check the highest 3P-CLR peaks Table S9.
#WCAN2 peaks 

xpehh_all$xpehh <- xpehh_all$xpehh*-1
cutoff <- quantile(xpehh_all$xpehh, probs = 0.995, na.rm = TRUE)
cutoff2 <- quantile(xpehh_all$xpehh, probs = 0.005, na.rm = TRUE)

#SC5D
max(xpehh_all$xpehh[xpehh_all$chr=="chr6" & xpehh_all$pos>=28457738 & xpehh_all$pos<=28459744])
ggplot(xpehh_all[xpehh_all$chr=="chr6" & xpehh_all$pos> 28000000 & xpehh_all$pos< 29000000], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  geom_hline(yintercept = cutoff2, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome"
  )


#twd40, brrc2, vip4 
max(xpehh_all$xpehh[xpehh_all$chr=="chr8" & xpehh_all$pos>=44666284 & xpehh_all$pos<=44678737])
max(xpehh_all$xpehh[xpehh_all$chr=="chr8" & xpehh_all$pos>=46043497 & xpehh_all$pos<=46056637])
max(xpehh_all$xpehh[xpehh_all$chr=="chr8" & xpehh_all$pos>=46056760 & xpehh_all$pos<=46060745])

ggplot(xpehh_all[xpehh_all$chr=="chr8" & xpehh_all$pos> 44043497 & xpehh_all$pos< 47043497], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome")


#missing data for AA7G12830, SNDN2, DET1 p1, DET1p2, APX2, AA7G12800, LAX2, CSP4: missing data

ggplot(xpehh_all[xpehh_all$chr=="chr7" & xpehh_all$pos> 9021563 & xpehh_all$pos< 11425637], aes(x = pos, y = xpehh) ) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome")+
  geom_vline(xintercept = 10730418)+
  geom_vline(xintercept = 10408014)



#WCAN1

###Check the highest 3P-CLR peaks Table S9.
#WCAN2 peaks 

xpehh_all$xpehh <- xpehh_all$xpehh*-1
cutoff <- quantile(xpehh_all$xpehh, probs = 0.995, na.rm = TRUE)
cutoff2 <- quantile(xpehh_all$xpehh, probs = 0.005, na.rm = TRUE)

##missing data for PCR4 and CASP ILR2 and APG1
max(xpehh_all$xpehh[xpehh_all$chr=="chr3" & xpehh_all$pos>=13795845 & xpehh_all$pos<=13799133])
max(xpehh_all$xpehh[xpehh_all$chr=="chr3" & xpehh_all$pos>=13811866 & xpehh_all$pos<=13818186])
max(xpehh_all$xpehh[xpehh_all$chr=="chr3" & xpehh_all$pos>=13816907 & xpehh_all$pos<=13818055])
max(xpehh_all$xpehh[xpehh_all$chr=="chr3" & xpehh_all$pos>=13829383 & xpehh_all$pos<=13830879])
ggplot(xpehh_all[xpehh_all$chr=="chr3" & xpehh_all$pos> 13295845 & xpehh_all$pos< 14299133], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome")+
  geom_vline(xintercept = 13795845)+
  geom_vline(xintercept = 13830879)


#AA7G10040 and CYP57 zip9 -WCAN1 peak
 
max(xpehh_all$xpehh[xpehh_all$chr=="chr7" & xpehh_all$pos>=6931385 & xpehh_all$pos<=7534874])

ggplot(xpehh_all[xpehh_all$chr=="chr7" & xpehh_all$pos> 6331385 & xpehh_all$pos< 7534874], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  geom_hline(yintercept = cutoff2, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome"
  )



########## WCAN1 vs CECAN

files <- sprintf("WCAN1_vs_CECAN.chr%s.xpehh.out", 1:8)

xpehh_all <- rbindlist(lapply(seq_along(files), function(i) {
  dt <- fread(files[i])
  dt[, chr := paste0("chr", i)]
  dt
}), use.names = TRUE, fill = TRUE)

# genome-wide histogram
hist(xpehh_all$xpehh, breaks = 100,
     main = "XP-EHH genome-wide",
     xlab = "XP-EHH")

# genome-wide quantiles
quantile(xpehh_all$xpehh,
         probs = c(0.99, 0.995, 0.999),
         na.rm = TRUE)

# upper 0.5% tail
cutoff <- quantile(xpehh_all$xpehh, probs = 0.995, na.rm = TRUE)
cutoff2 <- quantile(xpehh_all$xpehh, probs = 0.005, na.rm = TRUE)

tail <- xpehh_all[xpehh >= cutoff]
tail <- tail[order(chr, pos)]

# merge nearby tail SNPs within each chromosome
tail[, region_id := cumsum(c(TRUE, chr[-1] != chr[-.N] | diff(pos) > 20000))]

regions <- tail[, .(
  start = min(pos),
  end = max(pos),
  n_snps = .N,
  max_xpehh = max(xpehh),
  mean_xpehh = mean(xpehh),
  top_snp = pos[which.max(xpehh)]
), by = .(chr, region_id)]

regions[, width := end - start + 1]
regions <- regions[order(-n_snps, -max_xpehh)]

regions

head(regions)

##missing data for PCR4 and CASP ILR2 and APG1
max(xpehh_all$xpehh[xpehh_all$chr=="chr3" & xpehh_all$pos>=13795845 & xpehh_all$pos<=13799133])
max(xpehh_all$xpehh[xpehh_all$chr=="chr3" & xpehh_all$pos>=13811866 & xpehh_all$pos<=13818186])
max(xpehh_all$xpehh[xpehh_all$chr=="chr3" & xpehh_all$pos>=13816907 & xpehh_all$pos<=13818055])
max(xpehh_all$xpehh[xpehh_all$chr=="chr3" & xpehh_all$pos>=13829383 & xpehh_all$pos<=13830879])
ggplot(xpehh_all[xpehh_all$chr=="chr3" & xpehh_all$pos> 13295845 & xpehh_all$pos< 14299133], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome")+
  geom_vline(xintercept = 13795845)+
  geom_vline(xintercept = 13830879)


#AA7G10040 and CYP57 zip9 -WCAN1 peak

max(xpehh_all$xpehh[xpehh_all$chr=="chr7" & xpehh_all$pos>=6945018 & xpehh_all$pos<=6948658])
max(xpehh_all$xpehh[xpehh_all$chr=="chr7" & xpehh_all$pos>=6931385 & xpehh_all$pos<=6934874])

ggplot(xpehh_all[xpehh_all$chr=="chr7" & xpehh_all$pos> 6331385 & xpehh_all$pos< 7534874], aes(x = pos, y = xpehh)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  geom_hline(yintercept = cutoff2, linetype = "dashed") +
  facet_wrap(~ chr, scales = "free_x", ncol = 2) +
  theme_bw() +
  labs(
    x = "Position (bp)",
    y = "XP-EHH",
    title = "WCAN2 vs CECAN XP-EHH per chromosome"
  )+
  geom_segment(x=6934874, xend=6931385 , y =-1.4 ,yend=-1.4, size=2)+
  geom_segment(x=6931385, xend=6948658 , y =-1.4 ,yend=-1.4, size=10)+
  geom_segment(x=6945018, xend=6948658 , y =-1.4 ,yend=-1.4, size=10)+
  geom_vline(xintercept = 6931385)
