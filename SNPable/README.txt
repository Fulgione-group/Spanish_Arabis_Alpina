# create a SNPable mask based on kmer size of 150 (read length)


Extract all overlapping k-mer subsequences as read sequences: "run_splitfasta_alpinaV5.1_150mers.sh"

Align all reads to the genome with BWA. Other aligners would work if they do global alignment w.r.t. reads and give suboptimal hits: "bwa_aln_V5.1_150mers.sh"

Convert sai to sam: "sai_to_sam_V5.1_150mers.sh"

Generate rawMask: "gen_raw_mask_V5.1_150mers.sh"

Generate final mask: "apply_mask_s_V5.1_150mers.sh"
