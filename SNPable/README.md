# Create a SNPable mask based on k-mer size of 150 (read length)

Furthermore, we retained only uniquely mapping genomic regions inferred with [SNPable](https://lh3lh3.users.sourceforge.net/snpable.shtml),  
using a k-mer length of 150 and intermediate stringency (`r=0.5`).

- Extract all overlapping k-mer subsequences as read sequences:  
  `run_splitfasta_alpinaV5.1_150mers.sh`

- Align all reads to the genome with BWA. Other aligners would also work if they perform global alignment with respect to reads and report suboptimal hits:  
  `bwa_aln_V5.1_150mers.sh`

- Convert `.sai` to `.sam`:  
  `sai_to_sam_V5.1_150mers.sh`

- Generate raw mask:  
  `gen_raw_mask_V5.1_150mers.sh`

- Generate final mask:  
  `apply_mask_s_V5.1_150mers.sh`
