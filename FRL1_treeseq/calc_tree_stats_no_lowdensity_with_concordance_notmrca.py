import tskit
import itertools
import numpy as np
from collections import defaultdict

bedfile = "low_density_gaps_all.bed"
groupfile = "CAN_FR_FRIL1_no881.poplabels"

tree_seq_files = {
    "chr1": "relate_popsize_CAN_FR_msp_chr1.trees",
    "chr2": "relate_popsize_CAN_FR_msp_chr2.trees",
    "chr3": "relate_popsize_CAN_FR_msp_chr3.trees",
    "chr4": "relate_popsize_CAN_FR_msp_chr4.trees",
    "chr5": "relate_popsize_CAN_FR_msp_chr5.trees",
    "chr6": "relate_popsize_CAN_FR_msp_chr6.trees",
    "chr7": "relate_popsize_CAN_FR_msp_chr7.trees",
    "chr8": "relate_popsize_CAN_FR_msp_chr8.trees",
}

# load sample groups
groups = []

with open(groupfile) as f:
    header = next(f)
    for line in f:
        if line.strip() == "":
            continue
        fields = line.split()
        groups.append(fields[2])  # CAN / FR

CAN = {i for i, g in enumerate(groups) if g == "CAN"}
FR = {i for i, g in enumerate(groups) if g == "FR"}
N = len(CAN) + len(FR)

print("Loaded groups:", len(CAN), "CAN,", len(FR), "FR")

# load low-density gaps
gaps = defaultdict(list)

with open(bedfile) as bed:
    for line in bed:
        if line.strip() == "" or line.startswith("#"):
            continue
        fields = line.split()
        chrom = fields[0]
        start = float(fields[1])
        end = float(fields[2])
        gaps[chrom].append((start, end))

for chrom in gaps:
    gaps[chrom].sort()


def overlaps_gap(chrom, start, end):
    for gstart, gend in gaps.get(chrom, []):
        if gstart < end and gend > start:
            return True
    return False


def colless_index(tree):
    imbalance = 0
    for node in tree.nodes():
        children = list(tree.children(node))
        if len(children) == 2:
            n1 = tree.num_samples(children[0])
            n2 = tree.num_samples(children[1])
            imbalance += abs(n1 - n2)
    return imbalance


def root_balance(tree):
    children = tree.children(tree.root)

    if len(children) != 2:
        return np.nan, np.nan, np.nan

    sizes = [tree.num_samples(c) for c in children]
    small = min(sizes)
    large = max(sizes)
    balance = small / (small + large)

    return small, large, balance


def can_fr_concordance(tree):
    children = tree.children(tree.root)

    if len(children) != 2:
        return np.nan, np.nan

    A = set(tree.samples(children[0]))
    B = set(tree.samples(children[1]))

    option1 = len(A & CAN) + len(B & FR)
    option2 = len(A & FR) + len(B & CAN)

    same = max(option1, option2)

    concordance = same / N
    misclassified = N - same

    return concordance, misclassified


outfile = "All_trees_summary_2_no_lowdensity_with_CANFR_concordance.tsv"

with open(outfile, "w") as output_file:

    output_file.write(
        "chromosome\ttreeID\ttreeStart\ttreeEnd\t"
        "numMutations\ttimeRoot\t"
        "rootSmallSide\trootLargeSide\trootBalance\t"
        "CAN_FR_concordance\tCAN_FR_misclassified\t"
        "collessIndex\tmeanTMRCA\tmedianTMRCA\tmaxTMRCA\tvarTMRCA\n"
    )

    for chromosome, tree_seq_file in tree_seq_files.items():

        print("Loading", tree_seq_file)

        ts = tskit.load(tree_seq_file)
        samples = ts.samples()
        tree_id = 0

        if len(samples) != N:
            raise ValueError(
                f"Number of tree samples ({len(samples)}) does not match group file ({N})"
            )

        for tree in ts.trees():

            tree_start, tree_end = tree.interval

            if overlaps_gap(chromosome, tree_start, tree_end):
                tree_id += 1
                continue

            root_time = tree.time(tree.root)
            num_mutations = tree.num_mutations

            small, large, balance = root_balance(tree)
            concordance, misclassified = can_fr_concordance(tree)

            tmrcas = []

            for s1, s2 in itertools.combinations(samples, 2):
                mrca = tree.mrca(s1, s2)
                tmrcas.append(tree.time(mrca))

            mean_tmrca = np.mean(tmrcas)
            median_tmrca = np.median(tmrcas)
            max_tmrca = np.max(tmrcas)
            var_tmrca = np.var(tmrcas)

            coll = colless_index(tree)

            output_file.write(
                f"{chromosome}\t{tree_id}\t{tree_start}\t{tree_end}\t"
                f"{num_mutations}\t{root_time:.2f}\t"
                f"{small}\t{large}\t{balance:.4f}\t"
                f"{concordance:.4f}\t{misclassified}\t"
                f"{coll}\t{mean_tmrca:.2f}\t{median_tmrca:.2f}\t"
                f"{max_tmrca:.2f}\t{var_tmrca:.2f}\n"
            )

            tree_id += 1

print("Wrote", outfile)
