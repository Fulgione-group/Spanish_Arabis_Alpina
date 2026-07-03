#!/usr/bin/env python3

import tskit
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
        if line.strip():
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


def root_stats(tree):
    children = tree.children(tree.root)

    if len(children) != 2:
        return [np.nan] * 10

    c1, c2 = children

    n1 = tree.num_samples(c1)
    n2 = tree.num_samples(c2)

    small = min(n1, n2)
    large = max(n1, n2)
    balance = small / (small + large)

    root_time = tree.time(tree.root)

    child1_time = tree.time(c1)
    child2_time = tree.time(c2)

    root_child_branch1 = root_time - child1_time
    root_child_branch2 = root_time - child2_time

    mean_root_child_branch = (
        root_child_branch1 + root_child_branch2
    ) / 2

    if root_time > 0:
        child_root_ratio = mean_root_child_branch / root_time
    else:
        child_root_ratio = np.nan

    return (
        small,
        large,
        balance,
        root_time,
        child1_time,
        child2_time,
        root_child_branch1,
        root_child_branch2,
        mean_root_child_branch,
        child_root_ratio,
    )


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


outfile = "All_trees_summary_no_lowdensity_with_root_child_ratio.tsv"

with open(outfile, "w") as output_file:

    output_file.write(
        "chromosome\ttreeID\ttreeStart\ttreeEnd\tspan_bp\t"
        "numMutations\t"
        "timeRoot\tchild1Time\tchild2Time\t"
        "rootChildBranch1\trootChildBranch2\t"
        "meanRootChildBranch\tchildRootRatio\t"
        "rootSmallSide\trootLargeSide\trootBalance\t"
        "CAN_FR_concordance\tCAN_FR_misclassified\t"
        "collessIndex\n"
    )

    for chromosome, tree_seq_file in tree_seq_files.items():

        print("Loading", tree_seq_file)

        ts = tskit.load(tree_seq_file)

        if ts.num_samples != N:
            raise ValueError(
                f"{chromosome}: tree samples = {ts.num_samples}, group samples = {N}"
            )

        for tree in ts.trees():

            tree_start, tree_end = tree.interval
            span = tree_end - tree_start

            if overlaps_gap(chromosome, tree_start, tree_end):
                continue

            (
                small,
                large,
                balance,
                root_time,
                child1_time,
                child2_time,
                root_child_branch1,
                root_child_branch2,
                mean_root_child_branch,
                child_root_ratio,
            ) = root_stats(tree)

            concordance, misclassified = can_fr_concordance(tree)
            coll = colless_index(tree)

            output_file.write(
                f"{chromosome}\t{tree.index}\t"
                f"{tree_start:.1f}\t{tree_end:.1f}\t{span:.1f}\t"
                f"{tree.num_mutations}\t"
                f"{root_time:.6f}\t{child1_time:.6f}\t{child2_time:.6f}\t"
                f"{root_child_branch1:.6f}\t{root_child_branch2:.6f}\t"
                f"{mean_root_child_branch:.6f}\t{child_root_ratio:.6f}\t"
                f"{small}\t{large}\t{balance:.6f}\t"
                f"{concordance:.6f}\t{misclassified}\t"
                f"{coll}\n"
            )

print("Wrote", outfile)
