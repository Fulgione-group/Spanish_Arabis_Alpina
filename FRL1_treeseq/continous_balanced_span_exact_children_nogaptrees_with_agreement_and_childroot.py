import glob
import tskit
from collections import defaultdict

cutoff = 152
min_agreement = 0.95
bedfile = "low_density_gaps_all.bed"

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


def chrom_from_treefile(f):
    return f.replace("relate_popsize_CAN_FR_msp_", "").replace(".trees", "")


def overlaps_gap(chrom, start, end):
    for gstart, gend in gaps.get(chrom, []):
        if gstart < end and gend > start:
            return True
    return False


def split_agreement(split1, split2):
    a1, b1 = split1
    a2, b2 = split2

    n = len(a1) + len(b1)

    same = len(a1 & a2) + len(b1 & b2)
    swapped = len(a1 & b2) + len(b1 & a2)

    return max(same, swapped) / n


all_runs = []

for f in sorted(glob.glob("relate_popsize_CAN_FR_msp_chr*.trees")):

    print("Loading", f)

    chrom = chrom_from_treefile(f)
    ts = tskit.load(f)

    run_start = None
    run_end = None
    run_tree_index = None
    prev_split = None
    n_trees = 0
    min_agree_run = 1.0

    min_root_age = None
    max_root_age = None
    sum_root_age = 0.0

    min_root_child_branch = None
    max_root_child_branch = None
    sum_root_child_branch = 0.0

    for tree in ts.trees():

        tstart = tree.interval.left
        tend = tree.interval.right

        if overlaps_gap(chrom, tstart, tend):
            ok = False

        else:
            children = tree.children(tree.root)

            if len(children) != 2:
                ok = False

            else:
                s1 = frozenset(tree.samples(children[0]))
                s2 = frozenset(tree.samples(children[1]))

                if min(len(s1), len(s2)) >= cutoff:
                    split = (s1, s2)

                    root_age = tree.time(tree.root)

                    child_ages = [
                        tree.time(children[0]),
                        tree.time(children[1])
                    ]

                    root_child_branches = [
                        root_age - child_ages[0],
                        root_age - child_ages[1]
                    ]

                    mean_root_child_branch = sum(root_child_branches) / 2

                    ok = True
                else:
                    ok = False

        if ok:

            if run_start is None:
                run_start = tstart
                run_end = tend
                run_tree_index = tree.index
                prev_split = split
                n_trees = 1
                min_agree_run = 1.0

                min_root_age = root_age
                max_root_age = root_age
                sum_root_age = root_age

                min_root_child_branch = mean_root_child_branch
                max_root_child_branch = mean_root_child_branch
                sum_root_child_branch = mean_root_child_branch

            else:
                agree = split_agreement(prev_split, split)

                if agree >= min_agreement:
                    run_end = tend
                    n_trees += 1
                    prev_split = split
                    min_agree_run = min(min_agree_run, agree)

                    min_root_age = min(min_root_age, root_age)
                    max_root_age = max(max_root_age, root_age)
                    sum_root_age += root_age

                    min_root_child_branch = min(
                        min_root_child_branch,
                        mean_root_child_branch
                    )
                    max_root_child_branch = max(
                        max_root_child_branch,
                        mean_root_child_branch
                    )
                    sum_root_child_branch += mean_root_child_branch

                else:
                    all_runs.append(
                        (
                            f,
                            run_start,
                            run_end,
                            run_end - run_start,
                            n_trees,
                            run_tree_index,
                            round(min_agree_run, 4),
                            min_root_age,
                            max_root_age,
                            sum_root_age / n_trees,
                            min_root_child_branch,
                            max_root_child_branch,
                            sum_root_child_branch / n_trees
                        )
                    )

                    run_start = tstart
                    run_end = tend
                    run_tree_index = tree.index
                    prev_split = split
                    n_trees = 1
                    min_agree_run = 1.0

                    min_root_age = root_age
                    max_root_age = root_age
                    sum_root_age = root_age

                    min_root_child_branch = mean_root_child_branch
                    max_root_child_branch = mean_root_child_branch
                    sum_root_child_branch = mean_root_child_branch

        else:

            if run_start is not None:
                all_runs.append(
                    (
                        f,
                        run_start,
                        run_end,
                        run_end - run_start,
                        n_trees,
                        run_tree_index,
                        round(min_agree_run, 4),
                        min_root_age,
                        max_root_age,
                        sum_root_age / n_trees,
                        min_root_child_branch,
                        max_root_child_branch,
                        sum_root_child_branch / n_trees
                    )
                )

            run_start = None
            run_end = None
            run_tree_index = None
            prev_split = None
            n_trees = 0
            min_agree_run = 1.0

            min_root_age = None
            max_root_age = None
            sum_root_age = 0.0

            min_root_child_branch = None
            max_root_child_branch = None
            sum_root_child_branch = 0.0

    if run_start is not None:
        all_runs.append(
            (
                f,
                run_start,
                run_end,
                run_end - run_start,
                n_trees,
                run_tree_index,
                round(min_agree_run, 4),
                min_root_age,
                max_root_age,
                sum_root_age / n_trees,
                min_root_child_branch,
                max_root_child_branch,
                sum_root_child_branch / n_trees
            )
        )

all_runs = sorted(all_runs, key=lambda x: x[3], reverse=True)

outfile = "balanced_rootsplit_runs_95agreement_no_lowdensity_with_root_and_branch_ages.tsv"

with open(outfile, "w") as out:
    out.write(
        "chrom\tstart\tend\tspan_bp\tn_trees\tfirst_tree\t"
        "min_agreement\t"
        "min_root_age\tmax_root_age\tmean_root_age\t"
        "min_mean_root_child_branch\tmax_mean_root_child_branch\t"
        "mean_root_child_branch\n"
    )

    for r in all_runs:
        out.write(
            f"{r[0]}\t{r[1]:.1f}\t{r[2]:.1f}\t{r[3]:.1f}\t"
            f"{r[4]}\t{r[5]}\t{r[6]}\t"
            f"{r[7]}\t{r[8]}\t{r[9]}\t"
            f"{r[10]}\t{r[11]}\t{r[12]}\n"
        )

print("Wrote", outfile)
print("Found", len(all_runs), "runs")
