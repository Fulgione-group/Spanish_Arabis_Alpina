#!/usr/bin/env python3
"""
Simulate with msprime using your hard-coded demography (piecewise-exponential Ne),
then write a VCF.

Requires: msprime>=1.2, tskit
"""

import math
import argparse
import msprime
import pandas as pd
import demes

# ----------------------------
# global parameters (defaults)
# ----------------------------
MU = 7.0e-09
RECOMB_RATE = 1.096305050595e-08  # per bp per generation
SEQUENCE_LENGTH = 50_000_000

# defaults (will be overridden by CLI args if provided)
RANDOM_SEED_ANC = 1
RANDOM_SEED_MUT = 2
OUT_VCF = "simulated_3pclr_v3.vcf"
OUT_DEMES_YML = "demography.yml"

# ----------------------------
# split times in generations ago (from your file)
# ----------------------------
T_WCAN = 19547.47        # WCAN2_vs_WCAN1
T_CAN_ES17 = 45289.14    # CAN_vs_ES17
T_CAN_CECAN = 82204.56   # CAN_vs_CECAN
T_ANC_FR = 151050.3      # Anc_vs_FR01


# ----------------------------
# Ne(t) points (time in generations ago; Ne at that time)
# ----------------------------
FR_times = [0, 666.667, 926.33, 1287.13, 1788.46, 2485.06, 3452.98, 4797.91, 6666.67, 9263.3, 12871.3, 17884.6, 24850.6, 34529.8, 47979.1, 66666.7, 92633, 128713]
FR_Ne    = [24899.530395, 10021.847628, 8269.889498, 7496.049582, 7016.490155, 6753.316216, 7271.75418, 9886.835283, 16844.545063, 29089.555104, 49528.979406, 78086.322868, 103291.913276, 132087.135241, 141907.175678, 120468.864822, 90316.45078, 70435.828734]

CECAN_times = [0, 666.667, 926.33, 1287.13, 1788.46, 2485.06, 3452.98, 4797.91, 6666.67, 9263.3, 12871.3, 17884.6, 24850.6, 34529.8, 47979.1, 66666.7]
CECAN_Ne    = [3931.5285, 1917.207319, 2114.218543, 2792.344508, 4022.00843, 5656.319749, 8136.86864, 12815.353819, 21714.206302, 32069.372467, 41044.16352, 48634.820586, 52054.161314, 68503.942402, 81057.377275, 83883.328356]

ES17_times = [0, 666.667, 926.33, 1287.13, 1788.46, 2485.06, 3452.98, 4797.91, 6666.67, 9263.3, 12871.3, 17884.6, 24850.6, 34529.8]
ES17_Ne    = [6362.53738, 2906.588073, 2967.02449, 3669.751705, 4786.017172, 6075.459639, 8241.699372, 12991.298428, 23988.638981, 42764.649031, 67441.888697, 97157.936055, 113960.373699, 130089.085005]

WCAN1_times = [0, 666.667, 926.33, 1287.13, 1788.46, 2485.06, 3452.98, 4797.91, 6666.67, 9263.3, 12871.3, 17884.6]
WCAN1_Ne    = [18036.216723, 9189.521005, 9693.867659, 11523.525277, 13828.083732, 16183.743753, 19172.514284, 24698.920163, 36251.323173, 49693.885664, 62113.112948, 70397.945506]

WCAN2_times = [0, 666.667, 926.33, 1287.13, 1788.46, 2485.06, 3452.98, 4797.91, 6666.67, 9263.3]
WCAN2_Ne    = [14106.879361, 7419.861783, 8135.61087, 10508.506636, 14124.453384, 18775.605232, 24875.374374, 37603.597912, 61031.133202, 91616.705024]

# internal ancestors (demes-safe: start at their split time, not 0)
WCAN_ANC_times = [T_WCAN, 24850.6, 34529.8]   # NO T_CAN_ES17 here
WCAN_ANC_Ne    = [113449.155124, 91453.019177, 99117.855733]

CAN1_ANC_times = [T_CAN_ES17, 47979.1, 66666.7]   # NO T_CAN_CECAN here
CAN1_ANC_Ne    = [109836.964151, 112693.3254, 96712.547099]

CAN_ANC_times = [T_CAN_CECAN, 92633, 128713]   # no T_ANC_FR here
CAN_ANC_Ne    = [81818.5170234, 73131.453005, 49028.201533]

ROOT_times = [T_ANC_FR, 178846, 248506, 345298, 479791, 666667]
ROOT_Ne    = [47982.9106066, 46713.255353, 90234.30239, 146677.031844, 160208.655753, 205947.771645]


def add_piecewise_exponential_Ne(demography, pop, times, Nes):
    if times[0] < 0:
        raise ValueError(f"{pop}: times must be >= 0")
    if any(t2 <= t1 for t1, t2 in zip(times, times[1:])):
        raise ValueError(f"{pop}: times must be strictly increasing")

    demography.add_population(name=pop, initial_size=float(Nes[0]), growth_rate=0.0)

    for i in range(len(times) - 1):
        t0, t1 = float(times[i]), float(times[i + 1])
        Ne0, Ne1 = float(Nes[i]), float(Nes[i + 1])

        g = -math.log(Ne1 / Ne0) / (t1 - t0)

        demography.add_population_parameters_change(time=t0, population=pop, growth_rate=g)
        demography.add_population_parameters_change(
            time=t1, population=pop, growth_rate=0.0, initial_size=Ne1
        )


def make_individual_names(sample_sets):
    """
    One name per diploid individual (not per haplotype).
    Length must equal total individuals.
    """
    names = []
    for ss in sample_sets:
        pop = ss.population
        for i in range(ss.num_samples):
            names.append(f"{pop}_{i}")
    return names


def parse_args():
    parser = argparse.ArgumentParser(
        description="msprime simulation with variable VCF output and replicate-specific seeds"
    )

    parser.add_argument(
        "--out-vcf",
        default=OUT_VCF,
        help="Output VCF filename (default: %(default)s)"
    )

    parser.add_argument(
        "--rep",
        type=int,
        default=1,
        help="Replicate ID (used to offset seeds; default: %(default)s)"
    )

    parser.add_argument(
        "--seed-anc",
        type=int,
        default=None,
        help="Override ancestry random seed (default: derived from --rep)"
    )

    parser.add_argument(
        "--seed-mut",
        type=int,
        default=None,
        help="Override mutation random seed (default: derived from --rep)"
    )

    parser.add_argument(
        "--out-demes-yml",
        default=OUT_DEMES_YML,
        help="Where to write the demes YAML (default: %(default)s)"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # replicate-specific, reproducible seeds (unless explicitly overridden)
    global RANDOM_SEED_ANC, RANDOM_SEED_MUT, OUT_VCF
    OUT_VCF = args.out_vcf
    RANDOM_SEED_ANC = args.seed_anc if args.seed_anc is not None else (1 + args.rep)
    RANDOM_SEED_MUT = args.seed_mut if args.seed_mut is not None else (10_000 + args.rep)

    demography = msprime.Demography()

    # extant
    add_piecewise_exponential_Ne(demography, "FR", FR_times, FR_Ne)
    add_piecewise_exponential_Ne(demography, "CECAN", CECAN_times, CECAN_Ne)
    add_piecewise_exponential_Ne(demography, "ES17", ES17_times, ES17_Ne)
    add_piecewise_exponential_Ne(demography, "WCAN1", WCAN1_times, WCAN1_Ne)
    add_piecewise_exponential_Ne(demography, "WCAN2", WCAN2_times, WCAN2_Ne)

    # internal ancestors (exist at time 0; not sampled)
    add_piecewise_exponential_Ne(demography, "WCAN_ANC", WCAN_ANC_times, WCAN_ANC_Ne)
    add_piecewise_exponential_Ne(demography, "CAN1_ANC", CAN1_ANC_times, CAN1_ANC_Ne)
    add_piecewise_exponential_Ne(demography, "CAN_ANC", CAN_ANC_times, CAN_ANC_Ne)
    add_piecewise_exponential_Ne(demography, "ROOT", ROOT_times, ROOT_Ne)

    # splits
    demography.add_population_split(time=T_WCAN,      derived=["WCAN1", "WCAN2"], ancestral="WCAN_ANC")
    demography.add_population_split(time=T_CAN_ES17,  derived=["ES17", "WCAN_ANC"], ancestral="CAN1_ANC")
    demography.add_population_split(time=T_CAN_CECAN, derived=["CECAN", "CAN1_ANC"], ancestral="CAN_ANC")
    demography.add_population_split(time=T_ANC_FR,    derived=["FR", "CAN_ANC"], ancestral="ROOT")

    demography.sort_events()

    # sample sizes from your file
    samples = [
        msprime.SampleSet(215, population="FR", time=0),
        msprime.SampleSet(25,  population="CECAN", time=0),
        msprime.SampleSet(15,  population="ES17", time=0),
        msprime.SampleSet(97,  population="WCAN1", time=0),
        msprime.SampleSet(74,  population="WCAN2", time=0),
    ]
    sample_names = make_individual_names(samples)

    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=SEQUENCE_LENGTH,
        recombination_rate=RECOMB_RATE,
        ploidy=2,
        random_seed=RANDOM_SEED_ANC,
    )
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=RANDOM_SEED_MUT)

    # --- Write VCF ---
    with open(OUT_VCF, "w") as f:
        ts.write_vcf(f, contig_id="chr1", individual_names=sample_names)

    print(f"Wrote VCF: {OUT_VCF}")
    print(f"Seeds: anc={RANDOM_SEED_ANC} mut={RANDOM_SEED_MUT} rep={args.rep}")
    print(ts)
    print("demorgaphy summary")
    print(demography)
    dbg = demography.debug()
    print("dem_summary")
    print(dbg)

    graph = demography.to_demes()
    yml_txt = demes.dumps(graph)
    print(yml_txt)  # YAML text to stdout
    with open(args.out_demes_yml, "w") as f:
        f.write(yml_txt)


if __name__ == "__main__":
    main()
