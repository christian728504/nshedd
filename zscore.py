#!/usr/bin/env -S uv run --script
#
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "joblib>=1.5.3",
#     "numpy>=2.5.1",
#     "pybigwig>=0.3.25",
#     "pysam>=0.24.0",
# ]
# ///

"""
This is a uv-runnable command-line script that computes per-region z-scores of signal enrichment over a set of genomic intervals (a BED file).
Given a signal file plus peaks.bed, it reads the average/count signal in each region, either from a BigWig file (zscores → mean signal per region
via pyBigWig) or from a BAM file (bamzscores → counts reads whose 5'/3' position falls inside each region via pysam), parallelizing the
work across regions with joblib. It then log-transforms and standardizes those values into z-scores (zt/ztpm), assigning a floor of -10 to regions
with zero signal, and prints each BED line with its z-score appended, optionally filtering to regions above a threshold (or printing all when the
threshold is "NA"). The BAM path also normalizes to a TPM-like scale using total mapped reads from idxstats and can write out the total
in-peak read count.
"""

import sys

sys.dont_write_bytecode = True

import argparse
import os
import gzip
import statistics
import pyBigWig
import pysam
import math
import numpy

from joblib import Parallel, delayed


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx : min(ndx + n, l)]


def flatten(x):
    r = []
    for xx in x:
        r += xx
    return r


def read(bigwig, regions):
    with pyBigWig.open(bigwig) as bw:
        return [
            statistics.mean(
                bw.values(x.split()[0], int(x.split()[1]), int(x.strip().split()[2]))
            )
            for x in f
        ]


def tmean(bw, *args):
    try:
        return statistics.mean(bw.values(*args))
    except:
        return 0.0


def contains(read, start, end):
    position = read.reference_start if not read.is_reverse else read.reference_end
    return position is not None and position > start and position < end


def bamcount(chromosome, start, end, bam, total, i=None):
    if i is not None:
        print(i, file=sys.stderr)
    return len(
        [x for x in bam.fetch(chromosome, start, end) if contains(x, start, end)]
    )


def bamzscoress(regions, bam, total):
    b = pysam.AlignmentFile(bam, "rb")
    r = [
        bamcount(
            x.split()[0],
            int(x.split()[1]),
            int(x.strip().split()[2]),
            b,
            total,
            "%d / %d" % (i, len(regions)) if i % 1000 == 0 else None,
        )
        for i, x in enumerate(regions)
    ]
    b.close()
    return r


def bamzscores(bed, bam, j=1):
    with (gzip.open if bed.endswith(".gz") else open)(bed, "rt") as f:
        lines = [x.strip() for x in f]
    total = sum(
        [
            int(l.strip().split("\t")[2])
            for l in pysam.idxstats(bam)
            if len(l.strip().split("\t")) >= 3
        ]
    )
    bsetlen = int(math.ceil(len(lines) / float(j)))
    breakpoints = [bsetlen * x for x in range(j + 1)]
    jsets = [
        lines[breakpoints[i] : breakpoints[i + 1]] for i in range(len(breakpoints) - 1)
    ]
    jresults = Parallel(n_jobs=j)(delayed(bamzscoress)(x, bam, total) for x in jsets)
    allresults = []
    for x in jresults:
        allresults += x
    return ztpm(allresults), sum(allresults)


def ztpm(signal):
    treads = sum(signal)
    zeros = len([x for x in signal if x == 0])
    r = [x * 100000 / float(treads + zeros) for x in signal]
    m = numpy.mean([math.log(x) for x in r if x > 0])
    s = numpy.std([math.log(x) for x in r if x > 0])
    return [(math.log(x) - m) / s if x > 0 else -10 for x in r]


def zt(signal):
    signalmean = statistics.mean([math.log(x) for x in signal if x > 0.0])
    signalstd = statistics.stdev([math.log(x) for x in signal if x > 0.0])
    return [
        (math.log(x) - signalmean) / signalstd if x > 0.0 else -10.0 for x in signal
    ]


def tmeans(bigwig, group):
    bw = pyBigWig.open(bigwig)
    if bw is None:
        raise Exception("Error opening %s: no such file or directory." % bigwig)
    return [
        tmean(bw, x.split()[0], int(x.split()[1]), int(x.strip().split()[2]))
        for x in group
    ]


def zscores(bed, bigwig, j=1):
    with (gzip.open if bed.endswith(".gz") else open)(bed, "rt") as f:
        l = [x for x in f]
    signal = flatten(
        Parallel(n_jobs=48)(delayed(tmeans)(bigwig, x) for x in batch(l, 48))
    )
    return zt(signal), None


def assert_lt_nproc(value) -> int:
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "value passed to jobs could not be parsed into an int"
        )
    if value > os.cpu_count():
        raise argparse.ArgumentTypeError(
            "number of jobs (%d) exceeds CPU count (%d)" % (value, os.cpu_count())
        )
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("signal", type=str, help="BAM or BigWig file")
    parser.add_argument("regions", type=str, help="BED file of genomic intervals")
    parser.add_argument(
        "--threshold", "-t", default="NA", type=float, help="zscore filter threshold"
    )
    parser.add_argument(
        "--jobs",
        "-j",
        type=assert_lt_nproc,
        default=1,
        help="Number of jobs to parallelize",
    )
    parser.add_argument(
        "--reads-in-peaks",
        default=None,
        help="For BAM input, file to write the total in-peak read count. If it "
        "already holds per-region scores (last column), those values replace the "
        "z-scores, matched by rank.",
    )
    args = parser.parse_args()

    scoring_fn = bamzscores if args.signal.endswith(".bam") else zscores
    z, t = scoring_fn(args.regions, args.signal, args.jobs)

    if args.reads_in_peaks:
        with open(args.reads_in_peaks, "r") as f:
            sortedscores = sorted([float(x.strip().split()[-1]) for x in f])
        sortedindexes = sorted(range(len(z)), key=lambda i: z[i])
        sortedindexes = {v: k for k, v in enumerate(sortedindexes)}
    else:
        sortedscores = z
        sortedindexes = range(len(z))
    with (gzip.open if args.regions.endswith(".gz") else open)(args.regions, "rt") as f:
        for idx, line in enumerate(f):
            if sortedscores[sortedindexes[idx]] > args.threshold:
                print("%s\t%f" % (line.strip(), sortedscores[sortedindexes[idx]]))
    if args.signal.endswith(".bam") and args.reads_in_peaks:
        with open(args.reads_in_peaks, "w") as o:
            o.write(str(t) + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
