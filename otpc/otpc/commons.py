import argparse
import os
import pyfastx
import pysam
from subprocess import call
import pandas as pd
from datetime import datetime
from enum import Enum
import sys
from Bio.Seq import Seq

RED = '\033[31m'
GREEN = '\033[32m'
RESET = '\033[0m'

class Mtype(Enum):
    PROG = (GREEN, "PROGRESS")
    ERR = (RED, "ERROR")
    WARN = (RED, "WARNING")

def message(s, mtype) -> str:
    if mtype not in Mtype:
        raise Exception("Error while printing message")
    return f"{datetime.now()} {mtype.value[0]}{mtype.value[1]}{RESET} {s}"

def align(qfn, tfn, prefix, args) -> str:
    print(message(f"aligning query probes to target transcripts", Mtype.PROG))
    ofn = os.path.join(args.out_dir, f'{prefix}.bam' if args.bam else f'{prefix}.sam')
    if args.binary:
        aligner = args.binary
    else:
        aligner = "nucmer" if args.nucmer else "bowtie2"
    if args.nucmer: # nucmer flow
        raise NotImplementedError # TODO: implement
    else: # bt2 flow
        idx_fn = os.path.join(args.out_dir, 'target')
        cmd = f'{aligner}-build -q {tfn} {idx_fn} --threads {args.threads}'
        print(cmd)
        call(cmd, shell=True)
        if args.bam:
            cmd = f'{aligner} -f -a -N 1 --local -x {idx_fn} ' + \
                f'-U {qfn} --very-sensitive-local --threads {args.threads} ' + \
                f'| samtools view -b -o {ofn} -@ {args.threads}'
        else:
            cmd = f'{aligner} -f -a -N 1 --local -x {idx_fn} ' + \
                f'-U {qfn} --very-sensitive-local --threads {args.threads} -S {ofn}'
        print(cmd)
        call(cmd, shell=True)
    return ofn