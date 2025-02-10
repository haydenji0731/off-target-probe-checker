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