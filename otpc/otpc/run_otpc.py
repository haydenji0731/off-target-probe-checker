#!/usr/bin/env python

from otpc.commons import *
from otpc import rc

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--version', action='version', version='%(prog)s v0.0.1')
    parser.add_argument('-q', '--query', type=str, required=True, \
                    help="query probe sequences (fasta)")
    parser.add_argument('-t', '--target', type=str, required=True, \
                    help="target transcript sequences (fasta)")
    parser.add_argument('-a', '--annotation', type=str, required=True, \
                    help="target transcriptome annotation (gff or gtf)")
    parser.add_argument('--gtf', required=False, default=False, help="", \
                    action='store_true')
    parser.add_argument('--bam', required=False, default=False, help="", \
                    action='store_true')
    parser.add_argument('--keep-dot', required=False, default=False, help="", \
                    action='store_true')
    parser.add_argument('--schema', type=list, required=False, \
                    help="", default=['transcript', 'ID', 'Parent'])
    parser.add_argument('-o', '--out-dir', type=str, required=False, \
                    help="output directory (default: out)", default="out")
    parser.add_argument('-p', '--threads', type=int, required=False, \
                    help="number of threads (default: 1)", default=1)
    parser.add_argument('--fwd', required=False, help="", \
                    default=False, action='store_true')
    parser.add_argument('--nucmer', required=False, help="", \
                    default=False, action='store_true') # TODO: consider making this default
    parser.add_argument('-b', '--binary', type=str, required=False, \
                    help="", default=None)
    args = parser.parse_args()
    return args

def main() -> None:
    args = parse()
    if not os.path.exists(args.out_dir): os.makedirs(args.out_dir)
    if args.fwd:
        rc.main(args)
        new_qfn = os.path.join(args.out_dir, 'fwd_oriented.fa')
        align(new_qfn, args.target, 'aln', args)
    else:
        align(args.query, args.target, 'aln', args)
    

if __name__ == "__main__":
    main()