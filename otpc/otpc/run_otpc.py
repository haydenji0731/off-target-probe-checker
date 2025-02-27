#!/usr/bin/env python

from otpc.commons import *
from otpc import rc, detect

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
    parser.add_argument('--schema', type=lambda s: s.split(','), required=False, \
                    default=['transcript', 'ID', 'Parent', 'gene_name', 'transcript_type'])
    parser.add_argument('-o', '--out-dir', type=str, required=False, \
                    help="output directory (default: out)", default="out")
    parser.add_argument('-p', '--threads', type=int, required=False, \
                    help="number of threads (default: 1)", default=1)
    parser.add_argument('--fwd', required=False, help="", \
                    default=False, action='store_true')
    parser.add_argument('--nucmer', required=False, help="", \
                    default=False, action='store_true') # TODO: consider making this default
    parser.add_argument('-l', '--min-exact-match', required=False, help="", \
                    default=20, type=int)
    parser.add_argument('--skip-detect', required=False, help="", \
                    default=False, action='store_true')
    parser.add_argument('-b', '--binary', type=str, required=False, \
                    help="", default=None)
    parser.add_argument('-pad', '--padding', type=int, required=False, \
                    help="", default=5)
    parser.add_argument('--force', required=False, help="", \
                    default=False, action='store_true')
    parser.add_argument('--skip-index', required=False, help="", \
                    default=False, action='store_true')
    args = parser.parse_args()
    return args

def main() -> None:
    args = parse()
    if not os.path.exists(args.query) or not os.path.exists(args.target):
        print(message(f"cannot locate query and/or target files", Mtype.ERR))
        sys.exit(-1)
    if not os.path.exists(args.out_dir): os.makedirs(args.out_dir)
    # store parameters
    param_fn = os.path.join(args.out_dir, "params.json")
    store_params(args, param_fn)
    if args.fwd:
        rc.main(args)
    if not args.skip_detect:
        if args.fwd: 
            qfn = os.path.join(args.out_dir, 'fwd_oriented.fa')
            if not os.path.exists(qfn):
                print(message(f"cannot locate fwd_oriented.fa file", Mtype.ERR))
                sys.exit(-1)
        else:
            qfn = args.query
        bfn = align(qfn, args.target, "main", True, args) # second pass alignment is forward-only
        detect.main(qfn, bfn, args)

if __name__ == "__main__":
    main()