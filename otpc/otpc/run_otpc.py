#!/usr/bin/env python

from otpc.commons import *
from otpc import flip, track, stat

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--version', action='version', version='%(prog)s v0.0.1')

    # common args
    parser.add_argument('-o', '--out-dir', type=str, required=True, help="output dir")
    parser.add_argument('-p', '--threads', type=int, default=1, required=False, help="number of threads (default: 1)")
    parser.add_argument('--bam', action='store_true', default=False, required=False, help="")
    parser.add_argument('-b', '--binary', type=str, default=None, required=False, help="")
    parser.add_argument('--bowtie2', action='store_true', default=False, required=False, help="")
    parser.add_argument('--gtf', action='store_true', default=False, required=False, help="")
    parser.add_argument('-l', '--min-exact-match', required=False, help="", \
                default=20, type=int)
    parser.add_argument('--schema', type=lambda s: s.split(','), required=False, \
                default=['transcript', 'ID', 'Parent', 'gene_name', 'transcript_type'])
    parser.add_argument('--keep-dot', required=False, default=False, help="", \
                action='store_true')
    parser.add_argument('--force', required=False, help="", \
                default=False, action='store_true')
    parser.add_argument('--skip-index', required=False, help="", \
                default=False, action='store_true')
    
    subparsers = parser.add_subparsers(dest='module', \
                            help="[flip, track, stat]")

    # flip module
    parser_flip = subparsers.add_parser('flip', help="")
    parser_flip.add_argument('-i', '--in-file', type=str, required=True, \
                            help="input probe sequences (fasta)")
    parser_flip.add_argument('-a', '--src-annotation', type=str, required=True, \
                            help="source transcriptome annotation (gff or gtf)")
    parser_flip.add_argument('-f', '--src-fasta', type=str, required=True, \
                            help="source transcript sequences (fasta)")
    
    # track module
    parser_track = subparsers.add_parser('track', help="")
    parser_track.add_argument('-q', '--query', type=str, required=True, \
                            help="query probe sequences (fasta)")
    parser_track.add_argument('-t', '--target', type=str, required=True, \
                            help="target transcript sequences (fasta)")
    parser_track.add_argument('-a', '--annotation', type=str, required=True, \
                    help="target transcriptome annotation (gff or gtf)")
    parser_track.add_argument('-m', '--mode', type=str, choices=['pad', 'nm'], \
                            required=False, default='pad', help="")
    parser_track.add_argument('-pl', '--pad-length', type=int, required=False, \
                    help="", default=0)
    parser_track.add_argument('-x', '--max-mismatch', type=int, required=False, \
                    help="", default=1)
    
    # stat module
    parser_stat = subparsers.add_parser('stat', help="")
    parser_stat.add_argument('-i', '--in-file', type=str, required=True, \
                    help="track module results (i.e., probe2targets.csv)")
    parser_stat.add_argument('-q', '--query', type=str, required=True, \
                    help="query probe sequences (fasta)")
    parser_stat.add_argument('--exclude-pseudo', required=False, default=False, help="", \
                    action='store_true')
    parser_stat.add_argument('--pc-only', required=False, default=False, help="", \
                    action='store_true')
    parser_stat.add_argument('-s', '--syn-file', type=str, required=False,
                    help="", default=None)
    
    args = parser.parse_args()
    if args.module not in ['flip', 'track', 'stat']:
        parser.error(f"Invalid module {args.module}. Valid options are: flip, track, predict")
    return args

def check_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def check_flip_args(args) -> bool:
    check_dir(args.out_dir)
    return all(os.path.exists(pth) for pth in \
            [args.in_file, args.src_annotation, args.src_fasta])

def check_track_args(args) -> bool:
    check_dir(args.out_dir)
    return all(os.path.exists(pth) for pth in \
            [args.query, args.target, args.annotation])

def check_stat_args(args) -> bool:
    check_dir(args.out_dir)
    return all(os.path.exists(pth) for pth in \
            [args.in_file])

def main() -> None:
    args = parse()
    if args.module == 'flip':
        if not check_flip_args(args):
            print(message(f"cannot locate files", Mtype.ERR))
            sys.exit(-1)
        print(message(f"### FLIP ###", Mtype.PROG))
        param_fn = os.path.join(args.out_dir, "flip_params.json")
        store_params(args, param_fn)
        flip.main(args)
    elif args.module == 'track':
        if not check_track_args(args):
            print(message(f"cannot locate files", Mtype.ERR))
            sys.exit(-1)
        print(message(f"### TRACK ###", Mtype.PROG))
        param_fn = os.path.join(args.out_dir, "track_params.json")
        store_params(args, param_fn)
        track.main(args)
    elif args.module == 'stat':
        if not check_stat_args(args):
            print(message(f"cannot locate files", Mtype.ERR))
            sys.exit(-1)
        if args.pc_only and args.exclude_pseudo:
            print(message(f"cannot use both --pc-only and --exclude-pseudo flags", Mtype.ERR))
            sys.exit(-1)
        print(message(f"### STAT ###", Mtype.PROG))
        param_fn = os.path.join(args.out_dir, "stat_params.json")
        store_params(args, param_fn)
        stat.main(args)

if __name__ == "__main__":
    main()