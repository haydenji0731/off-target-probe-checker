from otpc.commons import *

def convert_md2bit(s):
    running = ""
    bit_s = ""
    for c in s:
        if c.isdigit():
            running += c
        else:
            if len(running) > 0:
                bit_s += '1' * int(running)
            bit_s += '0'
            running = ""
    if len(running) > 0:
        bit_s += '1' * int(running)
    return bit_s

def convert_cigar2bit(tup):
    bit_s = ""
    for x in tup:
        op, l = x
        if op == 0: # match
            bit_s += '1' * l
        elif op == 1 or op == 4: # soft clip or ins
            bit_s += '0' * l
    return bit_s

def convert_cigar2bit_del(tup, n):
    bit_s_lst = []
    bit_s = ""
    for x in tup:
        op, l = x
        if op == 0: # match
            bit_s += '1' * l
        elif op == 1 or op == 4: # soft clip or ins
            bit_s += '0' * l
        elif op == 2:
            temp = '0' * len(bit_s)
            bit_s += '0' * (n - len(bit_s))
            bit_s_lst.append(bit_s)
            bit_s = temp
    bit_s_lst.append(bit_s) # TODO: check if this correct
    return bit_s_lst

def detect_off_target(fn, qfa, pad, tinfos):
    unaligned = []
    ainfos = dict()
    with pysam.AlignmentFile(fn, 'rb') as fh:
        for brec in fh:
            qname = brec.query_name
            if brec.is_unmapped:
                unaligned.append(qname)
                continue
            elif brec.is_supplementary:
                continue
            else:
                tname = brec.reference_name
                qlen = len(qfa[qname].seq)
                crit_bvec = "0" * pad + "1" * (qlen - 2 * pad) + "0" * pad
                assert len(crit_bvec) == qlen # sanity check
                crit_dvec = int(crit_bvec, 2)
                if qname not in ainfos:
                    ainfos[qname] = set()
                cigar = brec.cigarstring
                if cigar == f'{qlen}M':
                    num_mismatch = int(brec.get_tag('NM'))
                    if num_mismatch == 0:
                        ainfos[qname].add(tinfos[tname]) # (gid, gname)
                    else:
                        md_tag = brec.get_tag('MD')
                        md_bvec = convert_md2bit(md_tag)
                        if crit_dvec & int(md_bvec, 2) == crit_dvec:
                            ainfos[qname].add(tinfos[tname])
                else:
                    if 'D' in cigar: # handle deletions separately
                        md_bvecs = convert_cigar2bit_del(brec.cigartuples, qlen)
                        hit = False
                        for bvec in md_bvecs:
                            if crit_dvec & int(bvec, 2) == crit_dvec:
                                hit = True
                                break
                        if hit:
                           ainfos[qname].add(tinfos[tname])
                    else:
                        md_bvec = convert_cigar2bit(brec.cigartuples)
                        if crit_dvec & int(md_bvec, 2) == crit_dvec:
                            ainfos[qname].add(tinfos[tname])
    return unaligned, ainfos

def write_results(ainfos, d):
    fn = os.path.join(d, 'probe2targets.tsv')
    with open(fn, 'w') as fh:
        fh.write('probe_id\tn_targets\ttarget_ids\ttarget_names\n')
        for qname in ainfos:
            temp = ainfos[qname]
            gids = [x[0] for x in temp]
            gnames = [x[1] for x in temp]
            gids_s = ','.join(gids)
            gnames_s = ','.join(gnames)
            fh.write(f'{qname}\t{len(gids)}\t[{gids_s}]\t[{gnames_s}]\n')

def main(qfn, bfn, args) -> None:
    qfa = pyfastx.Fasta(qfn)

    # check if t2g.csv file already exists
    fn = os.path.join(args.out_dir, 't2g.csv')
    if not os.path.exists(fn) or args.force:
        att_sep = ' ' if args.gtf else '='
        tinfos = build_tinfos(args.annotation, att_sep, args.schema, args.keep_dot)
        write_tinfos(args.out_dir, tinfos)
    else:
        tinfos = load_tinfos(fn)
    
    print(message(f"detecting potential off-target probe activities", Mtype.PROG))
    unaligned, ainfos = detect_off_target(bfn, qfa, args.padding, tinfos)
    write_lst(unaligned, os.path.join(args.out_dir, 'main.unmapped.lst'))
    write_results(ainfos, args.out_dir)
    print(message(f"finished", Mtype.PROG))