from otpc.commons import *

def check_tinfo_completeness(qinfos, tinfos) -> set:
    tset = set([x[0] for x in tinfos.values()])
    missing = set()
    for qname in qinfos:
        gid, gname = qinfos[qname]
        if gid not in tset:
            missing.add((gid, gname))
    print(message(f"{len(missing)} target genes missing from the input transcriptome", Mtype.PROG))
    return missing

def load_qinfos(fn):
    qfa = pyfastx.Fasta(fn)
    qinfos = dict()
    for q in qfa:
        temp = q.name.split("|")
        gid = temp[0]
        gname = temp[1]
        qinfos[q.name] = (gid, gname)
    return qinfos, qfa

def write_qinfos(out_dir, qinfos) -> None:
    fn = os.path.join(out_dir, 'probe_infos.csv')
    with open(fn, 'w') as fh:
        fh.write('probe_id,gene_id,gene_name\n')
        for x in qinfos:
            fh.write(f'{x},{qinfos[x][0]},{qinfos[x][1]}\n')

def load_bam(fn, qinfos, tinfos, is_bam) -> dict:
    ainfos = dict()
    fh = pysam.AlignmentFile(fn, 'rb') if is_bam else pysam.AlignmentFile(fn, 'r')
    for brec in fh:
        qname = brec.query_name
        if qname not in ainfos:
            ainfos[qname] = []
        if brec.is_unmapped: continue
        tname = brec.reference_name
        assert qname in qinfos # sanity check
        target_gid, _ = qinfos[qname]
        if target_gid == tinfos[tname][0]:
            ainfos[qname].append(brec.is_forward)
    fh.close()
    return ainfos

def rc(ainfos, qfa, out_dir):
    unaligned = []
    to_rc = []
    for qname in ainfos:
        if len(ainfos[qname]) == 0:
            unaligned.append(qname)
        elif all(ainfos[qname]):
            continue
        else:
            assert all(not x for x in ainfos[qname])
            to_rc.append(qname)
    print(message(f"{len(to_rc)} transcripts to reverse complement", Mtype.PROG))
    fn = os.path.join(out_dir, 'fwd_oriented.fa')
    with open(fn, 'w') as fh:
        for q in qfa:
            if q.name in to_rc:
                seq = Seq(q.seq)
                out_s = seq.reverse_complement()
            else:
                out_s = q.seq
            fh.write(f'>{q.name}\n{out_s}\n')
    return unaligned, to_rc

def main(args) -> None:
    bfn = align(args.query, args.target, 'temp', args)
    att_sep = ' ' if args.gtf else '='

    fn = os.path.join(args.out_dir, 't2g.csv')
    print(message(f"loading target transcriptome infos", Mtype.PROG))
    if os.path.exists(fn):
        tinfos = load_tinfos(fn)
    else:
        tinfos = build_tinfos(args.annotation, att_sep, args.schema, args.keep_dot)
        write_tinfos(args.out_dir, tinfos)
    
    print(message(f"loading query probes infos", Mtype.PROG))
    qinfos, qfa = load_qinfos(args.query)

    missing_t = check_tinfo_completeness(qinfos, tinfos)
    if len(missing_t) > 0: # if >= 1 target gene is missing in the txome
        write_lst(missing_t, os.path.join(args.out_dir, 'missing_in_target.lst'))
    
    print(message(f"parsing alignment results", Mtype.PROG))
    ainfos = load_bam(bfn, qinfos, tinfos, args.bam)
    unaligned, rced = rc(ainfos, qfa, args.out_dir)
    if len(unaligned) > 0:
        write_lst(unaligned, os.path.join(args.out_dir, 'temp.unmapped.lst'))
    if len(rced) > 0:
        write_lst(rced, os.path.join(args.out_dir, 'rev_cmped_probes.lst'))
    print(message(f"wrote forward oriented probes a file", Mtype.PROG))