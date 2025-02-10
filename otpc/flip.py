from commons import *

# TODO: complete use_name logic 

def align(args) -> str:
    fn = os.path.join(args.out_dir, 'temp.bam')
    if args.binary:
        aligner = args.binary
    else:
        aligner = "nucmer" if args.nucmer else "bowtie2"
    if args.nucmer:
        raise NotImplementedError
    else: # bt2
        cmd = f'{aligner}-build {args.target} --threads {args.threads}'
        print(cmd)
        call(cmd, shell=True)
        cmd = f'{aligner} -f -a -N 1 --local -x {args.target} \
            -U {args.query} --very-sensitive-local --threads {args.threads} \
            | samtools view -b -o {fn} -@ {args.threads}'
        print(cmd)
        call(cmd, shell=True)
    return fn

def att2dict(s, sep):
    temp = s.split(';')
    d = dict()
    for x in temp:
        kv = x.split(sep)
        if len(kv) != 2: continue
        k = kv[0].strip()
        v = kv[1].strip()
        d[k] = v
    return d

def load_tinfos(fn, att_sep, schema) -> dict:
    df = pd.read_csv(fn, sep='\t', header=None)
    df.columns = ['ctg', 'src', 'feat', 'start', 'end', 'score', 'strand', 'frame', 'att']
    tinfos = dict()
    ctr = 0
    for i, row in df.iterrows():
        if row['feat'] == schema[0]:
            ctr += 1
            att_d = att2dict(row['att'], att_sep)
            if att_d[1] not in att_d or att_d[2] not in att_d: 
                print(message(f"Invalid schema", Mtype.ERR))
                return None # terminate
            tid = att_d[1]
            gid = att_d[2].split('.')[0] # TODO: make this optional
            tinfos[tid] = gid
    print(message(f"loaded {ctr} transcripts", Mtype.PROG))
    return tinfos

def write_tinfos(out_dir, tinfos) -> None:
    fn = os.path.join(out_dir, 't2g.csv')
    with open(fn, 'w') as fh:
        fh.write('gene_id,transcript_id\n')
        for x in tinfos:
            fh.write(f'{x},{tinfos[x]}\n')

def check_tinfo_completeness(qinfos, tinfos, use_name) -> bool:
    tset = set(tinfos.values())
    missing = set()
    for qname in qinfos:
        gid, gname = qinfos[qname]
        if use_name:
            if gname not in tset: missing.add(gname)
        else:
            if gid not in tset: missing.add(gid)
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

def load_bam(fn, qinfos, tinfos, use_name) -> dict:
    ainfos = dict()
    with pysam.AlignmentFile(fn, 'rb') as fh:
        for brec in fh:
            qname = brec.query_name
            if qname not in ainfos:
                ainfos[qname] = []
            if brec.is_unmapped: continue
            tname = brec.reference_name
            assert qname in qinfos # sanity check
            target_gid, target_gname = qinfos[qname]
            if target_gid == tinfos[tname]:
                ainfos[qname].append(brec.is_forward)
    return ainfos

def flip(ainfos, qfa, out_dir):
    unaligned = []
    to_flip = []
    for qname in ainfos:
        if len(ainfos[qname]) == 0:
            unaligned.append(qname)
        elif all(ainfos[qname]):
            continue
        else:
            assert all(not x for x in ainfos[qname])
            to_flip.append(qname)
    print(message(f"{len(to_flip)} transcripts to reverse complement", Mtype.PROG))
    fn = os.path.join(out_dir, 'fwd_oriented.fa')
    with open(fn, 'w') as fh:
        for q in qfa:
            if q.name in to_flip:
                seq = Seq(q.seq)
                out_s = seq.reverse_complement()
            else:
                out_s = q.seq
            fh.write(f'>{q.name}\n{out_s}\n')
    return unaligned, to_flip

def main(args) -> None:
    bfn = align(args)
    att_sep = ' ' if args.gtf else '='
    tinfos = load_tinfos(args.annotation, att_sep, args.schema)
    qinfos, qfa = load_qinfos(args.query)
    missing_t = check_tinfo_completeness(qinfos, tinfos, args.use_name) # TODO: write out missing
    write_tinfos(args.out_dir, tinfos)
    ainfos = load_bam(bfn)
    unaligned, flipped = flip(ainfos, qfa, args.out_dir)
