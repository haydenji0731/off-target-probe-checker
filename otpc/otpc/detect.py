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

def convert_md2bit_nucmer(s, tstart):
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
        bit_s += '1' * (int(running) - tstart)
    return bit_s

def convert_md2bit_del(s):
    running = ""
    bit_s = ""
    ignore = False
    mismatch_info = []
    for c in s:
        if c.isdigit():
            ignore = False
            running += c
        else:
            if ignore: continue
            if c == '^':
                if len(running) > 0:
                    bit_s += '1' * int(running)
                running = ""
                ignore = True
                continue
            else:
                if len(running) > 0:
                    bit_s += '1' * int(running)
                mismatch_info.append(len(bit_s))
                bit_s += '0'
                running = ""
    if len(running) > 0:
        bit_s += '1' * int(running)
    return bit_s, mismatch_info

def convert_md2bit_nucmer_del(s, tstart):
    running = ""
    bit_s = ""
    ignore = False
    mismatch_info = []
    for c in s:
        if c.isdigit():
            ignore = False
            running += c
        else:
            if ignore: continue
            if c == '^':
                if len(running) > 0:
                    bit_s += '1' * int(running)
                running = ""
                ignore = True
                continue
            else:
                if len(running) > 0:
                    bit_s += '1' * int(running)
                mismatch_info.append(len(bit_s))
                bit_s += '0'
                running = ""
    if len(running) > 0:
        bit_s += '1' * (int(running) - tstart)
    return bit_s, mismatch_info

def convert_cigar2bit(tup):
    bit_s = ""
    left_clip = 0
    right_clip = 0
    ins_info = []
    for x in tup:
        op, l = x
        if op == 0: # match
            bit_s += '1' * l
        elif op == 1 or op == 4: # soft clip or ins
            if op == 4:
                if len(bit_s) == 0:
                    left_clip = l
                else:
                    right_clip = l
            if op == 1:
                ins_info.append((bit_s.count('1'), l))
            bit_s += '0' * l
    return bit_s, (left_clip, right_clip), ins_info

def convert_cigar2bit_del(tup, n, mismatch_info):
    bit_s_lst = []
    bit_s = ""
    running = 0
    mut_ctr = 0
    for x in tup:
        op, l = x
        if op == 0: # match
            temp = '1' * l
            # TODO: test if this behaves as expected
            # switch '1' to '0' if a mismatch is reported in this stretch of matches
            if len(mismatch_info) > 0 and mut_ctr != len(mismatch_info):
                temp_lst = list(temp)
                for m in mismatch_info:
                    if m >= running and m < running + l:
                        temp_lst[m - running] = '1'
                        mut_ctr += 1
                temp = ''.join(temp_lst)
            running += l
            bit_s += temp
        elif op == 1 or op == 4: # soft clip or ins
            bit_s += '0' * l
        elif op == 2:
            temp = '0' * len(bit_s)
            bit_s += '0' * (n - len(bit_s))
            bit_s_lst.append(bit_s)
            bit_s = temp
    bit_s_lst.append(bit_s)
    return bit_s_lst

def bitwise_and(s1, s2):
    out = ""
    assert len(s1) == len(s2)
    for i in range(len(s1)):
        if s1[i] == '1' and s2[i] == '1':
            out += '1'
        else:
            out += '0'
    return out

def char2sym(char):
    if char == '0':
        return 'X'
    return '='

def compress_bvec(bvec):
    out = []
    curr_char = bvec[0]
    ctr = 1
    for char in bvec[1:]:
        if char == curr_char:
            ctr += 1
        else:
            out.append(f"{char2sym(curr_char)}{ctr}")
            curr_char = char
            ctr = 1
    out.append(f"{char2sym(curr_char)}{ctr}")
    return ''.join(out)

def detect_off_target(fn, qfa, pad, tinfos, is_nucmer):
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
                cigar_tups = brec.cigartuples
                num_mismatch = int(brec.get_tag('NM'))
                md_tag = brec.get_tag('MD')
                tstart = brec.reference_start
                if cigar == f'{qlen}M':
                    if num_mismatch == 0:
                        ainfos[qname].add((tinfos[tname], f'={qlen}')) # (gid, gname)
                    else:
                        if is_nucmer:
                            md_bvec = convert_md2bit_nucmer(md_tag, tstart)
                        else:
                            md_bvec = convert_md2bit(md_tag)
                        if crit_dvec & int(md_bvec, 2) == crit_dvec:
                            ainfos[qname].add((tinfos[tname], compress_bvec(md_bvec)))
                else:
                    if 'D' in cigar: # handle deletions separately
                        # NOTE: no need to check num_mismatch == 0 as dels count as mismatches (i.e., nm > 0 guaranteed)
                        if is_nucmer:
                            md_bvec, mismatch_info = convert_md2bit_nucmer_del(md_tag, tstart)
                        else:
                            md_bvec, mismatch_info = convert_md2bit_del(md_tag)
                        cigar_bvecs = convert_cigar2bit_del(cigar_tups, qlen, mismatch_info)
                        hit = False
                        final_bvec = None
                        for bvec in cigar_bvecs:
                            if crit_dvec & int(bvec, 2) == crit_dvec:
                                final_bvec = bvec
                                hit = True
                                break
                        if hit:
                            assert final_bvec is not None
                            ainfos[qname].add((tinfos[tname], compress_bvec(final_bvec)))
                    else:
                        cigar_bvec, clip_info, ins_info = convert_cigar2bit(cigar_tups)
                        if num_mismatch == 0: # accounts for cases with just soft clips
                            final_bvec = cigar_bvec # NOTE: ins counts as a mismatch
                        else:
                            if is_nucmer:
                                md_bvec = convert_md2bit_nucmer(md_tag, tstart)
                            else:
                                md_bvec = convert_md2bit(md_tag)
                            if len(ins_info) > 0:
                                # convert md_bvec to account for soft-clipped and inserted bases
                                temp = ""
                                prev_p = None
                                for i in range(len(ins_info)):
                                    p, l = ins_info[i]
                                    if i == 0:
                                        temp += md_bvec[:p] + '0' * l
                                    else:
                                        temp += md_bvec[prev_p:p] + '0' * l
                                    prev_p = p
                                temp += md_bvec[ins_info[-1][0]:]
                                temp = '0' * clip_info[0] + temp + '0' * clip_info[1]
                                assert len(temp) == len(cigar_bvec) # sanity check
                                final_bvec = bitwise_and(cigar_bvec, temp)
                            else:
                                temp = '0' * clip_info[0] + md_bvec + '0' * clip_info[1]
                                assert len(temp) == len(cigar_bvec) # sanity check
                                final_bvec = bitwise_and(cigar_bvec, temp)
                        if crit_dvec & int(final_bvec, 2) == crit_dvec:
                            ainfos[qname].add((tinfos[tname], compress_bvec(final_bvec)))
    return unaligned, ainfos

def write_results(ainfos, d):
    fn = os.path.join(d, 'probe2targets.tsv')
    with open(fn, 'w') as fh:
        fh.write('probe_id\tn_targets\ttarget_ids\ttarget_names\tcigars\ttranscript_types\n') # Caleb: add transcript type
        for qname in ainfos:
            temp = [x[0] for x in ainfos[qname]]
            gids = [x[0] for x in temp]
            gnames = [x[1] for x in temp]
            cigars = [x[1] for x in ainfos[qname]]
            ttypes = [x[2] for x in temp] # Caleb: add transcript type
            gids_s = ','.join(gids)
            gnames_s = ','.join(gnames) # TODO: test if None gene_name values throw an error here
            cigar_s = ','.join(cigars)
            fh.write(f'{qname}\t{len(gids)}\t[{gids_s}]\t[{gnames_s}]\t[{cigar_s}]\t{ttypes}\n') # Caleb: add transcript type

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
    unaligned, ainfos = detect_off_target(bfn, qfa, args.padding, tinfos, args.nucmer)
    write_lst(unaligned, os.path.join(args.out_dir, 'main.unmapped.lst'))
    write_results(ainfos, args.out_dir)
    print(message(f"finished", Mtype.PROG))