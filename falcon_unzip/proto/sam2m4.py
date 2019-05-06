#! /usr/bin/env python2.7

import re
import sys
from ..io import AlignmentFile
import logging

LOG = logging.getLogger() # root, to inherit from sub-loggers

CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQ = 7
CIGAR_X = 8
CIGAR_B = 9

def is_primary(sam):
    return not (sam.is_secondary or sam.is_supplementary)

def open_sam_bam_for_reading(file_path):
    if (file_path.endswith('bam')):
        fp = AlignmentFile(file_path, 'rb', check_sq=False)
        return fp, True
    fp = AlignmentFile(file_path, 'r', check_sq=False)
    return fp, False

def open_sam_bam_for_writing(file_path, header):
    if (file_path.endswith('bam')):
        fp = AlignmentFile(file_path, 'wb', header=header)
    else:
        fp = AlignmentFile(file_path, 'wh', header=header)
    return fp

def pysam_to_m4(aln, ref_lens = None, skip_supplementary=True, skip_secondary=True):
    ret = None

    aln_query_len = (aln.query_alignment_end - aln.query_alignment_start)

    if aln.is_unmapped == True:
        return ret
    if aln.is_supplementary == True and skip_supplementary == True:
        return ret
    if aln.is_secondary == True and skip_secondary == True:
        return ret
    if len(aln.cigar) == 0:
        return ret
    if aln.query_sequence == None or len(aln.query_sequence) == 0:
        return ret
    assert aln_query_len >= 0, 'Negative value should not be possible: {}-{}=={}'.format(
            aln.query_alignment_end, aln.query_alignment_start, aln_query_len)
    if aln_query_len == 0:
        msg = 'Alignment length should not be zero for alignments marked as mapped! (qname = {}) Skipping.'.format(
            aln.qname)
        LOG.warning(msg)
        return ret

    ref_len = (ref_lens[aln.reference_name]) if ref_lens != None else 0

    m = 0
    matches = 0
    mismatches = 0
    ins = 0
    dels = 0

    for op, count in aln.cigar:

        if op == CIGAR_EQ:
            matches += count
            m += count
        elif op == CIGAR_X:
            mismatches += count
            m += count
        elif op == CIGAR_M:
            # assert(op != CIGAR_M and "Cannot calculate match rate from 'M' operations.'")
            matches += count
            mismatches += count
            m += count
        elif op == CIGAR_I:
            ins += count
        elif op == CIGAR_D:
            dels += count

    score = 0 - (mismatches + ins + dels)
    identity = 100.0 * float(matches) / float(aln_query_len)
    ret = [aln.qname, aln.reference_name, score, identity,
                0, aln.query_alignment_start, aln.query_alignment_end, aln.query_length,
                (1 if aln.is_reverse == True else 0), aln.reference_start, aln.reference_end, ref_len, 254, aln]

    return ret

def sam_to_m4(in_path, chr_id = None):
    ret = []

    fp_in, is_bam = open_sam_bam_for_reading(in_path)
    it_fp_in = fp_in.fetch() if not is_bam else fp_in

    ref_lens = {}
    for val in fp_in.header['SQ']:
        ref_lens[val['SN']] = val['LN']

    num_lines = 0
    for sam in it_fp_in:
        num_lines += 1
        if (num_lines % 10000 == 0):
            LOG.info('Processed {} lines, total {} alignments loaded so far.'.format((num_lines - 1), len(ret)))

        m4 = pysam_to_m4(sam, ref_lens)
        if m4 == None:
            continue
        if chr_id and len(m4) > 2 and m4[1] != chr_id:
            continue
        ret.append(m4)

    return ret
