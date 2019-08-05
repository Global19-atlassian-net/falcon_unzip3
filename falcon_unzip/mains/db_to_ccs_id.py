"""
This script rescues unphased reads.
"""
import argparse
import sys

def run(lookup, rid_to_phase, rid_to_ctg, output, ctg):
    """Convert a CCS read name to a DazDB ID. Also add in unphased reads at the end.
    """
    nameLookup = {}
    with open(lookup, 'r') as lf:
        for line in lf:
            lc = line.strip().split("\t")
            nameLookup[lc[1]] = lc[0]

    # 'output' has the same format as the rid_to_phase file. (See below.)

    with open(output, 'w') as out:
        with open(rid_to_phase, 'r') as mf:
            seen = write_rid_to_phase(out, mf, nameLookup)
        with open(rid_to_ctg, 'r') as rctg:
            write_rid_to_ctg(out, rctg, nameLookup, seen, ctg)


def write_rid_to_phase(out, mf, nameLookup):
    """Write phased reads.
    Read from 'mf', using only reads from 'nameLookup'.
    Write into 'out'.
    """
    seen = {}

    #rid_to phase file contains four columns (white space split): DAZdb_ID primary_ctg phase_block phase
    #000322297 000000F 9000023 1

    for line in mf:
        lc = line.strip().split(" ")
        if lc[0] in nameLookup:
            out.write("%s %s %s %s\n" % (nameLookup[lc[0]], lc[1], lc[2], lc[3]))
            seen[lc[0]] = 1
    return seen


def write_rid_to_ctg(out, rctg, nameLookup, seen, ctg):
    """Write unphased reads.
    Read from 'rctg', using only reads from 'nameLookup'.
    Write into 'out'
    """
    for line in rctg:
        lc = line.strip().split(" ")
        if lc[0] in nameLookup:
            if lc[0] not in seen:
                out.write("%s %s %s %s\n" % (nameLookup[lc[0]], ctg, -1, 0))


def parse_args(argv):
    description = 'Remap DB ids (%%08d) to the CCS read name'
    parser = argparse.ArgumentParser(
        description=description,
    )
    parser.add_argument(
        '--lookup',
        required=True,
        help='The CCS name ID lookup (%%08d <-> CCS)'
    )
    parser.add_argument(
        '--rid-to-phase',
        required=True,
        help='Rid to phase file',
    )
    parser.add_argument(
        '--rid-to-ctg',
        required=True,
        help='read name (ccs) to primary ctg'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output file name'
    )
    parser.add_argument(
        '--ctg',
        required=True,
        help='Ctg name'
    )

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
