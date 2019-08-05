"""
Convert a CCS read name to rid (DazDB ID). Also rescue unphased reads at the end.

E.g.

lookup:
    000000000	Sim/1/ccs
    000000023	Sim/24/ccs
    000000132	Sim/133/ccs
    000000173	Sim/174/ccs
    000001000	Sim/1380/ccs
    000001188	Sim/1568/ccs
    000000744	Sim/1124/ccs
    000000239	Sim/240/ccs
    000000914	Sim/1294/ccs
    000001037	Sim/1417/ccs
    ...
rid_to_ctg:
    Sim/133/ccs 000000F
    Sim/174/ccs 000000F
    Sim/1380/ccs 000000F
    Sim/1568/ccs 000000F
    Sim/24/ccs 000000F
    ...
rid_to_phase:
    Sim/133/ccs 000000F 3000002 1
    Sim/174/ccs 000000F 1000001 0
    Sim/1380/ccs 000000F -1 0
    Sim/1568/ccs 000000F 1000001 1
    ...
output:
    000000132 000000F 3000002 1
    000000173 000000F 1000001 0
    000001000 000000F -1 0
    000001188 000000F 1000001 1
    000000023 000000F -1 0
    ...
"""
import argparse
import sys

def run(lookup, rid_to_phase, rid_to_ctg, output, ctg):
    """Convert a CCS read name to rid (DazDB ID). Also add in unphased reads at the end.
    """
    nameLookup = {}  # CCS read name -> DB ID
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
    'mf' has CCS read-names instead of rids.
    Write into 'out', with actual rid instead of CCS read-name
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
    Write into 'out', with actual rid.
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
        help='The DazzDB ID -> CCS name (%%08d <-> CCS)'
    )
    parser.add_argument(
        '--rid-to-phase',
        required=True,
        help='Rid to phase file, but "rid" is actually CCS Read Name',
    )
    parser.add_argument(
        '--rid-to-ctg',
        required=True,
        help='read name (ccs) to primary ctg'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output file name. Same as rid-to-phase, but actually using "rid"'
    )
    parser.add_argument(
        '--ctg',
        required=True,
        help='Ctg name (for writing unphased reads)'
    )

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
