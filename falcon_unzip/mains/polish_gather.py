import argparse
import logging
import os
import sys
from .. import io

def run(all_ctgs_json,
        p_ctg_fn,
        h_ctg_fn):
    """
    Deserialize list of polished p&h ctgs then write them to the appropriate outputs
    """
    fns = io.deserialize(all_ctgs_json)
    io.touch(h_ctg_fn)
    io.touch(p_ctg_fn)
    for fn in fns:
        fn = '../../' + fn
        if is_haplotig(fn):
            call = 'cat {} >> {}'.format(fn, h_ctg_fn)
        else:
            call = 'cat {} >> {}'.format(fn, p_ctg_fn)
        io.syscall(call)
        
def is_haplotig(fn):
    """
    >>> is_haplotig('/a/foo_bar.fasta')
    True
    >>> is_haplotig('/a/foo.fasta')
    False
    """
    return '_' in os.path.basename(fn)

def parse_args(argv):
    description = 'Concatenate results of polishing task.'
    epilog = '.'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
    )
    parser.add_argument(
        '--all-ctgs-json', 
        required=True, 
        help="JSON formatted file containing the files to combine"
    )
    parser.add_argument(
        '--p-ctg-fn',
        required=True,
        help='Polished primary ctg. output file',
    )
    parser.add_argument(
        '--h-ctg-fn',
        required=True,
        help='Polished haplotig ctg. output file',
    )
    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))
    
if __name__ == '__main__':  # pragma: no cover
    main()
