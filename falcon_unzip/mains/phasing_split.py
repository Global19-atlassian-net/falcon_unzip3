import argparse
import logging
import os
import sys
from .. import io
from ..tasks.unzip import TASK_PHASING_RUN_SCRIPT

LOG = logging.getLogger()


def read_ctg_ids(stream):
    """For now, simply one ctg_id per line.
    """
    for line in stream:
        yield line.strip()

def run(base_dir, ctg_list_fn, rawread_ids_fn, pread_ids_fn, pread_to_contigs_fn, split_fn, bash_template_fn):
    LOG.info('Splitting ctg_ids from {!r} into {!r}.'.format(
        ctg_list_fn, split_fn))
    with open(bash_template_fn, 'w') as stream:
        stream.write(TASK_PHASING_RUN_SCRIPT)
    base_dir = os.path.abspath(base_dir)
    #basedir = os.path.dirname(os.path.abspath(scattered_fn))
    #rootdir = os.path.dirname(os.path.dirname(os.path.dirname(basedir))) # for now
    reads_dir = os.path.dirname(os.path.abspath(ctg_list_fn))
    jobs = list()
    for i, ctg_id in enumerate(read_ctg_ids(open(ctg_list_fn))):
        #job_name = 'phasing-{:03d}'.format(i) # wildcard value

        job = dict()
        job['input'] = dict( # These were generated by fetch_reads.py, so we find them by convention.
                ref_fasta = '{reads_dir}/{ctg_id}/ref.fa'.format(**locals()),
                read_fasta = '{reads_dir}/{ctg_id}/reads.fa'.format(**locals()),
                rawread_ids = rawread_ids_fn,
                pread_ids = pread_ids_fn,
                pread_to_contigs = pread_to_contigs_fn,
        )
        job['output'] = dict(
                rid_to_phase_out = 'rid_to_phase'.format(**locals()),
        )
        job['params'] = dict(
                ctg_id=ctg_id,
                base_dir=base_dir,
        )
        job['wildcards'] = dict(
                ctg_id=ctg_id, # This should match the wildcard used in the pattern elsewhere.
        )
        jobs.append(job)
    io.serialize(split_fn, jobs)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Split the tasks to phase the reads.'
    epilog = 'To learn about inputs and outputs in the serialized jobs, grep repo for TASK_RUN_PHASING_SCRIPT.'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--ctg-list-fn', required=True,
        help='Input. For now, each line is a ctg_id.',
    )
    parser.add_argument(
        '--rawread-ids-fn', required=True,
        help='Input propagated to phasing-run-script.',
    )
    parser.add_argument(
        '--pread-ids-fn', required=True,
        help='Input propagated to phasing-run-script.',
    )
    parser.add_argument(
        '--pread-to-contigs-fn', required=True,
        help='Input propagated to phasing-run-script.',
    )
    parser.add_argument(
        '--base-dir', required=True,
        help='Path to run-dir. (Parent of 3-unzip/) We need this because we have under-specified some inputs.')
    parser.add_argument(
        '--split-fn',
        help='Output. JSON list of units of work.')
    parser.add_argument(
        '--bash-template-fn',
        help='Output. Copy of known bash template, for use later.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()