"""
These are top-level "tasks", which ought to become
(or be included in) real tasks someday. This will
ease translation into WDL, etc.

Also, these are testable.
"""
from ..io import serialize
import logging

LOG = logging.getLogger(__name__)

def fai2ctgs(p_ctg_fai_fn, o_fn):
    CTGS = []

    with open(p_ctg_fai_fn) as f:
          for line in f:
            lc = line.strip().split("\t")
            if(lc[0][0] == "0"):
                CTGS.append(lc[0])

    serialize(o_fn, CTGS)

def getPH(combined_ph_fai_fn, readtoctg_fn):
    SIZES = {}

    with open(combined_ph_fai_fn) as fai:
        for line in fai:
            lc = line.strip().split("\t")
            SIZES[lc[0]] = lc[1]

    PH = {}

    with open(readtoctg_fn) as f:
        for line in f:
            if line.startswith('#'):
                continue
            lc = line.strip().split(" ")
            PH[lc[1]] = 1 + PH.get(lc[1], 0)

    SKIP = set()
    for (ctg, count) in PH.items():
        if not ctg in SIZES:
            LOG.warning("ctg {} is being skipped because it did not unzip".format(ctg))
            SKIP.add(ctg)
            continue
        if count < 5:
            LOG.warning("ctg {} is being skipped due to depth of coverage {} < 10 reads total".format(ctg, count))
            SKIP.add(ctg)
        if int(SIZES[ctg]) < 10000:
            LOG.warning("ctg {} is being skipped due to size {} < 10kbp".format(ctg, SIZES[ctg]))
            SKIP.add(ctg)

    for d in SKIP:
        del PH[d]

    return PH
