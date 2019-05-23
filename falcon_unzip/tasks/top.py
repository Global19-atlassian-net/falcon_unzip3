"""
These are top-level "tasks", which ought to become
(or be included in) real tasks someday. This will
ease translation into WDL, etc.

Also, these are testable.
"""

from ..io import serialize

def p_ctg_fai2ctgs(p_ctg_fai_fn, o_fn):
    CTGS = []

    with open(p_ctg_fai_fn) as f:
          for line in f:
            lc = line.strip().split("\t")
            if(lc[0][0] == "0"):
                CTGS.append(lc[0])

    serialize(o_fn, CTGS)
