#!/usr/bin/env python
from multiprocessing import Pool
import subprocess as sp
import shlex
import argparse
import os
import re
import sys

arid2phase = {}

# TODO: Kill all threads gracefully if one encounters an error.
"""This whole program can hang with something like this:

Exception in thread Thread-3:
Traceback (most recent call last):
  File "/mnt/software/p/python/2.7.13-UCS4/centos-6/lib/python2.7/threading.py", line 801, in __bootstrap_inner
    self.run()
  File "/mnt/software/p/python/2.7.13-UCS4/centos-6/lib/python2.7/threading.py", line 754, in run
    self.__target(*self.__args, **self.__kwargs)
  File "/mnt/software/p/python/2.7.13-UCS4/centos-6/lib/python2.7/multiprocessing/pool.py", line 389, in _handle_results
    task = get()
TypeError: ('__init__() takes at least 3 arguments (1 given)', <class 'subprocess.CalledProcessError'>, ())

But the most typical cause of error is a missing input, so just check that early.
"""

def filter_stage1(input_):
    db_fn, fn, max_diff, max_ovlp, min_ovlp, min_len = input_
    try:
        ignore_rtn = []
        current_q_id = None
        ave_idt = 0.0
        all_over_len = 0.0
        overlap_data = {"5p":0, "3p":0}
        overlap_phase = {"5p":set(), "3p":set()}
        q_id = None

        for l in sp.check_output(shlex.split("LA4Falcon -mo %s %s" % (db_fn, fn) ) ).splitlines():
            l = l.strip().split()
            q_id, t_id = l[:2]

            if q_id not in arid2phase:
                continue
            if t_id not in arid2phase:
                continue
            if arid2phase[t_id][0] != arid2phase[q_id][0]:
                continue

            if arid2phase[t_id][1] == arid2phase[q_id][1] and arid2phase[t_id][2] != arid2phase[q_id][2]:
                continue

            if q_id != None and q_id != current_q_id:

                left_count = overlap_data["5p"]
                right_count = overlap_data["3p"]

                if abs(left_count - right_count) > max_diff:
                    ignore_rtn.append( current_q_id )
                elif left_count > max_ovlp or right_count > max_ovlp:
                    ignore_rtn.append( current_q_id )
                elif left_count < min_ovlp or right_count < min_ovlp:
                    ignore_rtn.append( current_q_id )

                #remove unphased reads if they are sandwiched by the reads from the same phase
                if current_q_id not in arid2phase:
                    left_set = overlap_phase["5p"]
                    right_set = overlap_phase["3p"]
                    if len( left_set & right_set ) > 0:
                        ignore_rtn.append( current_q_id )

                overlap_data = {"5p":0, "3p":0}
                overlap_phase = {"5p":set(), "3p":set()}
                current_q_id = q_id
                ave_idt = 0.0
                all_over_len = 0.0

            overlap_len = -int(l[2])
            idt = float(l[3])
            q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
            t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])

            if idt < 90:
                continue

            if q_l < min_len or t_l < min_len:
                continue

            if l[-1] in ("contains", "overlap"):
                ave_idt += idt * overlap_len
                all_over_len += overlap_len
            if q_s == 0:
                overlap_data["5p"] += 1
                if t_id in arid2phase:
                    overlap_phase["5p"].add( arid2phase[t_id] )
            if q_e == q_l:
                overlap_data["3p"] += 1
                if t_id in arid2phase:
                    overlap_phase["3p"].add( arid2phase[t_id] )

        if q_id !=  None:
            left_count = overlap_data["5p"]
            right_count = overlap_data["3p"]
            if abs(left_count - right_count) > max_diff:
                ignore_rtn.append( current_q_id )
            elif left_count > max_ovlp or right_count > max_ovlp:
                ignore_rtn.append( current_q_id )
            elif left_count < min_ovlp or right_count < min_ovlp:
                ignore_rtn.append( current_q_id )

            #remove unphased reads if they are sandwiched by the reads from the same phase
            if current_q_id not in arid2phase:
                left_set = overlap_phase["5p"]
                right_set = overlap_phase["3p"]
                if len( left_set & right_set ) > 0:
                    ignore_rtn.append( current_q_id )


        return fn, ignore_rtn

    except (KeyboardInterrupt, SystemExit):
        return

def filter_stage2(input_):
    db_fn, fn, max_diff, max_ovlp, min_ovlp, min_len, ignore_set = input_
    try:
        contained_id = set()
        for l in sp.check_output(shlex.split("LA4Falcon -mo %s %s" % (db_fn, fn))).splitlines():
            l = l.strip().split()
            q_id, t_id = l[:2]


            if q_id not in arid2phase:
                continue
            if t_id not in arid2phase:
                continue

            if arid2phase[t_id][0] != arid2phase[q_id][0]:
                continue
            if arid2phase[t_id][1] == arid2phase[q_id][1] and arid2phase[t_id][2] != arid2phase[q_id][2]:
                continue


            q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
            t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])

            idt = float(l[3])
            if idt < 90:
                continue

            if q_l < min_len or t_l < min_len:
                continue

            if q_id in ignore_set:
                continue
            if t_id in ignore_set:
                continue
            if l[-1] == "contained":
                contained_id.add(q_id)
            if l[-1] == "contains":
                contained_id.add(t_id)
        return fn, contained_id

    except (KeyboardInterrupt, SystemExit):
        return

def filter_stage3(input_):

    db_fn, fn, max_diff, max_ovlp, min_ovlp, min_len, ignore_set, contained_set, bestn = input_

    overlap_data = {"5p":[], "3p":[]}
    try:
        ovlp_output = []
        current_q_id = None
        for l in sp.check_output(shlex.split("LA4Falcon -mo %s %s" % (db_fn, fn) )).splitlines():
            l = l.strip().split()
            q_id, t_id = l[:2]


            if q_id not in arid2phase:
                continue
            if t_id not in arid2phase:
                continue

            if arid2phase[t_id][0] != arid2phase[q_id][0]:
                continue
            if arid2phase[t_id][1] == arid2phase[q_id][1] and arid2phase[t_id][2] != arid2phase[q_id][2]:
                continue

            if current_q_id == None:
                current_q_id = q_id
                overlap_data = {"5p":[], "3p":[]}

            elif q_id != current_q_id:
                left = overlap_data["5p"]
                right = overlap_data["3p"]
                left.sort()
                right.sort()

                for i in xrange(len(left)):
                    inphase, score, m_range, ovlp = left[i]
                    ovlp_output.append(ovlp)
                    #print " ".join(ovlp), read_end_data[current_q_id]
                    if i >= bestn and m_range > 1000:
                        break

                for i in xrange(len(right)):
                    inphase, score, m_range, ovlp = right[i]
                    ovlp_output.append(ovlp)
                    #print " ".join(ovlp), read_end_data[current_q_id]
                    if i >= bestn and m_range > 1000:
                        break

                overlap_data = {"5p":[], "3p":[]}
                current_q_id = q_id

            if q_id in contained_set:
                continue
            if t_id in contained_set:
                continue
            if q_id in ignore_set:
                continue
            if t_id in ignore_set:
                continue

            overlap_len = -int(l[2])
            idt = float(l[3])
            q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
            t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])

            if idt < 90:
                continue
            if q_l < min_len or t_l < min_len:
                continue
            if q_s == 0:
                l.extend( [".".join(arid2phase.get(current_q_id, "NA")), ".".join(arid2phase.get(t_id, "NA"))])
                inphase = 1 if arid2phase.get(current_q_id, "NA") == arid2phase.get(t_id, "NA") else 0
                #nphase = 1 if arid2phase[current_q_id] == arid2phase[t_id] else 0
                overlap_data["5p"].append( (-inphase, -overlap_len,  t_l - (t_e - t_s),  l ) )
            elif q_e == q_l:
                l.extend( [".".join(arid2phase.get(current_q_id, "NA")), ".".join(arid2phase.get(t_id, "NA"))])
                inphase = 1 if arid2phase.get(current_q_id, "NA") == arid2phase.get(t_id, "NA") else 0
                #inphase = 1 if arid2phase[current_q_id] == arid2phase[t_id] else 0
                overlap_data["3p"].append( (-inphase, -overlap_len, t_l - (t_e - t_s), l ) )

        left = overlap_data["5p"]
        right = overlap_data["3p"]
        left.sort()
        right.sort()


        for i in xrange(len(left)):
            inphase, score, m_range, ovlp = left[i]
            ovlp_output.append(ovlp)
            #print " ".join(ovlp), read_end_data[current_q_id]
            if i >= bestn and m_range > 1000:
                break

        for i in xrange(len(right)):
            inphase, score, m_range, ovlp = right[i]
            ovlp_output.append(ovlp)
            #print " ".join(ovlp), read_end_data[current_q_id]
            if i >= bestn and m_range > 1000:
                break

        return fn, ovlp_output
    except (KeyboardInterrupt, SystemExit):
        return


def parse_args(argv):

    parser = argparse.ArgumentParser(description='a simple multi-processes LAS ovelap data filter')
    parser.add_argument('--n_core', type=int, default=4,
                        help='number of processes used for generating consensus')
    parser.add_argument('--fofn', type=str, help='file contains the path of all LAS file to be processed in parallel')
    parser.add_argument('--db', type=str, help='read db file path')
    parser.add_argument('--max_diff', type=int, help="max difference of 5' and 3' coverage")
    parser.add_argument('--max_cov', type=int, help="max coverage of 5' or 3' coverage")
    parser.add_argument('--min_cov', type=int, help="min coverage of 5' or 3' coverage")
    parser.add_argument('--min_len', type=int, default=2500, help="min length of the reads")
    parser.add_argument('--bestn', type=int, default=10, help="output at least best n overlaps on 5' or 3' ends if possible")
    parser.add_argument('--rid_phase_map', type=str, help="the file that encode the relationship of the read id to phase blocks", required=True)
    args = parser.parse_args(argv[1:])

    return args

def assert_exists(fn):
    if not os.path.exists(fn):
        msg = 'File does not exist: {!r}'.format(fn)
        raise Exception(msg)

def main(argv=sys.argv):
    args = parse_args(argv)

    max_diff = args.max_diff
    max_cov = args.max_cov
    min_cov = args.min_cov
    min_len = args.min_len
    bestn = args.bestn
    db_fn = args.db

    assert_exists(db_fn)

    with open(args.rid_phase_map) as f:
        for row in f:
            row = row.strip().split()
            arid2phase[row[0]] = (row[1], row[2], row[3]) #ctg_id, phase_blk_id, phase_id
    assert arid2phase, 'Empty rid_phase_map: {!r}'.format(args.rid_phase_map)
    exe_pool = Pool(args.n_core)

    file_list = open(args.fofn).read().split("\n")
    inputs = []
    for fn in file_list:
        if not fn: continue
        assert_exists(fn)
        inputs.append( (db_fn, fn, max_diff, max_cov, min_cov, min_len) )

    ignore_all = []
    for res in exe_pool.imap(filter_stage1, inputs):
        ignore_all.extend( res[1] )

    inputs = []
    ignore_all = set(ignore_all)
    for fn in file_list:
        if len(fn) != 0:
            inputs.append( (db_fn, fn, max_diff, max_cov, min_cov, min_len, ignore_all) )
    contained = set()
    for res in exe_pool.imap(filter_stage2, inputs):
        contained.update(res[1])
        #print res[0], len(res[1]), len(contained)

    #print "all", len(contained)
    inputs = []
    ignore_all = set(ignore_all)
    for fn in file_list:
        if len(fn) != 0:
            inputs.append( (db_fn, fn, max_diff, max_cov, min_cov, min_len, ignore_all, contained, bestn) )
    for res in exe_pool.imap(filter_stage3, inputs):
        for l in res[1]:
            print " ".join(l)
