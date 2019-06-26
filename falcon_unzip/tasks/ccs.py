from falcon_kit.pype import Dist
#from falcon_kit import pype_tasks
from falcon_kit.pype import (wrap_gen_task as gen_task, gen_parallel_tasks, Dist)
from .. import io
from . import top
import glob
import logging
import os

LOG = logging.getLogger(__name__)

TASK_POLISH="""
## gathering phased and unphased reads
samtools faidx {input.FA} {params.ctg} > ref.fasta
perl -lane 'print $F[0] if $F[1] eq "{params.ctg}" && $F[2] != -1' {input.READTOCTG} | sort | uniq > readnames.txt
samtools view -F 1796 {input.UBAM} {params.ctg} | cut -f 1 | sort | uniq >> readnames.txt
samtools fqidx -r readnames.txt {input.FQ} > reads.fastq
time pbmm2 align --sort -j {params.pypeflow_nproc} --preset CCS ref.fasta reads.fastq | samtools view -F 1796  > aln.sam
time racon -t {params.pypeflow_nproc} reads.fastq aln.sam ref.fasta > {output.POL}
# Cleanup 
rm *.sam *.fastq
"""

TASK_GATHER_UNPHASED="""
find {input.ctg*} > bam_list.txt
samtools merge -b bam_list.txt tmp.bam
samtools sort -o {output.UBAM} tmp.bam
samtools index {output.UBAM}
"""

TASK_PLACE_UNPHASED="""
python3 -m falcon_unzip.mains.polish_unphased_readmapping --ctg {params.ctg} --fai {input.fai} --out-read-names read_names.txt --out-ref-names ref_names.txt --read-to-ctg {input.readtoctg}
samtools fqidx -r read_names.txt {input.fq} > reads.fastq
samtools faidx -r ref_names.txt {input.fa}  > ref.fasta
pbmm2 align -j 1 --preset CCS --sort ref.fasta reads.fastq | samtools view -bS -F 1796 > aln.bam
"""

TASK_POLISH_PREAMBLE="""
cat {input.P} {input.H} > {output.COMBINED}
samtools faidx {output.COMBINED}
cat ../../3-unzip/2-htigs/chunk_*/uow-*/*ctg_edges* > {output.EDGES}
python3 -m falcon_unzip.mains.polish_read_to_ctg --rid-to-phase-fn {input.RIDTOPHASE} --edges-fn {output.EDGES} --lookup-fn {input.READNAMELOOKUP} --out {output.READTOCTG}

# Assume FQ is external. (TODO: Use a program so we can be more lenient.)
ln -sf {input.FQ} {output.FQO}
samtools fqidx {output.FQO}
"""

TASK_HASM_COLLECT_SCRIPT = """\
## prepare for quviering the haplotig
## (assume we are in 3-unzip/somewhere/)

# TODO: Stop using job_done.
python3 -m falcon_unzip.mains.graphs_to_h_tigs_2 combine --results-fn={input.results} --done-fn={output.job_done}

find ./0-phasing -name "phased_reads" | sort | xargs cat >| all_phased_reads
#find ./2-htigs -name "h_ctg_ids.*" | sort | xargs cat >| all_h_ctg_ids
#find ./2-htigs -name "p_ctg_edges.*" | sort | xargs cat >| all_p_ctg_edges
#find ./2-htigs -name "h_ctg_edges.*" | sort | xargs cat >| all_h_ctg_edges
#find ./2-htigs -name "p_ctg.*.fa" | sort | xargs cat >| all_p_ctg.fa
#find ./2-htigs -name "h_ctg.*.fa" | sort | xargs cat >| all_h_ctg.fa

if [[ ! -s all_p_ctg.fa ]]; then
    echo "Empty all_p_ctg.fa -- No point in continuing!"
    exit 1
fi

samtools faidx all_p_ctg.fa

# # Generate a GFA for only primary contigs and haplotigs.
# time python3 -m falcon_unzip.mains.unzip_gen_gfa_v1 --unzip-root . --p-ctg-fasta ./all_p_ctg.fa --h-ctg-fasta ./all_h_ctg.fa --preads-fasta {input.preads4falcon} >| ./asm.gfa

# # Generate a GFA of all assembly graph edges. This GFA can contain
# # edges and nodes which are not part of primary contigs and haplotigs
# time python3 -m falcon_unzip.mains.unzip_gen_gfa_v1 --unzip-root . --p-ctg-fasta ./all_p_ctg.fa --h-ctg-fasta ./all_h_ctg.fa --preads-fasta {input.preads4falcon} --add-string-graph >| ./sg.gfa
"""


TASK_GRAPH_TO_H_TIGS_SCRIPT = """\
asm_dir=$(dirname {input.falcon_asm_done})
hasm_dir=$(dirname {input.p_ctg})

python3 -m falcon_unzip.mains.graphs_to_h_tigs_2 --gathered-rid-to-phase={input.gathered_rid_to_phase} --base-dir={params.topdir} --fc-asm-path ${{asm_dir}} --fc-hasm-path ${{hasm_dir}} --ctg-id all --rid-phase-map {input.rid_to_phase_all} --fasta {input.preads4falcon}

# more script -- a little bit hacky here, we should improve

#WD=$PWD
# for f in `cat ../reads/ctg_list `; do mkdir -p $WD/$f; cd $WD/$f; python3 -m falcon_unzip.mains.dedup_h_tigs $f; done

for f in `cat ../reads/ctg_list `
do
    mkdir -p ./$f;
    if [ -s ./$f/h_ctg.$f.fa ]
    then
        grep ">" ./$f/h_ctg.$f.fa | sed "s/^>//" >| ./$f/h_ctg_ids.$f
    else
        rm -rf ./$f/h_ctg_ids.$f
        touch ./$f/h_ctg_ids.$f
    fi
done

touch {output.htigs_done}
"""

TASK_GRAPH_TO_H_TIGS_SPLIT_SCRIPT = """\
asm_dir=$(dirname {input.falcon_asm_done})
hasm_dir=$(dirname {input.p_ctg})

# TODO: Should we look at ../reads/ctg_list ?

python3 -m falcon_unzip.mains.graphs_to_h_tigs_2 split \
        --gathered-rid-to-phase={input.gathered_rid_json} --base-dir={params.topdir} \
        --fc-asm-path ${{asm_dir}} --fc-hasm-path ${{hasm_dir}} \
        --rid-phase-map {input.rid_to_phase_all} --fasta {input.preads4falcon} \
        --split-fn={output.split} --bash-template-fn={output.bash_template}

# The bash-template is just a dummy, for now.
"""


TASK_HASM_SCRIPT = """\
# TODO: Needs preads.db

rm -f ./ctg_paths
python3 -m falcon_unzip.mains.ovlp_filter_with_phase_strict \
        --fofn {input.las_fofn} --max-diff 120 --max-cov 120 --min-cov 1 \
        --n-core {params.pypeflow_nproc} --min-len 2500 --db {input.preads_db} \
        --rid-phase-map {input.rid_to_phase_all} > preads.p_ovl
python3 -m falcon_unzip.mains.phased_ovlp_to_graph preads.p_ovl --min-len 2500 > fc.log

if [[ ! -e ./ctg_paths ]]; then
    exit 1
fi

rm -f preads.p_ovl

# Create haplotigs in a safe manner.

ln -sf {input.preads4falcon} .

rm -f {output.p_ctg}

# Given sg_edges_list, utg_data, ctg_paths, preads4falcon.fasta,
# write p_ctg.fa and a_ctg_all.fa,
# plus p_ctg_tiling_path, a_ctg_tiling_path:
time python3 -m falcon_kit.mains.graph_to_contig

if [[ ! -e {output.p_ctg} ]]; then
    exit 1
fi
"""

TASK_PHASING_GATHER_SCRIPT = """\
# creates a master table of rid to phase
cat {input.ctg*} > {output.rid_to_phase_all}

# creates the needed gathering JSON
find {input.ctg*} | xargs -I [] readlink -f [] | python3 -m falcon_unzip.mains.gen_rid_gathered_json > {output.gathered_rid_json}
"""

TASK_READ_PHASING = """

        #TODO: break up command, and maybe remove some deps.

        python3 -m falcon_unzip.mains.phasing_make_het_call --ctg-id {params.ctg} --bam-fn {input.BAM} --fasta-fn {input.T} --vmap-fn het_calls.ctg.vmap --vpos-fn het_calls.ctg.vpos --q-id-map-fn het_calls.ctg.msgpack
        python3 -m falcon_unzip.mains.phasing_generate_association_table --ctg-id {params.ctg} --vmap=het_calls.ctg.vmap --atable=association_table.ctg.astab
        python3 -m falcon_unzip.mains.phasing_get_phased_blocks --vmap=het_calls.ctg.vmap --atable=association_table.ctg.astab --p-variant=phased_vars.ctg.phased.txt
        python3 -m falcon_unzip.mains.phasing_get_phased_reads --ctg-id={params.ctg} --vmap=het_calls.ctg.vmap --p-variant=phased_vars.ctg.phased.txt --q-id-map=het_calls.ctg.msgpack --phased-reads=phased_reads.ctg.phased.txt

        #reformats the data keeping the last, second..forth columns
        cat phased_reads.ctg.phased.txt | perl -lane 'print "$F[-1] $F[1] $F[2] $F[3]"' >|   phased_reads.ctg.phased.reads.reformat.txt

        mkdir -p proto

        #pulls the reference into the proto dir
        samtools faidx {input.T} {params.ctg} > proto/ref.fa

        python3 -m falcon_unzip.proto.main_augment_pb --wd ./proto/ --ctg-id {params.ctg}     --p-ctg {input.PCTG} --p-ctg-tiling-path {input.PTILE} --a-ctg {input.ACTG} --a-ctg-tiling-path {input.ATILE}  --p-variant-fn phased_vars.ctg.phased.txt --preads-sam {input.BAM}  --extracted-ctg-fasta {input.T} --rawread-bam {input.BAM}  --rid-phase-map phased_reads.ctg.phased.reads.reformat.txt  --out-updated-rid-phase_map rid_to_phase.tmp

        #grabs the names of all the CCS reads that associate with this primary contig.
        (
        set -vx +e
        grep -w {params.ctg} {input.RID_TO_CTG} >| rid_to_ctg.txt
        true
        )

        #converts the CCS read names into DAZDB read ids.
        python3 -m falcon_unzip.mains.db_to_ccs_id --lookup {input.readname_lookup} --rid-to-phase rid_to_phase.tmp --rid-to-ctg rid_to_ctg.txt --output {output.M} --ctg {params.ctg}

        python3 -m falcon_unzip.proto.extract_phased_preads --ctg-id {params.ctg} --preads ../../../../1-preads_ovl/db2falcon/preads4falcon.fasta --rid-phase-map {output.M} --out proto/preads.fasta

        time pbmm2 align --preset CCS --sort -j 1 proto/ref.fa proto/preads.fasta | samtools view  > proto/preads.sam
"""


TASK_READNAME_LOOKUP = """
       #TODO: move into a script rather than a hard to grok command line.
       #This command makes a two column file where one column is the CCS read name and the other is DAZDB ID.


       paste <(DBdump -hr ../../1-preads_ovl/build/preads.db | grep "^L" ) <(DBdump -hr ../../1-preads_ovl/build/preads.db | grep "^H" ) | perl -lane '$ln = sprintf("%09d", $. -1); @Z = split /\s+/, $_; print "$ln\t$Z[-1]/$Z[1]/ccs"' > {output.readname_lookup}
"""

TASK_MAP = """
        time pbmm2 align --sort --preset CCS -j {params.pypeflow_nproc} {input.T} {input.R} | samtools view -F 3840 -bS > {output.OBAMA}
        time samtools index {output.OBAMA}
        time samtools view -F 3844 {output.OBAMA}  | perl -lane '$F[2] =~ s/\-.*//; print "$F[0] $F[2]"' > {output.RID_TO_CTG}
"""
TASK_PREAMBLE = """
        cat {input.P} {input.A} > {output.FA}
        samtools faidx {output.FA}
        samtools faidx {input.P}
"""



def run_workflow(wf, config, unzip_config_fn):
    default_njobs = int(config['job.defaults']['njobs'])
    wf.max_jobs = default_njobs

    #LOG.info('config=\n {}'.format(config))
    Unzip_config = config['Unzip']

    falcon_asm_done_fn = './2-asm-falcon/falcon_asm_done'

    #ifastq_fn = '/home/zkronenberg/dump/hg002_chr6_dataset/hg002_chr6_size_filt.fq'
    ifastq_fn = Unzip_config['fastq']
    LOG.info('Input fastq="{}"'.format(ifastq_fn))

    asm_dir = './2-asm-falcon'
    p_ctg_fn = os.path.join(asm_dir, 'p_ctg.fa')
    a_ctg_fn = os.path.join(asm_dir, 'a_ctg.fa')
    p_tile_fn = os.path.join(asm_dir, 'p_ctg_tiling_path')
    a_tile_fn = os.path.join(asm_dir, 'a_ctg_tiling_path')

    # Typical job-dist configuration
    dist_low = Dist(
        job_dict=config['job.defaults'],
        use_tmpdir=False, # until we fix a bug in pypeflow
    )

    # high resource configuration
    dist_high = Dist(
        job_dict=config['job.high'],
        use_tmpdir=False, # until we fix a bug in pypeflow
    )

    # 1-proc jobs
    dist_one = Dist(
        NPROC = 1,
        use_tmpdir=False, # until we fix a bug in pypeflow
    )

    # 1-proc high mem jobs like htig
    dist_highmem = Dist(
        job_dict=config['job.highmem'],
        use_tmpdir=False, # until we fix a bug in pypeflow
    )

    # For strictly local jobs, use this.
    dist_local = Dist(
        local=True,
        NPROC=1,
        use_tmpdir=False, # until we fix a bug in pypeflow
    )

    p_ctg_fai_fn = "./2-asm-falcon/p_ctg.fa.fai"

    wf.addTask(gen_task(
            script=TASK_PREAMBLE,
            inputs={
                "P": p_ctg_fn,
                "A": a_ctg_fn,
            },
            outputs={
                "FA": "3-unzip/ctgs/concat.fa",
                "FAI": "3-unzip/ctgs/concat.fa.fai",
            },
            parameters={},
            dist=dist_one,
    ))

    wf.refreshTargets()

    top.fai2ctgs(p_ctg_fai_fn, 'CTGS.json')
    CTGS = io.deserialize('CTGS.json')  # currently in top-dir

    rid_to_ctg = "./3-unzip/mapping/rid_to_cgt.txt"

    wf.addTask(gen_task(
            script=TASK_MAP,
            inputs={
                "T": "3-unzip/ctgs/concat.fa",
                "R": ifastq_fn,
            },
            outputs={
                "OBAMA": "3-unzip/mapping/reads_mapped.sorted.bam",
                "OBAMSAI": "3-unzip/mapping/reads_mapped.sorted.bam.bai",
                "RID_TO_CTG" : rid_to_ctg,
            },
            parameters={},
            dist=dist_high,
    ))


    readname_lookup = "3-unzip/readnames/readname_lookup.txt"

    wf.addTask(gen_task(
            script=TASK_READNAME_LOOKUP,
            inputs={
                "OBAMSAI": "3-unzip/mapping/reads_mapped.sorted.bam.bai",
            },
            outputs={
                'readname_lookup'  : readname_lookup,
            },
            parameters={},
            dist=dist_one,
    ))

    collected = dict()

    for ctg in CTGS:
        rid_to_phase_fn = "3-unzip/0-phasing/{}/uow-fake/rid_to_phase".format(ctg)
        collected['ctg' + ctg] = rid_to_phase_fn
        wf.addTask(gen_task(
            script=TASK_READ_PHASING,
            inputs={
                "PCTG": p_ctg_fn,
                "ACTG": a_ctg_fn,
                "PTILE": p_tile_fn,
                "ATILE": a_tile_fn,
                "BAM": "3-unzip/mapping/reads_mapped.sorted.bam",
                "T": "3-unzip/ctgs/concat.fa",
                "readname_lookup" : readname_lookup,
                "RID_TO_CTG"      : rid_to_ctg,
            },
            outputs={
                "M": rid_to_phase_fn,
            },
            parameters={
                'ctg': ctg,
            },
            dist=dist_one,
        ))

    wf.refreshTargets()

    check(sorted(collected.values()), sorted(glob.glob('3-unzip/0-phasing/*/uow-fake/rid_to_phase')))

    concatenated_rid_to_phase_fn = "3-unzip/0-phasing/gathered-rid-to-phase/rid_to_phase.all"
    gathered_rid_to_phase_json   = "3-unzip/0-phasing/gathered-rid-to-phase/gathered.json"

    wf.addTask(gen_task(
        script=TASK_PHASING_GATHER_SCRIPT,
        inputs=collected,
        outputs={'rid_to_phase_all'  : concatenated_rid_to_phase_fn,
                 'gathered_rid_json' : gathered_rid_to_phase_json,
        },
        parameters={},
        dist=dist_local,
    ))

    p_las_fofn_fn =   './1-preads_ovl/las-merge-combine/las_fofn.json'
    hasm_p_ctg_fn    = './3-unzip/1-hasm/p_ctg.fa'
    preads_db_fn     = './1-preads_ovl/build/preads.db'
    preads4falcon_fn = './1-preads_ovl/db2falcon/preads4falcon.fasta'

    wf.addTask(gen_task(
            script=TASK_HASM_SCRIPT,
            inputs={
                'preads_db': preads_db_fn,
                'preads4falcon': preads4falcon_fn,
                'las_fofn': p_las_fofn_fn,
                'rid_to_phase_all': concatenated_rid_to_phase_fn,
                'gathered_rid_json' : gathered_rid_to_phase_json,
            },
            outputs={
                'p_ctg': hasm_p_ctg_fn,
            },
            parameters={},
            dist=dist_low,
    ))

    g2h_all_units_fn = './3-unzip/2-htigs/split/all-units-of-work.json'
    dummy_fn = './3-unzip/2-htigs/split/dummy.sh'
    wf.addTask(gen_task(
            script=TASK_GRAPH_TO_H_TIGS_SPLIT_SCRIPT,
            inputs={
                'falcon_asm_done': falcon_asm_done_fn,
                'preads4falcon': preads4falcon_fn,
                'rid_to_phase_all': concatenated_rid_to_phase_fn,
                'gathered_rid_json': gathered_rid_to_phase_json,
                'p_ctg': hasm_p_ctg_fn,
            },
            outputs={
                'split': g2h_all_units_fn,
                'bash_template': dummy_fn,
            },
            parameters={},
            dist = dist_highmem,
    ))

    wf.refreshTargets()

    TASK_GTOH_APPLY_UNITS_OF_WORK = """\
    python3 -m falcon_unzip.mains.graphs_to_h_tigs_2 apply --units-of-work-fn={input.units_of_work} --results-fn={output.results}

    #--bash-template-fn= # not needed
    """
    gathered_g2h_fn = './3-unzip/2-htigs/gathered/gathered.json'
    gen_parallel_tasks(
        wf,
        g2h_all_units_fn, gathered_g2h_fn,
        run_dict=dict(
            bash_template_fn=dummy_fn,
            script='DUMMY',
            inputs={
                'units_of_work': './3-unzip/2-htigs/chunks/{chunk_id}/some-units-of-work.json',
            },
            outputs={
                'results': './3-unzip/2-htigs/{chunk_id}/result-list.json',
            },
            parameters={},
        ),
        dist=dist_highmem,  # single-threads for now
        run_script=TASK_GTOH_APPLY_UNITS_OF_WORK,
    )


    job_done = './3-unzip/hasm_done'
    wf.addTask(gen_task(
            script=TASK_HASM_COLLECT_SCRIPT,
            inputs={
                'preads4falcon': preads4falcon_fn,
                'results': gathered_g2h_fn,
            },
            outputs={
                'job_done': job_done,
                'all_phased_reads': './3-unzip/all_phased_reads',
                'p_ctg_fa': './3-unzip/all_p_ctg.fa',
                'h_ctg_fa': './3-unzip/all_h_ctg.fa',
            },
            parameters={},
            dist=dist_one,
    ))

    combined_ph = "./4-polishing/input/combined_ph.fa"
    combined_ph_fai = "./4-polishing/input/combined_ph.fa.fai"
    combined_eg = "./4-polishing/input/combined_edges.txt"
    readtoctg   = "./4-polishing/input/read2ctg.txt"
    ofastq_fn   = "./4-polishing/input/preamble.fastq"

    wf.addTask(gen_task(
        script=TASK_POLISH_PREAMBLE,
        inputs={
            'P': './3-unzip/all_p_ctg.fa',
            'H': './3-unzip/all_h_ctg.fa',
            'RIDTOPHASE' : './3-unzip/0-phasing/gathered-rid-to-phase/rid_to_phase.all',
            'READNAMELOOKUP' : './3-unzip/readnames/readname_lookup.txt',
            'FQ'             : ifastq_fn,
        },
        outputs={
            'COMBINED' : combined_ph,
            'EDGES' : combined_eg,
            'READTOCTG' : readtoctg,
            "FAI"       : combined_ph_fai,
            "FQO"       : ofastq_fn,
        },
        parameters={},
        dist=dist_one,
    ))

    wf.refreshTargets()

    collected = dict()

    unzip_p_ctg_fai_fn = "3-unzip/all_p_ctg.fa.fai"
    top.fai2ctgs(unzip_p_ctg_fai_fn, 'PUCTGS.json')
    PUCTGS = io.deserialize('PUCTGS.json')  # currently in top-dir

    for ctg in PUCTGS:
        fn = '4-polishing/temp-unphased/{}/aln.bam'.format(ctg)
        collected['ctg'+ctg] = fn
        wf.addTask(gen_task(
            script=TASK_PLACE_UNPHASED,
            inputs={
                "fa"  : combined_ph,
                "fai" : combined_ph_fai,
                "fq"  : ofastq_fn,
                'readtoctg' : readtoctg,
            },
            outputs={
                "POL": fn,
            },
            parameters={
                'ctg': ctg,
            },
            dist=dist_one,
        ))

    wf.refreshTargets()

    merged_unphased = "4-polishing/merged-unphased/merged_unphased.sorted.bam"

    wf.addTask(gen_task(
        script=TASK_GATHER_UNPHASED,
        inputs=collected,
        outputs={
            "UBAM": merged_unphased,
        },
        parameters={},
        dist=dist_one,
    ))

    wf.refreshTargets()

    PH = top.getPH(combined_ph_fai, readtoctg)

    fns = list()
    LOG.info('len(PH)={}, {!r}'.format(len(PH), dist_low))
    for ctg in PH:
        fn = "4-polishing/temp-phased/{}/{}.polished.fa".format(ctg, ctg)
        fns.append(fn)
        wf.addTask(gen_task(
            script=TASK_POLISH,
            inputs={
                "FA" : combined_ph,
                "FQ" : ifastq_fn,
                "READTOCTG" : readtoctg,
                "UBAM"      : merged_unphased
            },
            outputs={
                "POL": fn,
            },
            parameters={
                'ctg': ctg,
            },
            dist=dist_low,
        ))

    wf.refreshTargets()

    # TODO: Someday this needs to be its own task.
    h_ctg_fn = 'h_ctg.polished.fa'
    p_ctg_fn = 'p_ctg.polished.fa'
    io.touch(h_ctg_fn)
    io.touch(p_ctg_fn)
    for fn in fns:
        if is_haplotig(fn):
            call = 'cat {} >> {}'.format(fn, h_ctg_fn)
        else:
            call = 'cat {} >> {}'.format(fn, p_ctg_fn)
        io.syscall(call)

def is_haplotig(fn):
    """
    >>> is_haplotig('/a/foo_bar.fa')
    True
    >>> is_haplotig('/a/foo.fa')
    False
    """
    return '_' in os.path.basename(fn)

def check(a, b):
    """Simple runtime equality checking.
    """
    if a != b:
        LOG.warning('a != b\n{!r} !=\n{!r}'.format(a, b))
    assert a == b
