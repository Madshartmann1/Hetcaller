# Snakefile
# Streaming + thresholds expansion + ROH/ROH_trv + R stub
# Run with:  hetcall --help  (wrapper below)
# Or:        snakemake -j {cores} --configfile config.yaml

shell.executable("/bin/bash")
import os, shlex, ast

def _get_thresholds(x):
    if isinstance(x, (list, tuple)):
        return [str(t) for t in x]
    return [s for s in (y.strip() for y in str(x).strip("[]").split(",")) if s]

PREFIX   = config.get("prefix", "sample01")
OUTDIR   = config.get("outdir", "results")
SCRIPTS  = config.get("scripts", "scripts")
MINDEPTH = int(config.get("min_depth", 10))
THRESHOLDS = _get_thresholds(config.get("thresholds", [0.05]))
ROH_MIN  = config.get("roh_min", 0.2)
ROH_MIN_STR = str(ROH_MIN)
ANGSD_THREADS = int(config.get("angsd_threads", 8))

# Prefer site R if available, allow override from config
RSCRIPT = config.get("Rscript", "/home/ctools/R-4.1.2/bin/Rscript" if os.path.exists("/home/ctools/R-4.1.2/bin/Rscript") else "Rscript")
# Preferred plot format (pdf or svg)
PLOT_FORMAT = config.get("plot_format", "pdf")

RF  = str(config.get("rf", "") or "")
REG = str(config.get("regions", "") or "")

def _coerce_rf_reg(rf: str, reg: str):
    # Clean up any quotes that might have been added
    rf = rf.strip().strip('"').strip("'")
    reg = reg.strip().strip('"').strip("'")
    
    # If rf looks like a region string, treat it as regions
    if rf and any(t in rf for t in (":", "-", ",", ";", "|")) and not reg:
        reg = rf
        rf = ""
    
    return rf, reg

RF, REG = _coerce_rf_reg(RF, REG)
RF_FLAG = f"-rf {shlex.quote(RF)}" if RF else (f"-r {shlex.quote(REG)}" if REG else "")

pos_gz    = f"{OUTDIR}/{PREFIX}.pos.gz"
counts_gz = f"{OUTDIR}/{PREFIX}.counts.gz"

def basecall_txt(th):
    return f"{OUTDIR}/{PREFIX}.basecall_{th}_mincov{MINDEPTH}.txt"

def roh_txt(th, trv=False):
    tag = "_trv" if trv else ""
    return f"{OUTDIR}/{PREFIX}_het{th}_thres{config['roh_min']}.ROH{tag}.txt"

def minorfreq_all():
    return [
        f"{OUTDIR}/{PREFIX}.minorfreq.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_AC.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_AG.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_AT.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_CG.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_CT.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_GT.txt",
    ]


rule all:
    input:
        pos_gz,
        counts_gz,
        expand(f"{OUTDIR}/{PREFIX}.basecall_{{thresh}}_mincov{MINDEPTH}.txt.gz", thresh=THRESHOLDS),
        expand(f"{OUTDIR}/{PREFIX}_het{{thresh}}_thres{ROH_MIN_STR}.ROH.txt.gz",     thresh=THRESHOLDS),
        expand(f"{OUTDIR}/{PREFIX}_het{{thresh}}_thres{ROH_MIN_STR}.ROH_trv.txt.gz", thresh=THRESHOLDS),
        f"{OUTDIR}/{PREFIX}.minorfreq.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_AC.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_AG.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_AT.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_CG.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_CT.txt",
        f"{OUTDIR}/{PREFIX}.minorfreq_GT.txt",
        f"{OUTDIR}/plots/{PREFIX}.plots.ok"


rule angsd_counts:
    output:
        pos_gz    = pos_gz,
        counts_gz = counts_gz
    threads: ANGSD_THREADS
    params:
        bam=config["bam"],
        outpref=f"{OUTDIR}/{PREFIX}",
        rf_flag=RF_FLAG
    log:
        f"{OUTDIR}/{PREFIX}.angsd.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}"
        /home/ctools/angsd-0.935/angsd -minq 20 -minmapq 20 -uniqueOnly 1 -remove_bads 1 \
              -docounts 1 -dumpCounts 4 \
              -i "{params.bam}" -out "{params.outpref}" \
              -setMinDepthInd {MINDEPTH} {params.rf_flag} -nThreads {threads} \
              > "{log}" 2>&1
        test -s "{output.pos_gz}" && test -s "{output.counts_gz}"
        """

rule basecall_txt:
    input:
        pos = pos_gz,
        counts = counts_gz
    output:
        txt = f"{OUTDIR}/{PREFIX}.basecall_{{thresh}}_mincov{MINDEPTH}.txt.gz"
    params:
        scripts = SCRIPTS
    shell:
        r"""
        set -euo pipefail
        paste \
          <( zcat "{input.pos}" ) \
          <( cat "{params.scripts}/header.txt"; \
             zcat "{input.counts}" | tail -n +2 | sh "{params.scripts}/Basecall.sh" - {wildcards.thresh} ) \
          | gzip -c > "{output.txt}"
        """

rule minorfreq:
    input:
        base=f"{OUTDIR}/{PREFIX}.basecall_{THRESHOLDS[0]}_mincov{MINDEPTH}.txt.gz"
    output:
        all  = f"{OUTDIR}/{PREFIX}.minorfreq.txt",
        AC   = f"{OUTDIR}/{PREFIX}.minorfreq_AC.txt",
        AG   = f"{OUTDIR}/{PREFIX}.minorfreq_AG.txt",
        AT   = f"{OUTDIR}/{PREFIX}.minorfreq_AT.txt",
        CG   = f"{OUTDIR}/{PREFIX}.minorfreq_CG.txt",
        CT   = f"{OUTDIR}/{PREFIX}.minorfreq_CT.txt",
        GT   = f"{OUTDIR}/{PREFIX}.minorfreq_GT.txt"
    shell:
        r"""
        set -euo pipefail
        # replicate your awk/grep pipeline, one pass
        zcat "{input.base}" | grep HET \
        | awk '{{split($0, arr); delete arr[1]; delete arr[2]; delete arr[3]; delete arr[8]; asort(arr); print arr[length(arr)-2],$8}}' \
        | sort | uniq -c > "{output.all}"

        grep -w AC "{output.all}" > "{output.AC}" || true
        grep -w AG "{output.all}" > "{output.AG}" || true
        grep -w AT "{output.all}" > "{output.AT}" || true
        grep -w CG "{output.all}" > "{output.CG}" || true
        grep -w CT "{output.all}" > "{output.CT}" || true
        grep -w GT "{output.all}" > "{output.GT}" || true
        """

rule roh:
    input:
        base = f"{OUTDIR}/{PREFIX}.basecall_{{thresh}}_mincov{MINDEPTH}.txt.gz"
    output:
        txt = f"{OUTDIR}/{PREFIX}_het{{thresh}}_thres{config['roh_min']}.ROH.txt.gz"
    params:
        scripts = SCRIPTS, win = 100000, step = 25000, thr_min = config["roh_min"]
    shell:
        r"""
        set -euo pipefail
        {{
          zcat "{input.base}" | python "{params.scripts}/ROH.py" /dev/stdin {params.win} {params.step} {params.thr_min}
        }} > /tmp/roh_$$.txt
        {{
          cat /tmp/roh_$$.txt
          grep "consecutive windows" /tmp/roh_$$.txt | awk '$1>=1{{print $1*$NF}}'  | awk '{{sum+=$1}}END{{print "BP in ROH at least 100kb: "sum*100000}}'
          grep "consecutive windows" /tmp/roh_$$.txt | awk '$1>=10{{print $1*$NF}}' | awk '{{sum+=$1}}END{{print "BP in ROH at least 1Mb: "sum*100000}}'
          grep "consecutive windows" /tmp/roh_$$.txt | awk '$1>=50{{print $1*$NF}}' | awk '{{sum+=$1}}END{{print "BP in ROH at least 5Mb: "sum*100000}}'
          echo "BP analysed $(grep -c Proportion /tmp/roh_$$.txt | awk '{{print $1*100000}}')"
          echo "Genome-wide Heterozygosity $(zcat "{input.base}" | grep -c HET | paste - <(zcat "{input.base}" | wc -l) | awk '{{print $1/($2-1)}}')"
        }} | gzip -c > "{output.txt}"
        rm -f /tmp/roh_$$.txt
        """

#transcversions only    
rule roh_trv:
    input:
        base=f"{OUTDIR}/{PREFIX}.basecall_{{thresh}}_mincov{MINDEPTH}.txt.gz"
    output:
        txt=f"{OUTDIR}/{PREFIX}_het{{thresh}}_thres{config['roh_min']}.ROH_trv.txt.gz"
    params:
        scripts=SCRIPTS, win=100000, step=25000, thr_min=config["roh_min"]
    shell:
        r"""
        set -euo pipefail
        {{
          zcat "{input.base}" | python "{params.scripts}/ROH_trv.py" /dev/stdin {params.win} {params.step} {params.thr_min}
        }} > /tmp/roh_trv_$$.txt
        {{
          cat /tmp/roh_trv_$$.txt
          grep "consecutive windows" /tmp/roh_trv_$$.txt | awk '$1>=1{{print $1*$NF}}'  | awk '{{sum+=$1}}END{{print "BP in ROH at least 100kb: "sum*100000}}'
          grep "consecutive windows" /tmp/roh_trv_$$.txt | awk '$1>=10{{print $1*$NF}}' | awk '{{sum+=$1}}END{{print "BP in ROH at least 1Mb: "sum*100000}}'
          grep "consecutive windows" /tmp/roh_trv_$$.txt | awk '$1>=50{{print $1*$NF}}' | awk '{{sum+=$1}}END{{print "BP in ROH at least 5Mb: "sum*100000}}'
          echo "BP analysed $(grep -c Proportion /tmp/roh_trv_$$.txt | awk '{{print $1*100000}}')"
          echo "Genome-wide Heterozygosity $(zcat "{input.base}" | grep -c HET_trv | paste - <(zcat "{input.base}" | wc -l) | awk '{{print $1/($2-1)}}')"
        }} | gzip -c > "{output.txt}"
        rm -f /tmp/roh_trv_$$.txt
        """


rule plot_qc:
    input:
        # put real plot inputs later; this keeps a dependency chain
        minorfreq_all()
    output:
        plots_ok = f"{OUTDIR}/plots/{PREFIX}.plots.ok",
        plot_file = f"{OUTDIR}/plots/{PREFIX}_minorfreq_plots.{PLOT_FORMAT}",
        summary_txt = f"{OUTDIR}/plots/{PREFIX}_minorfreq_summary.txt"
    params:
        scripts = SCRIPTS,
        minorfreq_file = f"{OUTDIR}/{PREFIX}.minorfreq.txt",
        output_prefix = f"{OUTDIR}/plots/{PREFIX}",
        rscript = RSCRIPT,
        plot_format = PLOT_FORMAT
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/plots"
        "{params.rscript}" "{params.scripts}/plot_minorfreq.R" "{params.minorfreq_file}" --output "{params.output_prefix}" --format "{params.plot_format}"
        touch "{output.plots_ok}"
        """

