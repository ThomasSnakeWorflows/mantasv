
import os
import shutil
from collections import defaultdict


DEFAULT_THREADS = 8
DEFAULT_MEM_PER_CPUS = 8


def get_threads(rule, default=DEFAULT_THREADS):
    cluster_config = snakemake.workflow.cluster_config
    if rule in cluster_config and "threads" in cluster_config[rule]:
        return cluster_config[rule]["threads"]
    if "default" in cluster_config and "threads" in cluster_config["default"]:
        return cluster_config["default"]["threads"]
    return default


def get_mem(rule, default):
    cluster_config = snakemake.workflow.cluster_config
    if rule in cluster_config and "mem" in cluster_config[rule]:
        return cluster_config[rule]["mem"]
    if "default" in cluster_config and "mem" in cluster_config["default"]:
        return cluster_config["default"]["mem"]
    return default


def get_manta_mem(cpus):
    cluster_config = snakemake.workflow.cluster_config
    if "runmanta" in cluster_config and "mem-per-cpu" in cluster_config["runmanta"]:
        mem = int(cpus) * int(cluster_config["runmanta"]["mem-per-cpu"])
        return mem
    if "default" in cluster_config and "mem-per-cpu" in cluster_config["default"]:
        mem = int(cpus) * int(cluster_config["default"]["mem-per-cpu"])
        return mem
    return cpus*DEFAULT_MEM_PER_CPUS


def get_bams(sample_file):
    with open(sample_file) as fin:
        bams = [os.path.abspath(line.rstrip()) for line in fin]
    return bams


bamfiles = get_bams(config['sample_file'])
reference = os.path.abspath(config['reference'])

rootmanta = os.environ['MANTA_ROOT']
manta_bin = ":".join([os.path.abspath(rootmanta + "/bin"),
                      os.path.abspath(rootmanta + "/libexec")])
os.environ["PATH"] = "%s:%s" % (manta_bin, os.environ["PATH"])

workdir: config['workdir']

svtypes = ['DEL', 'INS', 'INV', 'DUP']

localrules: regions, configmanta, combine, convert, splitvcf

rule all:
    input:
        expand("manta_{svtype}.vcf.gz", svtype=svtypes),
        expand("manta_{svtype}.vcf.gz.tbi", svtype=svtypes)

rule convert:
    input:
        vcf = "rundir/results/variants/diploidSV.vcf.gz",
        reference = reference,
    output:
        "rundir/results/variants/diploidSVconverted.vcf.gz"
    params:
        convertinv = "",
        samtoolsexe = shutil.which("samtools")
    shell:
        "convertInversion.py {params.samtoolsexe} {input.reference} {input.vcf}"
        " | bcftools sort -Oz -o {output}"


rule configmanta:
    input:
        bams = bamfiles,
        reference = reference,
        fai = reference+".fai",
        regions = "regions.bed.gz"
    output:
        mantarunner = "rundir/runWorkflow.py"
    params:
        rundir = "rundir"
    log:
        stdout = "logs/config.o",
        stderr = "logs/config.e"
    threads:
        1
    run:
        from subprocess import run
        command = "configManta.py "
        for bam in input.bams:
            command += "--bam %s " % bam
        command += "--referenceFasta %s " % input.reference
        command += "--runDir %s " % params.rundir
        command += "--callRegions %s " % input.regions
        print(command)
        result = run(command, shell=True)

rule regions:
    input:
        reference = reference,
        fai = reference+".fai"
    output:
        "regions.bed.gz.tbi",
        region = "regions.bed.gz"
    params:
        chromregex = config['chromregex']
    shell:
        "cat {input.fai} | egrep \'^{params.chromregex}\' "
        " | awk -v OFS='\\t' '{{ print $1,0,$2}}' | bgzip -c > {output.region}; "
        "tabix {output.region}"

rule runmanta:
    input:
        mantarunner = "rundir/runWorkflow.py"
    output:
        "rundir/results/variants/diploidSV.vcf.gz"
    threads:
        get_threads("runmanta")
    params:
        mem = get_manta_mem(get_threads("runmanta"))
    log:
        stdout = "logs/run.o",
        stderr = "logs/run.e"
    shell:
        " python2 {input} -j {threads} --memGb={params.mem} --quiet "
        " 1>{log.stdout} 2>{log.stderr} "


rule combine:
    input:
        diploid = "rundir/results/variants/diploidSVconverted.vcf.gz"
    output:
        vcf="manta.vcf.gz",
        tbi="manta.vcf.gz.tbi"
    params:
        resultdir = "rundir/results"
    shell:
        """
        bcftools view {input.diploid} -Oz -o {output.vcf}
        tabix {output.vcf}
        """
        #"../manta_combine.sh {params.resultdir} | bgzip -c > {output}"

rule splitvcf:
    input:
        "manta.vcf.gz"
    output:
        variants = "manta_{svtype}.vcf.gz",
        index = "manta_{svtype}.vcf.gz.tbi"
    shell:
        """
        bcftools view -i'INFO/SVTYPE=\"{wildcards.svtype}\"' {input} -Oz -o {output.variants}
        tabix {output.variants}
        """
