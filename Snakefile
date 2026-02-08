import os

FASTAS = config["fastas"].split()
REFERENCE = config["reference"]

# assign names automatically
SAMPLES = {}
for f in FASTAS:
    name = os.path.splitext(os.path.basename(f))[0]
    SAMPLES[name] = f

assert REFERENCE in SAMPLES

rule make_pathfile:
    output:
        "pathfile.txt"
    run:
        with open(output[0], "w") as f:
            # simple star tree rooted at reference
            leaves = ",".join(SAMPLES.keys())
            f.write(f"({leaves});\n")

            for name, fasta in SAMPLES.items():
                f.write(f"{name} {os.path.abspath(fasta)}\n")

rule cactus:
    input:
        pathfile="pathfile.txt"
    output:
        gfa="results/pan.gfa.gz"
    params:
        ref=REFERENCE
    shell:
        """
        cactus-pangenome work {input.pathfile} \
          --outDir results \
          --outName pan \
          --reference {params.ref};
        rm {input.pathfile}
        """

rule scan:
    input:
        ref_file=REFERENCE
    output:
        bed="results/ref-tan.bed"
    shell:
        "tantan -f4 {input.ref_file} > {output.bed}"

rule unzip:
    input:
        gfa_zipped="results/pan.gfa.gz"
    output:
        gfa_unziped="results/pan.gfa"
    shell:
        "gunzip {input.gfa_zipped}}"

rule analyze:
    input:
        gfa="results/pan.gfa",
        tandems="results/ref-tan.bed"
    output:
        log="results/repeats_aligned.txt"
    shell:
        """
        python processing.py \
          --graph {input.gfa} \
          --repeats {input.tandems} \
          --log {output.log} \
        """

rule all:
    output:
        "results/repeats_aligned.txt"