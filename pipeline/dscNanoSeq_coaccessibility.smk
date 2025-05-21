
base_var = 'DSCNANOSEQ_DIR'
#base_dir = '/path/to/Github/dscNanoSeq/'

if base_var not in os.environ:
    print(f"Error: Environment variable {base_var} is not set", file=sys.stderr)
    print("Please set this variable before running the program", file=sys.stderr)
    sys.exit(1)  # Exit with non-zero status code indicating error

base_dir = os.environ[base_var]
print(f"{base_var} = {base_dir}")

pip_dir = base_dir + 'pipeline/'
script_dir = base_dir + 'scripts/'

# peaks should be provided with format:
#  {sample}_peaks.bed

samples, = glob_wildcards("{sample,[^/]+}_fragment.bed.gz")
samples = set(samples)
print(samples)

def ks_result(wildcards):
    checkpoints.split_fragment.get(sample = wildcards.sample)
    chrom = glob_wildcards(f"{wildcards.sample}_fragment_split/{wildcards.sample}_fragment_{{chr}}.bed.gz").chr
    return expand(f"coaccessibility/{wildcards.sample}_{{chr}}_ks_results.txt", chr = chrom)

rule all:
    input:
        expand("{sample}_coaccessibility.txt", sample = samples),

checkpoint split_fragment:
    input:
        "{sample}_fragment.bed.gz",
    output:
        directory("{sample}_fragment_split"),
    shell: """
        if [ ! -d {output} ];then
            mkdir {output}
        fi

        for i in `zcat {input} | grep ^chr | awk '{{print $1}}' | uniq`
        do
            zcat {input} | awk -v c=$i '$1==c {{print $0}}' > "{output}/{wildcards.sample}_fragment_${{i}}.bed"
        done
        pigz {output}/*bed
    """

rule run_coaccessibility:
    input:
        frag = "{sample}_fragment_split/{sample}_fragment_{chr}.bed.gz",
        peak = lambda wildcards: wildcards.sample + '_peaks.bed',
        script = script_dir + 'coaccessibility.r',
    output:
        temp("coaccessibility/{sample}_{chr}_ks_results.txt"),
    params:
        "coaccessibility/{sample}_peak_{chr}_tmp.bed",
    threads: 10,
    shell: """
        awk -v c={wildcards.chr} '$1==c {{print $0}}' OFS='\\t' {input.peak} > {params}
        Rscript {input.script} {input.frag} {params} {threads} {output}
    """

rule collect_results:
    input:
        ks_result,
    output:
        "{sample}_coaccessibility.txt",
    shell: """
        echo -e "idx\\tnpeaks_right\\tncontrol_right\\tpeak_len_right\\tcontrol_len_right\\tpval_right\\tnpeaks_left\\tncontrol_left\\tpeak_len_left\\tcontrol_len_left\\tpval_left\\tchrom\\tstart\\tend" > {output}
        cat {input} >> {output}
    """