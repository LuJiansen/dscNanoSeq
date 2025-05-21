ref = "GRCh38" # r mm10, GRCh38_mm10_mixed
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
database_dir = base_dir + 'database/'

if(ref in ["GRCh38","hg38"]):
    ref_dir = database_dir + 'hg38/'
elif(ref in ["mm10","GRCm38"]):
    ref_dir = database_dir + 'mm10/'
elif(ref == "GRCh38_mm10_mixed"):
    ref_dir = database_dir + 'hg38_mm10_mixed/'
else:
    print("unknown reference!")

# GM12878
snp_vcf = database_dir + 'snp/illumina_PlatinumGenomes_2017_hg38_NA12878_PS.vcf.gz'
# DBA x C57
#snp_vcf = database_dir + 'snp/mm10_SNP/DBA_2J_SNP_FI1_PS_phased.vcf.gz'

chrsize = ref_dir + ref + '.chrom.sizes'
slop_length = 75

# get wildcards automatically
samples, = glob_wildcards("mapping/{sample,[^/]+}_duprm.bam")
samples = set(samples)

rule all:
    input:
        expand("mapping/{sample}_duprm_phased_HP1_flank.bed.gz", sample = samples),
        expand("mapping/{sample}_duprm_phased_HP2_flank.bed.gz", sample = samples),
        expand("mapping/{sample}_duprm_phased_HP1_slop.bw", sample = samples),
        expand("mapping/{sample}_duprm_phased_HP2_slop.bw", sample = samples),

rule haplotag:
    input:
        bam = "mapping/{sample}_duprm.bam",
        snp = snp_vcf,
        ref = ref_dir + ref + '.fa'
    output:
        "mapping/{sample}_phased.bam",
    threads: 5,
    log:
        "logs/{sample}_haplotyping.log",
    shell: """
        if [ ! -s {input.bam}.bai ];then
            samtools index {input.bam}
        fi
        
        whatshap haplotag \
            -o {output} \
            --reference {input.ref} \
            --ignore-read-groups \
            --tag-supplementary \
            --skip-missing-contigs \
            {input.snp} {input.bam} > {log}
    """

rule split_bam_HP1:
    input:
        rules.haplotag.output,
    output:
        "mapping/{sample}_duprm_phased_HP1.bam",
    shell: """
        cat <(samtools view -H {input}) <(samtools view {input} | grep "HP:i:1") |
        samtools view -o {output} - 
    """

rule split_bam_HP2:
    input:
        rules.haplotag.output,
    output:
        "mapping/{sample}_duprm_phased_HP2.bam",
    shell: """
        cat <(samtools view -H {input}) <(samtools view {input} | grep "HP:i:2") |
        samtools view -o {output} - 
    """

rule bam2bed_HP1:
    input:
        rules.split_bam_HP1.output,
    output:
        "mapping/{sample}_duprm_phased_HP1.bed.gz",
    params:
        "mapping/{sample}_duprm_phased_HP1.bed",
    shell: """
        bedtools bamtobed -i {input} | 
        awk '{{
            split($4, name, "_")
                    print $1,$2,$3,name[2],$5,$6
        }}' OFS='\\t' > {params}
        bgzip {params}
        tabix {output}
    """

rule bam2bed_HP2:
    input:
        rules.split_bam_HP2.output,
    output:
        "mapping/{sample}_duprm_phased_HP2.bed.gz",
    params:
        "mapping/{sample}_duprm_phased_HP2.bed",
    shell: """
        bedtools bamtobed -i {input} | 
        awk '{{
            split($4, name, "_")
                    print $1,$2,$3,name[2],$5,$6
        }}' OFS='\\t' > {params}
        bgzip {params}
        tabix {output}
    """

rule HP1_flank:
    input:
        bed = rules.bam2bed_HP1.output,
        chrsize = chrsize,
    output:
        "mapping/{sample}_duprm_phased_HP1_flank.bed.gz",
    shell: """
        zcat {input.bed} | grep ^chr |
        bedtools flank -b 1 -i stdin -g {input.chrsize} |
        pigz > {output}
    """

rule HP2_flank:
    input:
        bed = rules.bam2bed_HP2.output,
        chrsize = chrsize,
    output:
        "mapping/{sample}_duprm_phased_HP2_flank.bed.gz",
    shell: """
        zcat {input.bed} | grep ^chr |
        bedtools flank -b 1 -i stdin -g {input.chrsize} |
        pigz > {output}
    """

rule HP1_slop:
    input:
        bed = rules.HP1_flank.output,
        chrsize = chrsize,
    params:
        slop_length,
    output:
        "mapping/{sample}_duprm_phased_HP1_slop.bed.gz",
    shell: """
        zcat {input.bed} | bedtools slop -b {params} -i stdin -g {input.chrsize} | pigz > {output}
    """

rule bed2bw1:
    input:
        bed = rules.HP1_slop.output,
        chrsize = chrsize,
    params:
        "mapping/{sample}_duprm_phased_HP1_slop.bedGraph",
    output:
        bg = temp("mapping/{sample}_duprm_phased_HP1_slop.bedGraph.gz"),
        bw = "mapping/{sample}_duprm_phased_HP1_slop.bw",
    shell: """
        zcat {input.bed} | sort -k1,1 -k2,2 |
        bedtools genomecov -bg -i stdin -g {input.chrsize} > {params}
        bedGraphToBigWig {params} {input.chrsize} {output.bw}
        pigz {params}
    """

rule HP2_slop:
    input:
        bed = rules.HP2_flank.output,
        chrsize = chrsize,
    params:
        slop_length,
    output:
        "mapping/{sample}_duprm_phased_HP2_slop.bed.gz",
    shell: """
        zcat {input.bed} | bedtools slop -b {params} -i stdin -g {input.chrsize} | pigz > {output}
    """

rule bed2bw2:
    input:
        bed = rules.HP2_slop.output,
        chrsize = chrsize,
    params:
        "mapping/{sample}_duprm_phased_HP2_slop.bedGraph",
    output:
        bg = temp("mapping/{sample}_duprm_phased_HP2_slop.bedGraph.gz"),
        bw = "mapping/{sample}_duprm_phased_HP2_slop.bw",
    shell: """
        zcat {input.bed} | sort -k1,1 -k2,2 |
        bedtools genomecov -bg -i stdin -g {input.chrsize} > {params}
        bedGraphToBigWig {params} {input.chrsize} {output.bw}
        pigz {params}
    """