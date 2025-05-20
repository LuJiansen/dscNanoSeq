ref = "GRCh38" # r mm10, GRCh38_mm10_mixed
base_dir = '/mnt/d/GitHub/dscNanoSeq/'
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

bc_regex = '(?P<cell_1>.{10}){e<=2}(?P<discard_1>.{0,3})(?P<discard_2>TCGAG){e<=1}(?P<cell_2>.{10}){e<=2}(?P<umi_1>.{10})(?P<discard_3>CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG){e<=7}.*'
#bc_regex = '(?P<cell_1>.{{10}}){{e<=2}}(?P<discard_1>.{{0,3}})(?P<discard_2>TCGAG){{e<=1}}(?P<cell_2>.{{10}}){{e<=2}}(?P<umi_1>.{{10}})(?P<discard_3>CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG){{e<=7}}.*'

whitelist = database_dir + "mobi_whitelist_v2.txt"
mmi = ref_dir + ref + '.mmi'
chrsize = ref_dir + ref + '.chrom.sizes'
clip_length = 150
slop_length = 75
chunk_size = 500000

print(bc_regex)
# get wildcards automatically
samples, = glob_wildcards("raw_data/{sample,[^/]+}.fastq.gz")
samples = set(samples)

def filter_bams(wildcards):
    checkpoints.split_fastq.get(sample = wildcards.sample)
    chunks = glob_wildcards(f"raw_data/split/{wildcards.sample}/{wildcards.sample}.part_{{chunk}}.fastq.gz").chunk
    return expand(f"mapping/{wildcards.sample}_{{chunk}}_bc_q30.bam", chunk = chunks)

def merge_fastqs(wildcards):
    checkpoints.split_fastq.get(sample = wildcards.sample)
    chunks = glob_wildcards(f"raw_data/split/{wildcards.sample}/{wildcards.sample}.part_{{chunk}}.fastq.gz").chunk
    return expand(f"trim/{wildcards.sample}_{{chunk}}_bc_merged.fastq.gz", chunk = chunks)

def merge_fastq_stats(wildcards):
    checkpoints.split_fastq.get(sample = wildcards.sample)
    chunks = glob_wildcards(f"raw_data/split/{wildcards.sample}/{wildcards.sample}.part_{{chunk}}.fastq.gz").chunk
    return expand(f"trim/{wildcards.sample}_{{chunk}}_fastq_stats.txt", chunk = chunks)

def merge_map_stats(wildcards):
    checkpoints.split_fastq.get(sample = wildcards.sample)
    chunks = glob_wildcards(f"raw_data/split/{wildcards.sample}/{wildcards.sample}.part_{{chunk}}.fastq.gz").chunk
    return expand(f"mapping/{wildcards.sample}_{{chunk}}_bc.bam.stats", chunk = chunks)

def merge_filter_stats(wildcards):
    checkpoints.split_fastq.get(sample = wildcards.sample)
    chunks = glob_wildcards(f"raw_data/split/{wildcards.sample}/{wildcards.sample}.part_{{chunk}}.fastq.gz").chunk
    return expand(f"mapping/{wildcards.sample}_{{chunk}}_bc_q30.bam.stats", chunk = chunks)


rule all:
    input:
        # expand("bigwig/{sample}_flank.bw", sample = samples),
        # expand("bigwig/{sample}_slop.bw", sample = samples),
        expand("mapping/{sample}_duprm.bam", sample = samples),
        expand("arrow/{sample}.arrow", sample = samples),
        expand("trim/{sample}_fastq_stats.txt", sample = samples),
        expand("mapping/{sample}_bc.bam.stats", sample = samples),
        expand("mapping/{sample}_bc_q30.bam.stats", sample = samples),
        expand("mapping/{sample}_dedup.bam.stats", sample = samples),
        expand("trim/{sample}_all_bc_merged.fastq.gz", sample = samples),
        expand("trim/{sample}_raw_barcode_stats.txt", sample = samples),

# trim read1 sequences
checkpoint split_fastq:
    input:
        "raw_data/{sample}.fastq.gz",
    output:
        directory('raw_data/split/{sample}'),
    params:
        size = chunk_size,
    threads: 20, 
    shell: """
        seqkit split2 -s {params.size} -t dna -j {threads} -1 {input} -O {output}
    """

rule trim:
    input:
        "raw_data/split/{sample}/{sample}.part_{chunk}.fastq.gz",
    output:
        temp("trim/split/{sample}/{sample}.part_{chunk}_trimmed.fastq.gz"),
    threads: 5,
    log:
        "logs/{sample}_{chunk}_cutadpt.log",
    shell: """
        cutadapt \
            -q 7 -e 0.25 -j {threads} -m 100 \
            -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT...AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -o {output} \
            {input} \
            > {log}
    """

# demultiplex barcodes and UMIs in 5'
rule umitools_run1:
    input:
        fq = "trim/split/{sample}/{sample}.part_{chunk}_trimmed.fastq.gz",
        whitelist = whitelist,
    output:
        tagged = temp("trim/{sample}_{chunk}_tagged_1.fastq.gz"),
        filtered = temp("trim/{sample}_{chunk}_filtered_1.fastq.gz"),
    params:
        lambda wildcards: bc_regex,
    log:
        "logs/{sample}_{chunk}_umitools_run1.log",
    shell: """
        # umi_tools version 1.1.5
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern='{params}' \
            --whitelist={input.whitelist} \
            --filtered-out {output.filtered} \
            -I {input.fq} \
            -L {log} -S {output.tagged}
    """

# make reverse-complement of the run1 results,
# Since `umi_tools extract` is very slow when processing 3' barcodes 
rule reverse_comp_tagged:
    input:
        rules.umitools_run1.output.tagged,
    output:
        temp("trim/{sample}_{chunk}_tagged_rc_1.fastq.gz"),
    shell: """
        seqkit seq -r -p {input} | pigz > {output}
    """

rule reverse_comp_filtered:
    input:
        rules.umitools_run1.output.filtered,
    output:
        temp("trim/{sample}_{chunk}_filtered_rc_1.fastq.gz"),
    shell: """
        seqkit seq -r -p {input} | pigz > {output}
    """

# extract the 3' barcodes and UMIs
rule umitools_run2_tagged:
    input:
        fq = rules.reverse_comp_tagged.output,
        whitelist = whitelist,
    output:
        tagged=temp("trim/{sample}_{chunk}_tagged_rc_1_tagged_2.fastq.gz"),
        filtered=temp("trim/{sample}_{chunk}_tagged_rc_1_filtered_2.fastq.gz"),
    params:
        lambda wildcards: bc_regex,
    log:
        "logs/{sample}_{chunk}_umitools_run2_tagged.log",
    shell: """
        # umi_tools version 1.1.5
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern='{params}' \
            --whitelist={input.whitelist} \
            --filtered-out {output.filtered} \
            -I {input.fq} \
            -L {log} -S {output.tagged}
    """

rule umitools_run2_filtered:
    input:
        fq = rules.reverse_comp_filtered.output,
        whitelist = whitelist,
    output:
        tagged=temp("trim/{sample}_{chunk}_filtered_rc_1_tagged_2.fastq.gz"),
        filtered=temp("trim/{sample}_{chunk}_filtered_rc_1_filtered_2.fastq.gz"),
    params:
        lambda wildcards: bc_regex,
    log:
        "logs/{sample}_{chunk}_umitools_run2_filtered.log",
    shell: """
        # umi_tools version 1.1.5
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern='{params}' \
            --whitelist={input.whitelist} \
            --filtered-out {output.filtered} \
            -I {input.fq} \
            -L {log} -S {output.tagged}
    """

# For reads detected dual barcodes, we evaluate their accordance
rule dual_barcode_compare:
    input:
        rules.umitools_run2_tagged.output.tagged,
    output:
        passed = temp("trim/{sample}_{chunk}_tagged_1_tagged_2_accor.fastq.gz"),
        filtered = temp("trim/{sample}_{chunk}_tagged_1_tagged_2_discor.fastq.gz"),
    params:
        "trim/{sample}_{chunk}_tagged_1_tagged_2_discor.fastq",
    shell: """
        zcat {input} | awk -v filter={params} '{{
            if (NR % 4 == 1) {{
                split($1, name, " ")
                split(name[1], parts, "_")
                if (parts[2] == parts[4]) {{
                    print
                    getline; print
                    getline; print
                    getline; print
                }} else {{
                    print > filter
                    getline; print > filter
                    getline; print > filter
                    getline; print > filter
                }}
            }}
        }}' | pigz > {output.passed}
        pigz {params}
    """

# collect all the reads that:
# - have dual accordant barcodes
# - have run1 but not run2 barcodes
# - have run2 but not run1 barcodes
rule collect_reads:
    input:
        accor = rules.dual_barcode_compare.output.passed,
        run1_pass = rules.umitools_run2_tagged.output.filtered,
        run2_pass = rules.umitools_run2_filtered.output.tagged,
    output:
        temp("trim/{sample}_{chunk}_bc_merged.fastq.gz"),
    shell: """
        cat {input.accor} {input.run1_pass} {input.run2_pass} |
        seqkit seq -r -p | pigz > {output}
    """

rule fastq_stats:
    input:
        raw = "raw_data/split/{sample}/{sample}.part_{chunk}.fastq.gz",
        trimmed = rules.trim.output,
        run1_pass = rules.umitools_run1.output.tagged,
        run1_fail = rules.umitools_run1.output.filtered,
        run2_pass = rules.umitools_run2_filtered.output.tagged,
        run2_fail = rules.umitools_run2_filtered.output.filtered,
        accor = rules.dual_barcode_compare.output.passed,
        disaccor = rules.dual_barcode_compare.output.filtered,
        merged = rules.collect_reads.output,
    output:
        temp("trim/{sample}_{chunk}_fastq_stats.txt"),
    shell: """
        seqkit stats -a -T -j {threads} {input} > {output}
    """

rule merge_fastqs:
    input:
        merge_fastqs,
    output:
        "trim/{sample}_all_bc_merged.fastq.gz",
    shell: """
        cat {input} > {output}
    """

rule raw_barcode_sum:
    input:
        rules.merge_fastqs.output,
    output:
        "trim/{sample}_raw_barcode_stats.txt",
    shell: """
        seqkit seq -n {input} | awk -v FS='_' '{{print $2}}' | sort | uniq -c | awk '{{print $2,$1}}' OFS='\\t' > {output}
    """

rule collect_fastq_stats:
    input:
        merge_fastq_stats,
    output:
        "trim/{sample}_fastq_stats.txt",
    shell: """
        echo -e "file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\tQ1\\tQ2\\tQ3\\tsum_gap\\tN50\\tQ20(%)\\tQ30(%)" > {output}
        cat {input} | grep -P -v '^file\\tformat' >> {output}
    """

rule minimap2:
    input:
        fq = rules.collect_reads.output,
        ref = mmi,
    output:
        temp("mapping/{sample}_{chunk}_mm2.bam"),
    log:
        "logs/{sample}_{chunk}_mm2.log",
    threads: 20,
    shell: """
        minimap2 -t {threads} \
            --split-prefix mapping/{wildcards.sample}_{wildcards.chunk} \
            -ax map-ont \
            --secondary=no \
            {input.ref} \
            {input.fq} 2>{log} |
        samtools view -o {output} -
        """

# add barcodes to bam tag
rule add_barcode:
    input:
        rules.minimap2.output,
    output:
        temp("mapping/{sample}_{chunk}_bc.bam"),
    shell: """
        cat <(samtools view -H {input}) \
            <(samtools view {input} |
            awk -vOFS='\\t' '{{
                    split($1, name, "_")
                    print $0,"CB:Z:"name[2],"RX:Z:"name[3]
            }}') |
        samtools view -bS - > {output}
    """

rule map_stats:
    input:
        rules.add_barcode.output,
    output:
        temp("mapping/{sample}_{chunk}_bc.bam.stats"),
    shell: """
        bedtools bamtobed -i {input} | 
        awk '{{
            split($4, name, "_")
                    print $1,$2,$3,name[2],$5,$6
        }}' OFS='\\t' |
        awk '{{
            if($5>=20){{q20=1}}else{{q20=0}};if($5>30){{q30=1}}else{{q30=0}}
        }};{{
            print $4,$3-$2,q20,q30
        }}' OFS='\\t' | sort -k1 | 
        bedtools groupby -g 1 -c 1,2,2,2,2,2,3,4 -o count,sum,min,max,mean,median,sum,sum |
        awk '{{print $1,$2,$3,$4,$5,$6,$7,$8/$2,$9/$2}}' OFS='\\t' > {output}
    """

rule collect_map_stats:
    input:
        stats = merge_map_stats,
        script = script_dir + 'merge_bam_stats.r',
    output:
        "mapping/{sample}_bc.bam.stats",
    params:
        "_bc.bam.stats",
    shell: """
        Rscript {input.script} "mapping" {params} {output}
    """

# filter alignments, only unique mapped reads with clipping length < threshold were remained.
rule filter:
    input:
        rules.add_barcode.output,
    output:
        temp("mapping/{sample}_{chunk}_bc_q30.bam"),
    params:
        clip = clip_length,
        tmp = "tmp/{sample}_{chunk}",
    shell: """
        if [ ! -d tmp ];then
            mkdir tmp
        fi 
        samtools view -ShuF 2308 -q 30 -e "sclen < {params.clip} && hclen < {params.clip}" {input} |
        samtools sort -T {params.tmp} -o {output} - 
        """

rule filter_stats:
    input:
        rules.filter.output,
    output:
        temp("mapping/{sample}_{chunk}_bc_q30.bam.stats"),
    shell: """
        bedtools bamtobed -i {input} | 
        awk '{{
            split($4, name, "_")
                    print $1,$2,$3,name[2],$5,$6
        }}' OFS='\\t' |
        awk '{{
            if($5>=20){{q20=1}}else{{q20=0}};if($5>30){{q30=1}}else{{q30=0}}
        }};{{
            print $4,$3-$2,q20,q30
        }}' OFS='\\t' | sort -k1 | 
        bedtools groupby -g 1 -c 1,2,2,2,2,2,3,4 -o count,sum,min,max,mean,median,sum,sum |
        awk '{{print $1,$2,$3,$4,$5,$6,$7,$8/$2,$9/$2}}' OFS='\\t' > {output}
    """

rule collect_filter_stats:
    input:
        stats = merge_filter_stats,
        script = script_dir + 'merge_bam_stats.r',
    output:
        "mapping/{sample}_bc_q30.bam.stats",
    params:
        "_bc_q30.bam.stats",
    shell: """
        Rscript {input.script} "mapping" {params} {output}
    """

rule merge_bams:
    input:
        filter_bams,
    output:
        list = "mapping/{sample}_bam_list",
        bam = temp("mapping/{sample}_bc_q30.bam"),
    threads: 20,
    shell: """
        ls {input} > {output.list}
        samtools merge -@ {threads} -b {output.list} -o {output.bam}
    """

# dedup with umi_tools, since UMIs tend to be overestimated due to the sequencing 
# errors of Nanopore, and each genomic loci usaually was single copy in scATAC,
# so we only consider the mapped coordinates of alignments when deduplication.
rule umi_tools_dedup:
    input:
        rules.merge_bams.output.bam,
    output:
        "mapping/{sample}_duprm.bam",
    log:
        "logs/{sample}_dedup.log",
    shell: """
        # umi_tools version 1.1.5
        if [ ! -s {input}.bai ];then
            samtools index {input}
        fi

        umi_tools dedup \
            --extract-umi-method=tag \
            --per-cell \
            --ignore-umi \
            --cell-tag=CB \
            --temp-dir=mapping \
            -I {input} \
            -L {log} -S {output}
    """

# chr start end cellname MAPQ strand
rule bamtobed:
    input:
        rules.umi_tools_dedup.output,
    output:
        "fragment/{sample}_fragment.bed",
    shell: """
        bedtools bamtobed -i {input} | 
        awk '{{
            split($4, name, "_")
                    print $1,$2,$3,name[2],$5,$6
        }}' OFS='\\t' > {output}
    """

rule bgzip_frag:
    input:
        bed = rules.bamtobed.output,
    output:
        gz = "fragment/{sample}_fragment.bed.gz",
        idx = "fragment/{sample}_fragment.bed.gz.tbi",
    shell: """
        bgzip {input}
        tabix {output.gz}
    """

rule dedup_stats:
    input:
        rules.bgzip_frag.output.gz,
    output:
        "mapping/{sample}_dedup.bam.stats",
    shell: """
        echo -e "barcode\\tnFrags\\tdatasize\\tmin_len\\tmax_len\\tmean_len\\tmedian_len\\tMAQ20\\tMAQ30" > {output}
        zcat {input} |
        awk '{{
            if($5>=20){{q20=1}}else{{q20=0}};if($5>30){{q30=1}}else{{q30=0}}
        }};{{
            print $4,$3-$2,q20,q30
        }}' OFS='\\t' | sort -k1 | 
        bedtools groupby -g 1 -c 1,2,2,2,2,2,3,4 -o count,sum,min,max,mean,median,sum,sum |
        awk '{{print $1,$2,$3,$4,$5,$6,$7,$8/$2,$9/$2}}' OFS='\\t' >> {output}
    """

rule create_arrow:
    input:
        frag = rules.bgzip_frag.output.gz,
        script = script_dir + 'create_arrow.r',
    output:
        "arrow/{sample}.arrow",
    params:
        ref = ref,
    shell: """
        if [ ! -d arrow ];then
            mkdir -p arrow
        fi
        
        Rscript {input.script} {input.frag} {params.ref} {wildcards.sample}
        mv "{wildcards.sample}.arrow" arrow
    """

rule flank:
    input:
        bed = rules.bgzip_frag.output.gz,
        chrsize = chrsize,
    output:
        "fragment/{sample}_flank.bed",
    shell: """
        zcat {input.bed} | grep ^chr |
        bedtools flank -b 1 -i stdin -g {input.chrsize} > {output}
    """

rule bed2bw:
    input:
        bed = rules.flank.output,
        chrsize = chrsize,
    params:
        "bigwig/{sample}_flank.bedGraph",
    output:
        bg = temp("bigwig/{sample}_flank.bedGraph.gz"),
        bw = "bigwig/{sample}_flank.bw",
    shell: """
        cat {input.bed} | sort -k1,1 -k2,2 |
        bedtools genomecov -bg -i stdin -g {input.chrsize} > {params}
        bedGraphToBigWig {params} {input.chrsize} {output.bw}
        pigz {params}
    """

rule slop:
    input:
        bed = rules.flank.output,
        chrsize = chrsize,
    params:
        slop_length,
    output:
        "fragment/{sample}_slop.bed",
    shell: """
        bedtools slop -b {params} -i {input.bed} -g {input.chrsize} > {output}
    """

rule bed2bw2:
    input:
        bed = rules.slop.output,
        chrsize = chrsize,
    params:
        "bigwig/{sample}_slop.bedGraph",
    output:
        bg = temp("bigwig/{sample}_slop.bedGraph.gz"),
        bw = "bigwig/{sample}_slop.bw",
    shell: """
        cat {input.bed} | sort -k1,1 -k2,2 |
        bedtools genomecov -bg -i stdin -g {input.chrsize} > {params}
        bedGraphToBigWig {params} {input.chrsize} {output.bw}
        pigz {params}
    """