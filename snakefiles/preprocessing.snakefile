
indir = config["indir"]
outdir = config["outdir"]
sample_ids = open(config["sample_ids"]).read().strip().split("\n")

def get_output(wildcards):
    assert config["layout"] in ["paired","single"]
    if config["layout"] == "paired":
        trimmed = expand('{outdir}/trimmed-pe/{sample_id}_{mate}.fastq.gz',sample_id=sample_ids,outdir=outdir,mate=["1","2"])
        qc_steps = expand('{outdir}/{step}-pe/{sample_id}_{mate}_fastqc.html',sample_id=sample_ids,outdir=outdir,mate=["1","2"],step=["qc-raw","qc-trimmed"])
        removed =  expand('{outdir}/cleaned-pe/{sample_id}_{mate}.fastq.gz',sample_id=sample_ids,outdir=outdir,mate=["1","2"])
    else:
        trimmed = expand('{outdir}/trimmed-se/{sample_id}.fastq.gz',sample_id=sample_ids,outdir=outdir)
        qc_steps = expand('{outdir}/{step}-se/{sample_id}_fastqc.html',sample_id=sample_ids,outdir=outdir,step=["qc-raw","qc-trimmed"])
        removed = expand('{outdir}/cleaned-se/{sample_id}.fastq.gz',sample_id=sample_ids,outdir=outdir)
    requested = removed  + qc_steps #trimmed + qc_steps + removed
    if config["profiling"]:
        abundance = expand("{outdir}/metaphlan/abundance/{sample_id}.txt",sample_id=sample_ids,outdir=outdir)
        requested += abundance
    if config["infer_strandness"]:
        strandness =  expand("{outdir}/metaphlan/strandness/{sample_id}.txt",sample_id=sample_ids,outdir=outdir)
        requested += strandness
    return requested
        

rule all:
    input:
        get_output

### Quality control before reads trimming
rule qc_raw_pe:
    input:
        fastq = indir + '/{sample_id}_{mate}.fastq.gz',
    output:
        report = '{outdir}/qc-raw-pe/{sample_id}_{mate}_fastqc.html',
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        fastqc -o {wildcards.outdir}/qc-raw-pe {input.fastq} 
        """


### Quality control before reads trimming
rule qc_raw_se:
    input:
        fastq = indir + '/{sample_id}.fastq.gz',
    output:
        report = '{outdir}/qc-raw-se/{sample_id}_fastqc.html',
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        fastqc -o {wildcards.outdir}/qc-raw-se {input.fastq} 
        """

### Reads trimming
rule trimming_pe:
    input:
        fastq_1 = indir + '/{sample_id}_1.fastq.gz',
        fastq_2 = indir + '/{sample_id}_2.fastq.gz'
    output:
        fastq_1 = '{outdir}/trimmed-pe/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/trimmed-pe/{sample_id}_2.fastq.gz',
        report_1 = '{outdir}/log/{sample_id}/trimming_statistics_1.txt',
        report_2 = '{outdir}/log/{sample_id}/trimming_statistics_2.txt',
    threads: 2
    log: '{outdir}/log/{sample_id}/trimming.txt'
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        export PATH=$PATH:/apps/home/lulab_jinyunfan/qhsky1/tools/TrimGalore-0.6.6
        trim_galore --phred33 --paired \
        --cores {threads}  -o {wildcards.outdir}/trimmed-pe \
         --basename {wildcards.sample_id} {input.fastq_1} {input.fastq_2} > {log} 2>&1
        mv {wildcards.outdir}/trimmed-pe/{wildcards.sample_id}_val_1.fq.gz {output.fastq_1}
        mv {wildcards.outdir}/trimmed-pe/{wildcards.sample_id}_val_2.fq.gz {output.fastq_2}
        mv {wildcards.outdir}/trimmed-pe/{wildcards.sample_id}_1.fastq.gz_trimming_report.txt {output.report_1}
        mv {wildcards.outdir}/trimmed-pe/{wildcards.sample_id}_2.fastq.gz_trimming_report.txt {output.report_2}
        """


rule trimming_se:
    input:
        fastq = indir + '/{sample_id}.fastq.gz'
    output:
        fastq = "{outdir}/trimmed-se/{sample_id}.fastq.gz",
        report = "{outdir}/log/{sample_id}/trimming_statistics.txt"
    threads: 2
    log: '{outdir}/log/{sample_id}/trimming.txt'
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        export PATH=$PATH:/apps/home/lulab_jinyunfan/qhsky1/tools/TrimGalore-0.6.6
        trim_galore --phred33 --cores {threads} -o {wildcards.outdir}/trimmed-se --basename {wildcards.sample_id} {input.fastq} > {log} 2>&1
        mv {wildcards.outdir}/trimmed-se/{wildcards.sample_id}_trimmed.fq.gz {output.fastq}
        mv {wildcards.outdir}/trimmed-se/{wildcards.sample_id}.fastq.gz_trimming_report.txt {output.report}
        """


### Quality control for cleaned reads
rule qc_trimmed_pe:
    input:
        fastq = '{outdir}/trimmed-pe/{sample_id}_{mate}.fastq.gz',
    output:
        report = '{outdir}/qc-trimmed-pe/{sample_id}_{mate}_fastqc.html'
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        fastqc -o {wildcards.outdir}/qc-trimmed-pe {input.fastq} 
        """

rule qc_trimmed_se:
    input:
        fastq = '{outdir}/trimmed-se/{sample_id}.fastq.gz',
    output:
        report = '{outdir}/qc-trimmed-se/{sample_id}_fastqc.html'
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        fastqc -o {wildcards.outdir}/qc-trimmed-se {input.fastq} 
        """

### remove human reads

rule remove_unwanted_sequences_pe:
    input:
        fastq_1 = '{outdir}/trimmed-pe/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/trimmed-pe/{sample_id}_2.fastq.gz',
    output:
        fastq_1 = '{outdir}/cleaned-pe/{sample_id}_1.fastq.gz',
        fastq_2 = '{outdir}/cleaned-pe/{sample_id}_2.fastq.gz',
        bam = '{outdir}/bam/unwanted-sequences/{sample_id}.bam'
    threads: 4
    params:
        prefix = config["unwanted_sequences"] if "unwanted_sequences" in config else ""
    log: '{outdir}/log/remove-unwanted-reads/{sample_id}.log'
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        prefix={params.prefix}
        if [ ${{#prefix}} -gt 0 ];then
          bowtie2 --no-unal -p {threads} -1 {input.fastq_1} -2 {input.fastq_2}  --un-conc-gz {wildcards.outdir}/cleaned-pe/{wildcards.sample_id}_%.fastq.gz --no-discordant --end-to-end -x {params.prefix} 2> {log} | samtools view -b > {output.bam}
        else
          ln -s $PWD/{input.fastq_1} $PWD/{output.fastq_1}
          ln -s $PWD/{input.fastq_2} $PWD/{output.fastq_2}
          touch {output.bam}
        fi
        """

rule remove_unwanted_sequences_se:
    input:
        fastq = '{outdir}/trimmed-se/{sample_id}.fastq.gz',
    output:
        fastq = '{outdir}/cleaned-se/{sample_id}.fastq.gz',
        bam = '{outdir}/bam/unwanted-sequences/{sample_id}.bam'
    threads: 4
    params:
        prefix = config["unwanted_sequences"] if "unwanted_sequences" in config else ""
    log: '{outdir}/log/{sample_id}/remove-human-reads.log'
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        prefix={params.prefix}
        if [ ${{#prefix}} -gt 0 ];then
          bowtie2 --no-unal -p {threads} -U {input.fastq} --un-gz {wildcards.outdir}/cleaned-se/{wildcards.sample_id}.fastq.gz --end-to-end -x $prefix 2> {log} | samtools view -b > {output.bam}
        else
          ln -s $PWD/{input.fastq} $PWD/{output.fastq}
          touch {output.bam}
        fi
        """


rule metaphlan:
    input:
        fastq = '{outdir}/cleaned-pe/{sample_id}_1.fastq.gz' if config["layout"] == "paired" else '{outdir}/cleaned-se/{sample_id}.fastq.gz',
    output:
        abundance = "{outdir}/metaphlan/abundance/{sample_id}.txt",
        bt2output = "{outdir}/metaphlan/hits/{sample_id}.txt.gz",
        bam = "{outdir}/metaphlan/bam/{sample_id}.bam"
    threads: 4
    params: 
        bt2db = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases"
    log: '{outdir}/log/metaphlan/{sample_id}.txt' 
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/biobakery3/bin
        metaphlan {input.fastq} --min_ab 0.001  --input_type fastq  --bowtie2db {params.bt2db} --bowtie2out {wildcards.outdir}/metaphlan/hits/{wildcards.sample_id}.txt --output_file {output.abundance} --nproc {threads} --tmp_dir tmp --samout {wildcards.outdir}/metaphlan/bam/{wildcards.sample_id}.sam
        gzip {wildcards.outdir}/metaphlan/hits/{wildcards.sample_id}.txt
        samtools view -Sb {wildcards.outdir}/metaphlan/bam/{wildcards.sample_id}.sam  > {output.bam}
        rm {wildcards.outdir}/metaphlan/bam/{wildcards.sample_id}.sam
        """

rule get_detectable_genera:
    input:
        abundance = "{outdir}/metaphlan/abundance/{sample_id}.txt"
    output:
        genera_list = "{outdir}/metaphlan/genera/{sample_id}.txt"
    log: '{outdir}/log/get-detectable-genus/{sample_id}.txt'
    shell:
        """
        export PATH=~/qhsky1/miniconda/envs/biobakery3/bin:$PATH
        scripts/get-detectable-genus.py -i {input.abundance} -o  {output.genera_list} > {log} 2>&1
        """

rule get_genera_recurrence:
    input:
        genera_lists = expand("{outdir}/metaphlan/genera/{sample_id}.txt",outdir=outdir,sample_id=sample_ids)
    output:
        recurrence = "{outdir}/metaphlan/genera.recurrence.txt"
    shell:
        """
        cat {input.genera_lists} | sort | uniq -c | awk 'BEGIN{{FS=" ";OFS="\t";}}{{print $2,$1}}' | sort -k 2 -nr  > {output.recurrence}
        """    


rule infer_strandness:
    input:
        bam = "{outdir}/metaphlan/bam/{sample_id}.bam"
    output:
        fraction = "{outdir}/metaphlan/strandness/{sample_id}.txt"
    log: '{outdir}/log/infer-strandness/{sample_id}.txt'
    shell:
        """
        export PATH=~/qhsky1/miniconda/envs/bioinfo-env/bin:$PATH
        scripts/infer-strandness.py -i {input.bam} -o {output.fraction} > {log} 2>&1
        """
