
mtx2mgx = {}
with open(config["mtx_to_mgx"]) as f:
    for line in f:
        mtx_id, mgx_id = line.strip().split("\t")
        mtx2mgx[mtx_id] = mgx_id

sample_ids = open(config["sample_ids"]).read().strip().split("\n")
for sample_id in sample_ids:
    assert sample_id in mtx2mgx

outdir = config["outdir"]


def get_output(wildcards):
    bam = expand("{outdir}/bam/paired-contigs/{sample_id}.bam",sample_id=sample_ids,outdir=outdir)
    tx = expand('{outdir}/stringtie.30/paired-contigs/bed.annotated/{sample_id}.bed',sample_id=sample_ids,outdir=outdir)
    bw = expand('{outdir}/coverage/paired-contigs/{sample_id}.bigwig',sample_id=sample_ids,outdir=outdir)
    return bam + tx + bw #+ count

rule all:
    input:
        get_output

def get_contig_bt2idx(wildcards):
    bt2idx = []
    for i in ["1","2","3","4","rev.1","rev.2"]:
        mgx_id = mtx2mgx[wildcards.sample_id]
        path = f"{outdir}/assembly/{mgx_id}/final.contigs.{i}.bt2"
        bt2idx.append(path)
    return bt2idx


rule mapping_paired_contigs:
    input:
        fastq = ['{outdir}/cleaned-pe/{sample_id}_1.fastq.gz','{outdir}/cleaned-pe/{sample_id}_2.fastq.gz'] if config['layout'] == 'paired' else ['{outdir}/cleaned-se/{sample_id}.fastq.gz','{outdir}/cleaned-se/{sample_id}.fastq.gz'],
        bt2idx = get_contig_bt2idx
    output:
        bam = "{outdir}/bam/paired-contigs/{sample_id}.bam",
        unmapped = "{outdir}/unmapped-pe/paired-contigs/{sample_id}_1.fastq.gz" if config['layout'] == 'paired' else "{outdir}/unmapped-se/paired-contigs/{sample_id}.fastq.gz"
    log: "{outdir}/log/{sample_id}/mapping.paired.mgx.txt"
    params:
        layout = config["layout"]
    threads: 4
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        bt2idx_prefix=$(echo {input.bt2idx[0]} | sed 's/.1.bt2//g')
        if [ {params.layout} = "paired" ];then
          bowtie2 -p {threads} -1 {input.fastq[0]} -2 {input.fastq[1]} -x $bt2idx_prefix \
          --sensitive-local --no-discordant --no-unal \
          --un-conc-gz {outdir}/unmapped-pe/paired-contigs/{wildcards.sample_id}_%.fastq.gz 2> {log} \
          | samtools sort --output-fmt BAM --threads {threads} -o {output.bam}
        else
          bowtie2 -p {threads} -U {input.fastq[0]} -x $bt2idx_prefix \
          --sensitive-local --no-unal --un-gz {outdir}/unmapped-se/paired-contigs/{wildcards.sample_id}.fastq.gz 2> {log} \
          | samtools sort --output-fmt BAM --threads {threads} -o {output.bam}
        fi
        samtools index {output.bam}
        """

rule stringtie30:
    input:
        bam = '{outdir}/bam/paired-contigs/{sample_id}.bam'
    output:
        gff = '{outdir}/stringtie.30/paired-contigs/gff/{sample_id}.gff'
    threads: 2
    params:
        lib = {"no":"","forward":"--fr","reverse":"--rf"}[config["strandness"]]
    log: "{outdir}/log/{sample_id}/stringtie.30.txt"
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        stringtie -m 30 -s 2 {params.lib} -o {output.gff} -p {threads} {input.bam} > {log} 2>&1
        """    

rule format_stringtie30:
    input:
        gff = '{outdir}/stringtie.30/paired-contigs/gff/{sample_id}.gff'
    output:
        bed = '{outdir}/stringtie.30/paired-contigs/bed/{sample_id}.bed'
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        scripts/gff2bed.py --name transcript_id --value cov --gff {input.gff} -f transcript --bed {output.bed}
        """

def get_cds_intervals(wildcards):
    mgx_id = mtx2mgx[wildcards.sample_id]
    return f"{outdir}/annotation/{mgx_id}/metaprodigal.cds.pfam.bed"

def get_contig_index(wildcards):
    mgx_id = mtx2mgx[wildcards.sample_id]
    return  f"{outdir}/assembly/{mgx_id}/final.contigs.fa.fai"
    
rule annotate_transcripts:
    input:
        tx_bed = '{outdir}/stringtie.30/paired-contigs/bed/{sample_id}.bed',
        gene_bed = get_cds_intervals,
        faidx = get_contig_index
    output:
        bed = '{outdir}/stringtie.30/paired-contigs/bed.annotated/{sample_id}.bed'
    shell:
        """
        scripts/annotate-intervals.py --gene {input.gene_bed} --bed {input.tx_bed} --output {output.bed} --contig {input.faidx}
        """
    
rule get_coverage:
    input:
        bam = '{outdir}/bam/paired-contigs/{sample_id}.bam'        
    output:
        chromsize = '{outdir}/coverage/paired-contigs/{sample_id}.chrom.size',
        bedgraph = '{outdir}/coverage/paired-contigs/{sample_id}.bedgraph',
        bigwig = '{outdir}/coverage/paired-contigs/{sample_id}.bigwig'
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        samtools view -H {input.bam} | awk 'BEGIN{{OFS="\\t";FS="\\t";}}$1~"@SQ"{{seq_id=substr($2,4,length($2)-3);seq_length=substr($3,4,length($3)-3);print seq_id,seq_length}}' > {output.chromsize}
        bedtools genomecov -ibam {input.bam} -bg -split | LC_ALL=C sort -k1,1 -k2,2n > {output.bedgraph}
        bedGraphToBigWig {output.bedgraph} {output.chromsize} {output.bigwig} 
        """


