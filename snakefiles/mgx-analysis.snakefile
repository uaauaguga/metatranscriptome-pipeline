outdir = config["outdir"]


def get_output(wildcards):
    sample_ids = open(config["sample_ids"]).read().strip().split("\n")
    assert config["layout"] in ["paired","single"]
    contigs = expand("{outdir}/assembly/{sample_id}/final.contigs.fa",outdir=outdir,sample_id=sample_ids)
    bt2idx = expand("{outdir}/assembly/{sample_id}/final.contigs.1.bt2",outdir=outdir,sample_id=sample_ids)
    requested = contigs + bt2idx
    if config["gene_prediction"]:
        cds = expand("{outdir}/annotation/{sample_id}/metaprodigal.cds.gff",sample_id=sample_ids,outdir=outdir)
        requested += cds
    if config["gene_annotation"]:
        cds_pfam =  expand("{outdir}/annotation/{sample_id}/metaprodigal.cds.pfam.bed",sample_id=sample_ids,outdir=outdir)
        requested += cds_pfam
    if "rna_scanning" in config:
        cm_hits = expand("{outdir}/annotation/{sample_id}/{model}.rna.gff",sample_id=sample_ids,outdir=outdir,model=config["rna_scanning"])
        requested += cm_hits
    if "mapping" in config:
        bam = expand("{outdir}/bam/paired-contigs/{sample_id}.bam",sample_id=sample_ids,outdir=outdir)
        requested += bam
    #requested += expand("{outdir}/annotation/{sample_id}/dust.interval",sample_id=sample_ids,outdir=outdir)
    return requested


rule all:
    input:
        get_output


rule megahit:
    input:
        fastq = ['{outdir}/cleaned-pe/{sample_id}_1.fastq.gz','{outdir}/cleaned-pe/{sample_id}_2.fastq.gz'] if config['layout'] == 'paired' else ['{outdir}/cleaned-se/{sample_id}.fastq.gz','{outdir}/cleaned-se/{sample_id}.fastq.gz']
    output:
        contig = "{outdir}/assembly/{sample_id}/final.contigs.fa",
    log: "{outdir}/log/megahit/{sample_id}.txt"
    params: 
        layout = config["layout"]
    threads: 64
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        if [ {params.layout} = "paired" ];then
          echo {input.fastq}
          megahit -1 {input.fastq[0]} -2 {input.fastq[1]} --num-cpu-threads {threads} --out-dir {wildcards.outdir}/assembly/{wildcards.sample_id}_ > {log} 2>&1
        else
          megahit -r {input.fastq[0]} --num-cpu-threads {threads} --out-dir {wildcards.outdir}/assembly/{wildcards.sample_id}_ > {log} 2>&1
        fi
        rm -r {wildcards.outdir}/assembly/{wildcards.sample_id} && mv {wildcards.outdir}/assembly/{wildcards.sample_id}_ {wildcards.outdir}/assembly/{wildcards.sample_id}
        samtools faidx {output.contig}
        """


rule build_bowtie_index:
    input:
        fasta = "{outdir}/assembly/{sample_id}/final.contigs.fa"
    output:
        bt2idx =  ["{outdir}/assembly/{sample_id}/final.contigs." + c + ".bt2" for c in ["1","2","3","4","rev.1","rev.2"]]
    params:
        prefix = "{outdir}/assembly/{sample_id}/final.contigs"
    threads: 1 
    log:
        "{outdir}/log/{sample_id}-bowtie2-index.txt"
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo/bin
        bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1
        """

def get_contig_bt2idx(wildcards):
    bt2idx = []
    for i in ["1","2","3","4","rev.1","rev.2"]:
        mgx_id = wildcards.sample_id
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
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
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



rule gene_prediction:
    input:
        fasta = "{outdir}/assembly/{sample_id}/final.contigs.fa"
    output:
        gff = "{outdir}/annotation/{sample_id}/metaprodigal.cds.gff",
        fna = "{outdir}/annotation/{sample_id}/metaprodigal.cds.fna",
        faa = "{outdir}/annotation/{sample_id}/metaprodigal.orf.faa"
    log:
        "{outdir}/log/{sample_id}/metaprodigal.txt"
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        prodigal -p meta -f gff -i {input.fasta} -o {output.gff} -a {output.faa} -d {output.fna} > {log} 2>&1
        """

'''
rule dustmasker:
    input:
        fasta = "{outdir}/assembly/{sample_id}/final.contigs.fa"
    output:
        iv = "{outdir}/annotation/{sample_id}/dust.interval"
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        dustmasker -in {input.fasta} -out {output.iv} -level 20 -infmt fasta -outfmt interval
        """
'''


rule protein_scanning:
    input:
        fasta = "{outdir}/annotation/{sample_id}/metaprodigal.orf.faa",
        hmm = "reference/hmm/Pfam-A.hmm"
    output:
        tbl = "{outdir}/annotation/{sample_id}/pfam.domain.tbl",
    threads: 4
    log:
        "{outdir}/log/{sample_id}/scan.protein.txt"
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        hmmsearch --cpu 8 -Z 1000 --noali --cut_ga --tblout {output.tbl} {input.hmm} {input.fasta} > {log} 2>&1
        """

rule domain_annotation:
    input:
        gff = "{outdir}/annotation/{sample_id}/metaprodigal.cds.gff",
        tbl = "{outdir}/annotation/{sample_id}/pfam.domain.tbl"
    output:
        gff = "{outdir}/annotation/{sample_id}/metaprodigal.cds.pfam.gff"
    log: 
        "{outdir}/log/{sample_id}/cds.annotation.txt"
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        scripts/annotate-pfam-domain.py -i {input.gff} -t {input.tbl}  -o {output.gff} > {log} 2>&1
        """

rule reformat_CDS:
    input:
        gff = "{outdir}/annotation/{sample_id}/metaprodigal.cds.pfam.gff"
    output:
        bed = "{outdir}/annotation/{sample_id}/metaprodigal.cds.pfam.bed"
    shell:
        """
        scripts/gff2bed.py --gff {input.gff} --bed {output.bed} --name pfam_id --feature CDS --value conf
        """



rule rna_scanning:
    input:
        fasta = "{outdir}/assembly/{sample_id}/final.contigs.fa",
        cm = "reference/cm/{model}.cm"
    output:
        tbl = "{outdir}/annotation/{sample_id}/{model}.rna.tbl",
    threads: 4
    params:
        cutga = lambda wildcards: "--cut_ga" if wildcards.model == "Rfam" else ""
    log: 
        "{outdir}/log/{sample_id}/cmsearch.{model}.rna.txt"
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        cmsearch --hmmonly --cpu 8 -Z 1000 --noali {params.cutga} --tblout {output.tbl} {input.cm} {input.fasta} > {log} 2>&1 
        """

rule tbl2gff:
    input:
        tbl = "{outdir}/annotation/{sample_id}/{model}.rna.tbl"
    output:
        gff = "{outdir}/annotation/{sample_id}/{model}.rna.gff"
    shell:
        """
        export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin
        scripts/tbl-to-gff.py -i {input.tbl} -o {wildcards.outdir}/annotation/{wildcards.sample_id}/{wildcards.model}.rna.tmp.gff
        cat {wildcards.outdir}/annotation/{wildcards.sample_id}/{wildcards.model}.rna.tmp.gff | sort -k1,1 -k2,4n -k3,5n > {output.gff}
        rm {wildcards.outdir}/annotation/{wildcards.sample_id}/{wildcards.model}.rna.tmp.gff
        """




