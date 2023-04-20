# Bacteria noncoding transcripts analysis

- This pipeline takes paired metagenomic sequencing (MGX) and metatranscriptomic sequencing (MTX) data as input, and annotate intergenic regions with RNA expression
- Support paired end / single end data

## Preprocess MGX reads

- trim adaptor, run fastqc
- remove unwanted sequence (optional)
- run metaphlan for taxonomy profiling

```bash
snakemake --snakefile snakefiles/preprocessing.snakefile --configfile config/preprocessing/test-pe.mgx.yaml
snakemake --snakefile snakefiles/preprocessing.snakefile --configfile config/preprocessing/test-se.mgx.yaml
```

## Preprocess MTX reads

- trim adaptor, run fastqc
- remove unwanted sequence (optional)
- run metaphlan for taxonomy profiling
- infer strandness of the RNA library based on marker gene mapping result

```bash
snakemake --snakefile snakefiles/preprocessing.snakefile --configfile config/preprocessing/test-pe.mtx.yaml
snakemake --snakefile snakefiles/preprocessing.snakefile --configfile config/preprocessing/test-se.mtx.yaml
```

## Assemble MGX reads and annotate contigs

- assemble MGX reads with megahit
- run prodigal for gene prediction
- annotate predicted gene with Pfam & hmmsearch
- search contig for known noncoding RNA with Rfam & cmsearch (optional)

```bash
# paired end data
snakemake --snakefile snakefiles/mgx-analysis.snakefile --configfile config/mgx-analysis/test-pe.yaml
# single end data
snakemake --snakefile snakefiles/mgx-analysis.snakefile --configfile config/mgx-analysis/test-se.yaml
```

## Mapping MTX reads

- Map MTX reads to paired contigs
- Assemble transcripts with stringtie
- Annotate transcripts in a gene centric manner

```bash
# paired end data
snakemake --snakefile snakefiles/mtx-analysis.with.mgx.snakefile --configfile config/mtx-analysis-with-mgx/test-pe.yaml
# single end data
snakemake --snakefile snakefiles/mtx-analysis.with.mgx.snakefile --configfile config/mtx-analysis-with-mgx/test-se.yaml
```




