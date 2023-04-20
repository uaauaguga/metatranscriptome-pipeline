# Bacteria noncoding transcripts analysis

- This pipeline takes paired metagenomic sequencing (MGX) and metatranscriptomic sequencing (MTX) data as input, and annotate intergenic regions with RNA expression

- Preprocess MGX reads
```
snakemake --snakefile snakefiles/preprocessing.snakefile --configfile config/preprocessing/test-pe.mgx.yaml
snakemake --snakefile snakefiles/preprocessing.snakefile --configfile config/preprocessing/test-se.mgx.yaml
```

- Preprocess MTX reads

```
snakemake --snakefile snakefiles/preprocessing.snakefile --configfile config/preprocessing/test-pe.mtx.yaml
snakemake --snakefile snakefiles/preprocessing.snakefile --configfile config/preprocessing/test-se.mtx.yaml
```

- Assemble MGX reads and annotate contigs

```
snakemake --snakefile snakefiles/mgx-analysis.snakefile --configfile config/mgx-analysis/test-pe.yaml
snakemake --snakefile snakefiles/mgx-analysis.snakefile --configfile config/mgx-analysis/test-se.yaml
```

- Mapping MTX reads

```
snakemake --snakefile snakefiles/mtx-analysis.with.mgx.snakefile --configfile config/mtx-analysis-with-mgx/test-pe.yaml
snakemake --snakefile snakefiles/mtx-analysis.with.mgx.snakefile --configfile config/mtx-analysis-with-mgx/test-se.yaml
```




