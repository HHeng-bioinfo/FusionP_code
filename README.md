# Plasmid Fusion Analysis Pipeline

## Overview

This R script analyzes plasmid fusion events. It generates circular genome visualizations showing gene annotations and BLAST mapping results to identify potential plasmid fusion events.


## Requirements

### R Version

-   R \>= 4.3.0

### Required R Packages

``` r
circlize
GenomicRanges
data.table
ggplot2
grid
IRanges
paletteer
```

### Installation

``` r
install.packages(c("circlize", "data.table", "ggplot2", "grid", "paletteer"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "IRanges"))
```

## Input Files

### Required Input Files
0.  **Genome file**:
-   `plasmid_query.fa` - Each plasmid recovered and contigs named with filename as the beginning as "filename_number"
-   `plasmid_collection.fa` - The collection of all query plasmids

1.  **Plasmid annotation files**:

-   `OMAP_pls.tsv` - Processed plasmid data by OMAP-KP
-   `OMAP_pls_raw.tsv` - Raw plasmid sequences by OMAP-KP

2.  **BLAST results**:

-   `blast.out` - BLAST alignment results 

``` shell
makeblastdb -in plasmid_collection.fa -dbtype nucl  -parse_seqids -out plasmid_collection_db

blastn -query plasmid_query.fa  -db plasmid_collection_db  -evalue 1e-10  -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp slen qlen' -dust no -soft_masking false -max_target_seqs 100000000 | awk -v filename="$filename" '{print filename"\t"$0}' > blast.out

``` 

3.  **Annotation files**:

-   `mob_type.tsv` - MOB-typing results with replicon types and mobility predictions
-   `bakta_anno.tsv` - Bakta gene annotations tsv file
``` shell
mob_typer --infile  plasmid_query.fa  --out_file mob_type.tsv
bakta --output ./ --prefix plasmid_query   --meta  --keep-contig-headers  --locus-tag plasmid_query    --skip-rrna --skip-trna --skip-plot --force  plasmid_query.fa

``` 

4.  **Metadata files**:
-   `bact_meta.tsv` - data.table with strain, ST, KL, Year, and Area
-   `meta` - PDG000000012.1904.metadata.tsv (from NCBI Pathogen Watch)
-   `plasmids.repr.clu` - Reference plasmid cluster file (from Klety)
-   `CARD_AMR_clustered.csv` - Carbapenemase gene definitions (from Kleborate)

### Input File Formats

#### BLAST output format (tab-separated):

```         
filename qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen
```

#### MOB-typing output:

```         
sample_id rep_type(s) predicted_mobility size
```

#### Bakta annotation:

```         
#Sequence Id Type Start Stop Strand Locus Tag Gene Product DbXrefs
```

#### bact_meta tsv:

```         
strain ST KL Year Area
```


## Output Files

### Generated Files

1.  **blast_dt. tsv** - Filtered and processed BLAST alignments
2.  **blast_map.tsv** - Summary of BLAST mappings with coverage statistics
3.  **blast_map_f.tsv** - Filtered BLAST mappings (complete genomes, coverage \>= 80%)
4.  **evidence.txt** - List of strains with fusion evidence
5.  **blast_dt_list.txt** - Detailed BLAST data for plotted comparisons
6.  **hv\_[cluster]\_[strain].pdf** - Circular genome plots showing:

-   Reference plasmid gene annotations
-   BLAST alignment coverage
-   Gene type color coding (resistance, virulence, plasmid genes, IS elements)
-   Temporal information (color-coded by year)


## Notes

-   The script focuses on PC_499 cluster (pLVPK virulence plasmid) fusion events
-   Complete genomes and single-contig plasmids are also prioritized for reference
-   The script uses IRanges for accurate coverage calculation accounting for overlaps

