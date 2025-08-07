This contains data and code for the manuscript entitled

**A phylogenomic and metagenomic meta-analysis of bacterial diversity in the phyllosphere lifts a veil on Hyphomicrobiales dark matter**

Preprint available on bioRxiv : TBA

———————————————————————————————————

**Genome_and_metagenome_analysis**

Script: *Script-Genome-Analysis-version-For-Github.R*

This part is for the definition of the core genome of phyllosphere-associated Hyphomicrobiales using reference genomes and metagenomes, and the inference of a genome-based phylogenetic tree for this group.

———————————————————————————————————

**Bioproject-16S-rRNA-formating**

Script: *Script-BioProjects-Data-Formating-For-Github.R*

This part is for read trimming and ASV definition (with dada2) for each Bioproject, raw ASV table abundance and rarefaction. 

———————————————————————————————————
**Bioproject-16S-rRNA-analysing**

Script: *Script-BioProjects-Data-Analysing-For-Github.R*

This part is for the taxonomic classification of  ASVs for each Bioproject, comparison or classification with SILVA 138.2 and reference genomes as a reference database and calculation of abundance for major bacteria taxon and especially Hyphomicrobiales for each BioProject and each category within BioProjects.

———————————————————————————————————
**RAxML-tree-formating**

Script: *Convert-Tree.R*

This part is to convert RAxML output tree files in readable format for FigTree. 

- 16S_rRNA : tree built from complete 16S rRNA sequences extracted from reference genomes and metagenome

- rpoB : tree built from complete rpoB sequences extracted from reference genomes and metagenomes

- Lichenibacterium_is_1174-901-12: tree built fro reference 16S rRNA sequence from Lichenihabitantaceae (NCBI) and sequences labelled as 1174-901-12 in SILVA 138.1 database

———————————————————————————————————
**Comparison-rpoB-16SrRNA-barcoding**
Script: *Comparison-rpoB-16SrRNA-For-Github.R*
This part is for the taxonomic classification of ASVs from rpoB barcoding (BioProject PRJNA729807 = Leducq et al 2022, mBio) and comparison with relative abundances of phyllosphere-associated Hyphomicrobiales estimated from 16S rRNA genérale barcoding with the same samples.
