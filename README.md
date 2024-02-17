# lrsv_oysters
## Long Read Structural Variations for Oysters: A collection of executed scripts for Crassostrea gigas PacBio reads.


## Brief Background
A suite of scripts utilised for processing the Pacific oyster PacBio reads is available.


## Citation
Hyungtaek Jung et al. 2024: Unraveling the biotechnological potential of Crassostrea gigas: comparative genomics & structural variations, [BioRxiv](https://www.biorxiv.org/XXXX).


## Contents:
+ [LICENSE](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#license) 
+ [GETTING STARTED](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#getting-started)
+ [FAQ](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#faq)
+ [WIKI PAGE](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#wiki-page)
+ [AUTHORS](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#authors)
+ [COPYRIGHT](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#copyright)



## LICENSE
The current scripts are available under the MIT license and incorporate various open-source software. 


## GETTING STARTED
The present collection of scripts is primarily designed for analysing PacBio reads of the Pacific oyster.


### Tested Files
- C. gigas1: Genome assembly [GCA_011032805.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_011032805.1), Project number [PRJNA598006](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA598006) and Raw data [SRA_598006](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=598006) in NCBI
- C. gigas2: Genome assembly [GCF_902806645.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_902806645.1), Project number [PRJEB35351](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB35351) and Raw data [SRA_596972](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=596972) in NCBI


### General Usage
Each program has provided its own description.


### LRA and CuteSV
- For generating read mapping, statistics, and identifying structural variations in PacBio reads:
 	+ Requirement: LRA, Samtools, CuteSV
	+ Input: A fastq file.
	+ Output: A summary report of read mapping statistics and a VCF file detailing structural variations.

Usage
```
module load lra
module load samtools
module load cutesv

# LRA Indexing and Mapping for PacBio Data
lra index -CLR /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
lra align -CLR /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta /Cra_gig2/Cra_gig2_PBallRN.fastq -t 8 -p s > LRA_Cg1RCg2PB.sam

# Samtools: Conversion of SAM to BAM, Including Mapping Statistics
samtools faidx /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
samtools view -bT /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta LRA_Cg1RCg2PB.sam > LRA_Cg1RCg2PB.bam
samtools sort -@8 -T $TMPDIR/LRA_Cg1RCg2PB.sorted -O bam -o LRA_Cg1RCg2PB.sorted.bam LRA_Cg1RCg2PB.bam
samtools index LRA_Cg1RCg2PB.sorted.bam
samtools flagstat LRA_Cg1RCg2PB.sorted.bam &> LRA_Cg1RCg2PB_stats.txt

# CuteSV: Aligning Cg2 PacBio Data Against the Cg1 Reference Genome
# BAM File Generation via LRA and Samtools
# Utilise --genotype with CuteSV to Generate a Genotype Score, a Prerequisite for SURVIVOR Analysis
cuteSV --genotype LRA_Cg1RCg2PB.sorted.bam /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta Cg1RCg2PB_LRACuteSV.vcf /Cgig1/LR_MappingSV/LRA --threads 8 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5

```

- For specific parameter details, refer to the official websites of [LRA](https://github.com/ChaissonLab/LRA), [Samtools](https://github.com/samtools/samtools), and [CuteSV](https://github.com/tjiangHIT/cuteSV).


### Minimap2 and CuteSV
- For generating read mapping, statistics, and identifying structural variations in PacBio reads:
 	+ Requirement: Minimap2, Samtools, CuteSV
	+ Input: A fastq file.
	+ Output: A summary report of read mapping statistics and a VCF file detailing structural variations.

Usage
```
module load minimap2
module load samtools
module load cutesv

# Minimap2 Indexing and Mapping for PacBio Data
minimap2 -ax map-pb /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta /Cgig2/Cgig2_PBallRN.fastq > MM2_Cg1RCg2PB.sam

# Samtools: Conversion of SAM to BAM, Including Mapping Statistics
samtools faidx  /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
samtools view -bT /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta MM2_Cg1RCg2PB.sam > MM2_Cg1RCg2PB.bam
samtools sort -@8 -T $TMPDIR/MM2_Cg1RCg2PB.sorted -O bam -o MM2_Cg1RCg2PB.sorted.bam MM2_Cg1RCg2PB.bam
samtools index MM2_Cg1RCg2PB.sorted.bam
samtools flagstat MM2_Cg1RCg2PB.sorted.bam &> MM2_Cg1RCg2PB_stats.txt

# CuteSV: Aligning Cg2 PacBio Data Against the Cg1 Reference Genome
# BAM File Generation via Minimap2 and Samtools
# Utilise --genotype with CuteSV to Generate a Genotype Score, a Prerequisite for SURVIVOR Analysis
cuteSV --genotype MM2_Cg1RCg2PB.sorted.bam /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta Cg1RCg2PB_MM2CuteSV.vcf /Cgig1/LR_MappingSV/MM2 --threads 8 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5

```

- For specific parameter details, refer to the official websites of [Minimap2](https://github.com/lh3/minimap2), [Samtools](https://github.com/samtools/samtools), and [CuteSV](https://github.com/tjiangHIT/cuteSV).


### NGMLR and CuteSV
- For generating read mapping, statistics, and identifying structural variations in PacBio reads:
 	+ Requirement: NGMLR, Samtools, CuteSV
	+ Input: A fastq file.
	+ Output: A summary report of read mapping statistics and a VCF file detailing structural variations.

Usage
```
module load ngmlr
module load samtools
module load cutesv

# NGMLR Indexing and Mapping for PacBio Data
ngmlr -t 8 -r /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta -q /Cgig2/Cra_gig2_PBallRN.fastq -o NgmlR_Cg1RCg2PB.sam

# Samtools (specific version of 1.9 due to its comparability): Conversion of SAM to BAM, Including Mapping Statistics
samtools faidx /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
samtools view -bT /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta NgmlR_Cg1RCg2PB.sam > NgmlR_Cg1RCg2PB.bam
samtools sort -@8 -T $TMPDIR/NgmlR_Cg1RCg2PB.sorted -O bam -o NgmlR_Cg1RCg2PB.sorted.bam NgmlR_Cg1RCg2PB.bam
samtools index NgmlR_Cg1RCg2PB.sorted.bam
samtools flagstat NgmlR_Cg1RCg2PB.sorted.bam &> NgmlR_Cg1RCg2PB_stats.txt

# CuteSV: Aligning Cg2 PacBio Data Against the Cg1 Reference Genome
# BAM File Generation via NGMLR and Samtools
# Utilise --genotype with CuteSV to Generate a Genotype Score, a Prerequisite for SURVIVOR Analysis
cuteSV --genotype NgmlR_Cg1RCg2PB.sorted.bam /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta Cg1RCg2PB_NgMcuteSV.vcf /Cgig1/LR_MappingSV/NGMLR --threads 8 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5

```

- For specific parameter details, refer to the official websites of [NGMLR](https://github.com/philres/ngmlr), [Samtools](https://github.com/samtools/samtools), and [CuteSV](https://github.com/tjiangHIT/cuteSV).


### Vulcan and CuteSV
- For generating read mapping, statistics, and identifying structural variations in PacBio reads:
 	+ Requirement: Vulcan, Samtools, CuteSV
	+ Input: A fastq file.
	+ Output: A summary report of read mapping statistics and a VCF file detailing structural variations.

Usage
```
module load vulcan
module load samtools
module load cutesv

# Vulcan Indexing and Mapping for PacBio Data
vulcan -r /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta -i /Cgig2/Cra_gig2_PBallRN.fastq -w /Cgig1/LR_MappingSV/Vulcan/ -o VLCout_Cg1RCg2PB -t 8 -clr

# Samtools (specific version of 1.9 due to its comparability): Conversion of SAM to BAM, Including Mapping Statistics
samtools faidx /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
samtools view -bT /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta VLCout_Cg1RCg2PB.sam > VLCout_Cg1RCg2PB.bam
samtools sort -@8 -T $TMPDIR/VLCout_Cg1RCg2PB.sorted -O bam -o VLCout_Cg1RCg2PB.sorted.bam VLCout_Cg1RCg2PB.bam
samtools index VLCout_Cg1RCg2PB.sorted.bam
samtools flagstat VLCout_Cg1RCg2PB.sorted.bam &> VLCout_Cg1RCg2PB_stats.txt

# CuteSV: Aligning Cg2 PacBio Data Against the Cg1 Reference Genome
# BAM File Generation via Vulcan and Samtools v1.9
# Utilise --genotype with CuteSV to Generate a Genotype Score, a Prerequisite for SURVIVOR Analysis
cuteSV --genotype VLCout_Cg1RCg2PB.sorted.bam /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta Cg1RCg2PB_VLCcuteSV.vcf /Cgig1/LR_MappingSV/Vulcan --threads 8 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5

```

- For specific parameter details, refer to the official websites of [Vulcan](https://gitlab.com/treangenlab/vulcan), [Samtools](https://github.com/samtools/samtools), and [CuteSV](https://github.com/tjiangHIT/cuteSV).


### Winnowmap2 and CuteSV
- For generating read mapping, statistics, and identifying structural variations in PacBio reads:
 	+ Requirement: Winnowmap2, Samtools, CuteSV
	+ Input: A fastq file.
	+ Output: A summary report of read mapping statistics and a VCF file detailing structural variations.

Usage
```
module load winnowmap2
module load samtools
module load cutesv

# Winnowmap2 Indexing and Mapping for PacBio Data
vulcan -r /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta -i /Cgig2/Cra_gig2_PBallRN.fastq -w /Cgig1/LR_MappingSV/Vulcan/ -o VLCout_Cg1RCg2PB -t 8 -clr

# Step 1
meryl count k=15 output merylDB /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
meryl print greater-than distinct=0.9998 merylDB > Cg1Ref_repetitive_k15.txt

# Step 2
winnowmap -W Cg1Ref_repetitive_k15.txt -ax map-pb /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta /Cgig2/Cra_gig2_PBallRN.fastq > WnM_Cg1RCg2PB.sam

# Samtools (specific version of 1.9 due to its comparability): Conversion of SAM to BAM, Including Mapping Statistics
samtools faidx /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
samtools view -bT /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta WnM_Cg1RCg2PB.sam > WnM_Cg1RCg2PB.bam
samtools sort -@8 -T $TMPDIR/WnM_Cg1RCg2PB.sorted -O bam -o WnM_Cg1RCg2PBB.sorted.bam WnM_Cg1RCg2PBB.bam
samtools index WnM_Cg1RCg2PB.sorted.bam
samtools flagstat WnM_Cg1RCg2PB.sorted.bam &> WnM_Cg1RCg2PB_stats.txt

# CuteSV: Aligning Cg2 PacBio Data Against the Cg1 Reference Genome
# BAM File Generation via Winnowmap2 and Samtools
# Utilise --genotype with CuteSV to Generate a Genotype Score, a Prerequisite for SURVIVOR Analysis
cuteSV --genotype WnM_Cg1RCg2PB.sorted.bam /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta Cg1RCg2PB_WnMcuteSV.vcf /Cgig1/LR_MappingSV/WinnowMap --threads 8 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5

```

- For specific parameter details, refer to the official websites of [Winnowmap2](https://github.com/marbl/Winnowmap), [Samtools](https://github.com/samtools/samtools), and [CuteSV](https://github.com/tjiangHIT/cuteSV).


### SURVIVOR
- For simulating/evaluating SVs, merging, summarising, and comparing SVs within and among samples (here different aligners):
 	+ Requirement: SURVIVOR
	+ Input: A vcf file.
	+ Output: A summary report detailing the merging and summarization of statistics.

Usage
```
module load survivor

# Merging VCF files
SURVIVOR merge LRASV_files 1000 2 1 1 0 30 LRASV_Cg1Rmerged_Def.vcf
SURVIVOR merge MM2SV_files 1000 2 1 1 0 30 MM2SV_Cg1Rmerged_Def.vcf
SURVIVOR merge NgMSV_files 1000 2 1 1 0 30 NgMSV_Cg1Rmerged_Def.vcf
SURVIVOR merge VLCSV_files 1000 2 1 1 0 30 VLCSV_Cg1Rmerged_Def.vcf
SURVIVOR merge WnMSV_files 1000 2 1 1 0 30 WnMSV_Cg1Rmerged_Def.vcf

# Checking for Overlapping VCF Files with a PERL Script
# Ensure to run the overlap.txt check each time to prevent overlapp.txt issues
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' sample_merged.vcf | sed -e 's/\(.\)/\1 /g' > sample_merged_overlap.txt 

# Plotting
SURVIVOR genComp CuteSV_Cg1Rmerged_Def.vcf CuteSV_Cg1Rmerged_Def.mat.txt

```

- For specific parameter details, refer to the official website of [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR).


## FAQ

We encourage users to use the [Issues](https://github.com/OZTaekOppa/lrsv_oysters/issues).


## WIKI PAGE

Please see the GitHub page.


## AUTHORS

**Hyungtaek Jung** and et al. 2024.


## COPYRIGHT

The full scripts are distributed under the MIT license. 
