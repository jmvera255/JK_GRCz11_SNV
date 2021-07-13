Short Nucleotide Variant Discovery From Zebrafish Whole-exome Sequencing
================
Jess Vera
2021-06-28

-   [Overview](#overview)
-   [Software](#software)
-   [Prepare Accessory Files](#prepare-accessory-files)
-   [Prepare BAM File For Variant
    Calling](#prepare-bam-file-for-variant-calling)
    -   [Mark PCR Duplicate Reads](#mark-pcr-duplicate-reads)
    -   [Modify BAM file read group header for easier downstream
        analysis](#modify-bam-file-read-group-header-for-easier-downstream-analysis)
    -   [Clean up and index BAM file](#clean-up-and-index-bam-file)
    -   [Base Quality Score
        Recalibration](#base-quality-score-recalibration)
-   [Call Variants](#call-variants)
    -   [Filter “Raw” Variant Calls](#filter-raw-variant-calls)
-   [Homozygosity Score](#homozygosity-score)
    -   [Define Genomic Blocks](#define-genomic-blocks)
    -   [Merge VCF Files](#merge-vcf-files)
    -   [Counts Variants Per Block](#counts-variants-per-block)
    -   [Calculate and Plot Homozygosity
        Score](#calculate-and-plot-homozygosity-score)
-   [Identifying Candidate Variants](#identifying-candidate-variants)

# Overview

The following steps were performed to generate [variant call
format](https://github.com/samtools/hts-specs/blob/master/VCFv4.1.pdf)
files (aka VCF) that detail putative recessive SNPs and indels in
ENU-mutagenized zebrafish by sequencing and comparing a pool of mutants
to a pool of their WT siblings. GATK has established “best practices”
pipelines for germline short variant discovery which was used to inform
the below steps, more details are available at
<https://gatk.broadinstitute.org/hc/en-us/articles/360035535932>.

# Software

All work was performed using the following open source software:

-   GATK 4.2.0.0
-   Samtools 1.12
-   Picardtools 2.25.4
-   bcftools 1.12

Along with accessory scripts in the `scripts/` directory of the
associated GitHub repository at
<https://github.com/jmvera255/JK_GRCz11_SNV.git>.

# Prepare Accessory Files

A zebrafish genome sequence fasta file and several accessory files
(`.dict` and `.fai`) will be needed to complete the below steps. These
accessory files can be made using `Picardtools` and `Samtools`. More
details are available from GATK’s website at
<https://gatk.broadinstitute.org/hc/en-us/articles/360035531652>.

# Prepare BAM File For Variant Calling

Sequencing reads are aligned to a reference genome. These alignments are
written to a SAM file and subsequently converted to a binary version
called a BAM file. Below are the steps outlined in GATK best practices
guide used to prepare the alignment data for variant calling, that goal
of which is to reduce false positive calls. More details are available
at [Data pre-processing for variant
discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912).

The below commands are setup to to run in a bash script in which the
variable `SAMPLE` has been predefined, for example
`SAMPLE='E10-4-2-WT'`.

## Mark PCR Duplicate Reads

    # Mark duplicates
    java -jar picard.jar MarkDuplicates \
    INPUT=$1 \
    METRICS_FILE=$SAMPLE.markdups.metrics.txt \
    OUTPUT=temp.$SAMPLE.markdups.bam

## Modify BAM file read group header for easier downstream analysis

    # modify read group of processed bam file
    java -jar picard.jar AddOrReplaceReadGroups \
    I=temp.$SAMPLE.markdups.bam \
    O=$SAMPLE.markdups.bam \
    RGLB=$SAMPLE RGSM=$SAMPLE RGPL=Illumina RGPU=$SAMPLE

## Clean up and index BAM file

    rm temp.$SAMPLE.markdups.bam
    samtools index $SAMPLE.markdups.bam

## Base Quality Score Recalibration

This is a two-step process that detects systematic errors made by the
sequencing machine when it estimated the accuracy of each base call. To
do so, this step requires a true positive set of possible variants used
to exclude regions around known polymorphisms from the analysis. I
downloaded `danio_rerio.vcf.gz` from
<ftp://ftp.ensembl.org/pub/release-102/variation/vcf/danio_rerio/danio_rerio.vcf.gz>
on 2021-05-14 which contains all germline variations from the current
Ensembl release for *Danio rerio*.

    # estimate base recalibration metrics
    gatk BaseRecalibrator \
    -I $SAMPLE.markdups.bam \
    -R Danio_rerio.GRCz11.dna.primary_assembly.fa -O $SAMPLE.bqsr.txt \
    --known-sites danio_rerio.vcf.gz

    # apply recalibration
    gatk ApplyBQSR \
    -R Danio_rerio.GRCz11.dna.primary_assembly.fa -I $SAMPLE.markdups.bam \
    --bqsr-recal-file $SAMPLE.bqsr.txt \
    -O $SAMPLE.bqsr.bam

# Call Variants

    # get germline variants
    gatk HaplotypeCaller \
    -R Danio_rerio.GRCz11.dna.primary_assembly.fa \
    -I $SAMPLE.bqsr.bam \
    -O $SAMPLE.unfiltered.vcf.gz \
    --native-pair-hmm-threads 2 \
    -bamout $SAMPLE.haplotypes.bam

## Filter “Raw” Variant Calls

The variants called by HaplotypeCaller will contain (likely many) false
positives. As I did not have a “true positive” variant data set for
these zebrafish samples I did not perform variant quality score
recalibration before hard filtering variants based on metrics like read
depth and quality score. The cutoffs used below were selected after
calculating various summary statistics as well as visualizing
distributions of read depth, quality, etc. for each vcf file. **Before
filtering one should do a thorough assessment of the “raw” vcf results
to determine appropriate filtering cutoffs as these will likely vary
from project to project and across sequencing runs.**

> One option for visualizing variant stats is to use `bcftools query` to
> collect variables in a table format and to plot the distributions via
> custom scripts. For example, the following will collect variant
> quality, root-mean-square quality, allelic read depth, and total read
> depth for each variant and plot the distributions using R:
>
>     bcftools query \
>     -f '%QUAL \t %INFO/MQ \t [%AD{1}] \t %INFO/DP \n' \
>     $SAMPLE.unfiltered.vcf.gz > $SAMPLE.query_stats.txt
>
>     Rscript scripts/plot-stats.R $SAMPLE $SAMPLE.query_stats.txt
>
> Another good options is `bcftools plot-vcfstats` which uses Python.

The below command will filter out SNPs within 10bp of an indel or other
other variant type and filter SNPs and indels that have &lt; 4 reads
supporting the variant and which have an RMS quality score &lt;= 25.
These were the filters applied to the fin regeneration mutant samples.

    bcftools filter -O z -g 10 \
    -i 'INFO/MQ >25 & FORMAT/AD[0:1] > 3' \
    $SAMPLE.unfiltered.vcf.gz \
    > $SAMPLE.filtered.vcf.gz

    # index the vcf file
    bcftools index -t $SAMPLE.filtered.vcf.gz

# Homozygosity Score

[Leshchiner et al. *Genome Res. 2012. 22:
1541-1548*](http://www.genome.org/cgi/doi/10.1101/gr.135541.111) defined
homozygosity score as the “ratio of heterozygous SNP calls between
control and mutant pools multiplied by the number of informative
homozygous SNP calls in the mutant pool”. While not clear, I
interpretted “informative homozygous SNP calls” to mean homozygous SNVs
unique in the mutant relative to the WT sibling.

## Define Genomic Blocks

To calculate the homozygosity score for distinct blocks of the GRCz11
genome, I first defined the genomic blocks using a custom Perl script:

    # get chromosome only fai
    grep -v -P "^K|^M" Danio_rerio.GRCz11.dna.primary_assembly.fa.fai \
    > GRCz11.chroms.fa.fai

    # get 10 kb blocks
    perl scripts/chrom_blocks.pl \
    -i GRCz11.chroms.fa.fai \
    -k 10 \
    > GRCz11.10kb.regions.txt

## Merge VCF Files

The filtered VCF files for a given pair of mutant and WT siblings are
merged to allow identification of homozygous SNVs unique to the mutants.

    # merge WT and Mut filtered vcf files
    bcftools merge -m both -O z -o merged.vcf.gz \
    WT.filtered.vcf.gz Mut.filtered.vcf.gz

    # index the merged vcf file
    bcftools index -t merged.vcf.gz

## Counts Variants Per Block

**Note** - the step used here for collecting SNV calls within distinct
genomic blocks is very inefficient. It is very likely that other, more
efficient options exists. I was however able to utilize a
high-throughput computing approach to break up the task and thus
drastically reduce how long this step required.

    bash scripts/count_blocks.sh merged.vcf.gz GRCz11.10kb.regions.txt \
    > 10kb.counts.txt

## Calculate and Plot Homozygosity Score

    Rscript scripts/calc_homozygous_score.R samples.txt

Where `samples.txt` is a list of different sample names corresponding to
different `counts.txt` files created in the previous step. **Note, this
script makes a number of assumptions about the different window sizes
profiled in the previous step and the organization of a number of
associated files. You will need to revise the script accordingly to work
in the context of your analysis.**

After this step you will have a list of the top score regions for each
WT/Mutant pair for each of the predefined genomic blocks.

# Identifying Candidate Variants

`bcftools view` can be used to select and filter variants of interest
within the top scoring genomic regions identified above or any other
regions of interest.

SnpEff was used to annotation variants within regions of interest
following SnpEff documentation, for example:

    java -jar snpEff.jar GRCz11.99 candidate_variants.vcf \
    > candidate_variants.annotated.vcf

More details are available at
<http://pcingola.github.io/SnpEff/se_introduction/>.
