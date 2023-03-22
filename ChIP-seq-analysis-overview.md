ChIP-seq analysis overview
================
Manan Shah
16/03/2023

## Experimental Design

The design of your experiment largely depends on your end goal.

If your goal is just to identify the binding sites of your TF or
location of hisone modifications, while not recommended, it is possible
to just do 1 biological replicate.

However, if you would like to compare and do differential binding
analysis (i.e., compare the binding of a TF between 2 conditions or
treatments etc), you need **at least 2 biological replicates**.

There are a range of controls you can use for ChIP-seq

  - Input (AKA whole cell extract/WCE) - This is essentially just
    sonicated genomic DNA. This helps account for any sonication, PCR or
    sequencing biases. It is the most commonly used and accepted control
    for ChIP-seq. Ideally, having a matched input for each IP you are
    doing is recommended. However, if using the same cell line and not
    doing anything that is likely to significantly change the chromatin
    landscape, you can get away with doing less inputs.
  - IgG - Used more commonly for ChIP-qPCR. Typically not recommended
    for ChIP-seq, as an IgG IP should have little to no DNA which means
    when sequenced, there will be a lot of PCR duplicates making it seem
    like regions are enriched. There are also issues over what species
    to use IgG from.
  - For a TF, a KO of the TF in the same cell line is perhaps the best
    control. However, unless this is already part of your experiment,
    this can be very time consuming, difficult or even impossible for
    some TFs. It also has similar drawbacks to using IgG as a control.
    [Here](https://academic.oup.com/nar/article/47/11/5735/5494779) is a
    paper that argues for this control.

### Sequencing parameters

Sequencing length - As in general, longer the better, however, for
ChIP-seq reads probably don’t need to be too long, anything above 50bp
should be fine.

Paired-end vs single-end - In general, we typically do single-end
sequencing for ChIP-seq. Paired-end data for ChIP-seq doesn’t provide
that much extra information and is general not worth the cost. However,
it does have its benefits such as easier to detect PCR duplicates,
ability to get exact fragment sizes. The one scenario where you should
consider paired-end sequencing is when you are most interested in
similar or repetitive regions.

Sequencing depth - depends on your TF. But typically for TFs, aim for
**\>30 million reads**. For histone modifications, aim for **\>50
million reads**.

## Quality Control

The first aspect to look at for ChIP-seq quality control is the raw
reads themselves. Like RNA-seq and most high throughput sequencing data,
the raw reads will be fastq files. And FastQC is a great tool to use to
assess the quality of these reads.

*Installing FastQC*

``` bash
conda install -c bioconda fastqc
```

Running fastqc on all files with the extension .fastq.gz in the current
folder:

``` bash
fastqc *.fastq.gz
```

After running FastQC, it produces a html file for each fastq file with a
report on various metrics. The main ones to keep an eye out for in
ChIP-seq analysis are:

  - Base Sequence Quality - there a quality score (phred score) for each
    base in the read which is based on the chance that base is called
    wrong. Typically with modern sequencing machines and protocols, you
    generally don’t have too much of an issue with this (and Ramaciotti
    and other sequencing services generally have a minimum threshold for
    this as well). However, it is still good to check, you will likely
    notice a dip in the quality towards the start and end of the reads
    (particularly the ends). If the quality is dropping off massively
    towards the end, that may be an issue, however trimming should take
    care of that and in most cases the reads will be salvageable.

*Installing trim-galore*

``` bash
conda install -c bioconda trim-galore
```

Unlike RNA-seq, read quality and adapter trimming is done regularly for
ChIP-seq data. Trim-galore\! is a wrapper around the read trimmer
cutadapt and is useful as it has a lot of quality of life features such
as automatic adapter sequence detection, and fastqc automatically being
run on the trimmed files.

Trimming a read file to remove adapter sequences, bases below a Phred
quality of 20, and reads that are less than 35bp after trimming using
trim-galore:

``` bash
trim-galore {fastq.gz} -q 20 --fastqc_args {output directory for fastqc file} --length 35 -o {output directory for trimmed fastq file} --stringency 1 --cores {cores} --basename {samplename}
```

trim-galore will also automatically run fastqc on the trimmed files for
you.

## Aligning ChIP-seq reads

After trimming, it is time to align your reads to the genome. This
generally the most computationally intensive component of ChIP-seq
analysis so if running on your own computer, it may take many hours per
file.

There are two aligners mostly used for ChIP-seq (short read aligners):

  - BWA
  - Bowtie2

Both should produce very similar results. I generally use Bowtie2 with
the –very-sensitive parameter. For typical ChIP-seqs, I will generally
use the default values for all the other parameter. However, if you are
interested in repetitive or similar regions (such as HBG1/2), there are
some parameters you should look into here.

Installing Bowtie2

``` bash
conda install -c bioconda bowtie2
```

You will need an index for every genome you want to align reads to.
Luckily, the bowtie2 website provides pre-built indexes for the most
commonly used genomes.

For single-end reads:

``` bash
Bowtie2 –p {cores} –x {index directory} –U read1.fq,read2.fq,read3.fq,read4.fq --very-sensitive –S {sample name}.sam
```

For paired-end reads:

``` bash
Bowtie2 –p {cores} –x {index directory} –1 read1_1.fq,read2_1.fq… -2 read1_2.fq,read2_2.fq… --very-sensitive –S {sample name}.sam
```

This produces a sequence alignment map (SAM) file. These files are
tpyically very large so one of the first things we will want to do is to
compress them into binary alignment map (BAM) files. We will use
samtools to do this.

Installing samtools

``` bash
conda install -c bioconda samtools
```

Along with compressing the file into a BAM file, we will also filter the
file for any poorly mapped reads (below a mapping quality (mapq) of \<
10), sort the reads by coordinate (can also be sorted by other
attributes), and index the file (which is needed for visualising bam
files in IGV and for other purposes).

``` bash
samtools view –u -@ {cores} –q 10 {alignedfile.sam} | samtools sort –O bam –o {sortedfile.bam} -@ {cores} -
```
