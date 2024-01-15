---
layout: page
title: Bioinformatic pipelines for genomic sequencing data
permalink: /tutorials/bioinfo-week2/
---

### What you will need:

- A web browser <br>

### Getting started


Before diving deeper into phylogenetic analysis, this week we will take a step back and look at bioinformatic pipelines for working with raw sequencing data. The [sandbox.bio][sandbox] tutorials we'll run through are short and can be run directly in your browser, so no need to install anything.

[sandbox]: <sandbox.bio>


### Command line basics

If you are not familiar with the basics of working from the command line, I highly recommend the [sandbox.bio command-line tutorial][terminal-basics]. The tutorial will walk you through the basics of working with files and folders from the command line.

[terminal-basics]: <https://sandbox.bio/tutorials?id=terminal-basics>

### Aligning reads to a reference

Most modern sequencing platforms like Illumina and PacBio produce a large number of sequencing reads covering different regions of a target genome. At first we won't necessarily know which region of the genome each read corresponds to, so the first step of most sequencing pipelines is to "map" or align reads to a previously sequenced reference genome. The [Sequence alignment with bowtie2 tutorial][alignment] will show you how can map raw sequencing reads in a *fastq* file to a reference genome for the lambda phage.

[alignment]: <https://sandbox.bio/tutorials/bowtie2-intro>

***Brain stretcher:*** For most well-studied pathogens we generally have access to previously assembled reference genomes. But what would would we do if we were sequencing a pathogen *de novo* for the first time?

### Working with SAM/BAM files

After mapping, aligned sequence reads are generally stored as SAM (Sequence Alignment/Map) or BAM (Binary Alignment Map) format files. The [BAM parsing with samtools tutorial][samtools] will show you how to work with these files including sorting and indexing the aligned reads for downstream use.

[samtools]: <https://sandbox.bio/tutorials?id=samtools-intro>

***Hint:*** While not mentioned in the tutorial, you can also use samtools to inspect the depth of coverage at each genomic position i.e. the number of sequence reads that map to each position:

```
samtools depth <sorted_bam_file.bam> > <log_file.tsv>
```

***Recommendation:*** The sandbox tutorial uses IGV to visualize the mapped reads directly in your browser. I would also recommend [Tablet][tablet] for visualizing larger assemblies and alignments on your own machine.

[tablet]: <https://ics.hutton.ac.uk/tablet/>


### Variant calling

After aligning our sequencing reads, we typically want to know how our reads differ one another and the reference sequence. For instance, we may want to discover variants like single nucleotide polymorphisms (SNPs) for building phylogenetic trees. The [Variant Calling tutorial][variants] will walk you through how you can call variants and distinguish real variants from sequencing errors using bcftools. Then you will get to combine what you’ve learned to decode a secret message in the sequenced DNA!!

[variants]: <https://sandbox.bio/tutorials/dna-secrets>

### More tutorials to check out

The [DNA sequencing QC tutorial][quality] demonstrates how to evaluate the quality of sequencing reads and filter out low quality data.

[quality]: <https://sandbox.bio/tutorials/fastp-intro>

The [Viral Amplicon Sequencing tutorial][amplicon] demonstrates how to add quality control and primer trimming steps to a pipeline for determining the consensus sequence of a given pathogen strain.

[amplicon]: <https://sandbox.bio/tutorials/viral-amplicon>

If you really want to get advanced, you can check out the [Genomic intervals with bedtools tutorial][bedtools] to see how to locate and annotate specific genomic features like genes and primer binding sites in genomes

[bedtools]: <https://sandbox.bio/tutorials/bedtools-intro>
