# Slide-Seq Pipeline

This is the new alignment pipeline Slide-Seq data.

## Goals

 - Demultiplex sequencing data from Illumina binary files
 - Match sequenced barcodes to a known set of slideseq puck barcodes to get spatial information
 - Align reads to a genome using STAR (maybe STARsolo)
 - Count features and summarize as a gene-cell matrix
 - Do all of this fast, efficiently, and reproducibly. When possible we'd like to disentangle the different steps so that we can perform them separately.

## Explicit Non-Goals

 - Won't work in arbitary environments. This pipeline is built to run at the Broad on the local computing infrastructure.
 - Doesn't allow much flexibility in the workflow. This makes the entire process simpler to configure and design.

## Unclear goals

 - Create genome references. Previous version used a slightly customized reference, can we switch to a standard STAR reference?

## Requirements

The pipeline is designed to work from an Anaconda3 environment on the UGER cluster. We provide an `environment.yml` file with the requirements.

```shell
conda env create -n [environment name] -f environment.yml
```

The following tools are also needed but are available on UGER as dotkits:

1) `Picard`: - `use Picard-Tools`
1) `Java` - `use Java-1.8`
1) `Samtools` - `use Samtools`

Finally, we use `Drop-seq tools` from [the McCarroll lab](https://github.com/broadinstitute/Drop-seq), which can be
 

## Build genome reference

The Slide-seq tools need a genome reference in specific format for alignment and analysis. You could build a genome reference based on input fasta and gtf files. 

Command: 
```
build_reference.py manifest_file
```

Check `example.buildreference.txt` for manifest file format

## Run the Slide-seq tools

Notice: 
1) Check `example.manifest.txt` for manifest file format
4) See `user_doc.txt` for detailed usage of the Slide-seq tools. 

