# Slide-seq Pipeline

This is the new alignment pipeline for Slide-seq data. For more information about Slide-seq method itself, see [Rodriques & Stickels et al.](http://doi.org/10.1126/science.aaw1219) and [Stickels et al.](https://doi.org/10.1038/s41587-020-0739-1).

## Quickstart

### Install

`conda install` uses a fair amount of RAM, so you should do this from an interactive session and not on the login node.

```shell
use UGER  # to make the ish command available
ish -l h_vmem=4G  # create an interactive session with extra memory
git clone https://github.com/MacoskoLab/slideseq-tools.git  # clone this repository
cd slideseq-tools
conda env create -f environment.yml  # creates an environment named `slideseq`
```

You can name the environment something else by adding `-n [env_name]`. There's no reason to do this, but I figured out how to make it work with an arbitary env name, so I wanted to mention it.

### Submitting a flowcell

Again, you should do this from an interactive session, as downloading the worksheet is a little too strenuous for the login node.

```shell
use UGER  # for submitting jobs to the cluster
use Google-Cloud-SDK  # for accessing the Google worksheet
conda activate slideseq
submit_flowcell FLOWCELL [FLOWCELL...]
```

This is will submit a set of jobs to process the flowcell(s). You will get emails when the jobs start and end.

 - One job per lane of the flowcell, for demuxing
 - For each lane, an array of alignment jobs for each library sequenced
 - An array of processing jobs, one for each library

If you have already demuxed the flowcell and are just rerunning alignment, you can use the `--no-demux` flag. If you are just rerunning post-alignment processing, you can use the `--no-align` flag. **Note** these options will overwrite existing data in the library folders!

### Build a reference

This can be performed from the login node.

```shell
use UGER
conda activate slideseq
build_ref --genome-name [GENOME_NAME] --reference-fasta path/to/genome.fasta --reference-gtf path/to/genes.gtf
```

This will submit a job to build the reference. You will get emails when the jobs start and end.

When the job is complete, the reference can be specified in the Google sheet as `/broad/macosko/reference/[GENOME_NAME]/[GENOME_NAME].fasta`

Existing references:

 - Mouse: `/broad/macosko/reference/GRCm39.103/GRCm39.103.fasta`
 - Human: `/broad/macosko/reference/GRCh38.102/GRCh38.102.fasta`

## Overview

### Goals

 - Demultiplex sequencing data from Illumina binary files
 - Match sequenced barcodes to a known set of slideseq puck barcodes to get spatial information
 - Align reads to a genome using STAR
 - Count features and summarize as a gene-cell matrix
 - Do all of this fast, efficiently, and reproducibly. When possible we'd like to disentangle the different steps so that we can perform them separately.

### Explicit Non-Goals

 - Won't work in arbitary environments. This pipeline is built to run at the Broad on the local computing infrastructure. Hopefully it is clear what is happening and the pipeline can be translated, but it is not a priority for us.
 - Doesn't allow much flexibility in the workflow. This makes the entire process simpler to configure and design.

### Currently unsupported or not yet implemented

- Plot of read 1 base distribution
- UP distance plots
- Bead types other than `180402`
- Multiple alignments of one library (e.g. `exonic` and also `exonic+intronic`)
- Probably other stuff...

## Requirements

The pipeline is designed to work from an Anaconda3 environment on the UGER cluster. We provide an `environment.yml` file with the requirements.

### Other requirements

The following tools are also needed:

 - `Picard` (from [the Broad](https://github.com/broadinstitute/picard))
 - `Drop-seq tools` (from [the McCarroll lab](https://github.com/broadinstitute/Drop-seq))
 - `Java-1.8` for the above
 - `Google-Cloud-SDK` to access our metadata tracking sheet

If you are running the pipeline on UGER, you do not need to install any of these tools.

TODO: currently the locations of these tools are hard-coded. They should be in a config file instead.

## Config tips

### Setting up Google-Cloud-SDK

You might need to authenticate the first time you run the pipeline. This just sets up your Google credentials on UGER so that it knows you have access to the worksheet.

```shell
use Google-Cloud-SDK
gcloud auth login

# ... follow instructions
```

### Modifying `.my.bashrc`

You can add the following to `.my.bashrc` to avoid typing them each time:

```shell
use UGER
use Google-Cloud-SDK
```

This might make login a tiny bit slower, but saves typing and errors from forgetting.
