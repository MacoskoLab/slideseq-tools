# Slide-seq tools
Tools for analyzing Slide-seq data including building genome reference, aligning reads to reference genome, generating feature-barcode matrices, performing gene expression analysis and matching data from in situ sequencing and indexing of barcodes with short read sequencing data.

## Requirement

Several public tools need to be pre-installed:
1) `Drop-seq tools`: https://github.com/broadinstitute/Drop-seq
2) `Picard`: https://broadinstitute.github.io/picard/
3) `STAR`: https://github.com/alexdobin/STAR
4) `Java`
5) `Samtools`
6) `gcc/g++`
7) `Python` (prefer 3.6 or above)

Several Python packages need to be installed for calculation and ploting:
1) `numpy`
2) `pandas`
3) `plotnine`
4) `matplotlib`
    
## Build genome reference

The Slide-seq tools need a genome reference in specific format for alignment and analysis. You could build a genome reference based on input fasta and gtf files. 

Command: 
```
build_reference.py manifest_file
```

Check `example.buildreference.txt` for manifest file format

## Run the Slide-seq tools

Add below commands into `run.sh` and `build_reference.sh` or your `bashrc` file (command might be different on your system):
```
use Java-1.8
use .samtools-1.7
use Python-3.6
```

Compile CMatcher (command might be different on your system): 
```
g++ -std=c++11 -o cmatcher cmatcher.cpp
```

Submit a request to the pipeline: 
```
python run_pipeline.py manifest_file
```

Notice: 
1) Check `example.manifest.txt` for manifest file format
2) An email from slideseq@gmail.com will be sent to you if email_address is specified in the manifest file when the submission is received, the workflow finishes, and/or any job fails.
3) In order to speed up the process of NovaSeq data and NovaSeq S4 data, the Slide-seq tools split each lane into a few slices, run the alignment steps on the slices parallelly and combine the alignment outputs together. 
4) See `user_doc.txt` for detailed usage of the Slide-seq tools. 

