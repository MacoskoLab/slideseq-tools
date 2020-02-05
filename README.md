# slideseq-tools
A pipeline for processing Slide-seq data, which aligns reads to reference genome, generates feature-barcode matrices, performs gene expression analysis and matches data from in situ sequencing and indexing of barcodes with short read sequencing data.

## Requirement

The pipeline needs several public tools pre-installed:
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

The Slide-seq pipeline needs a genome reference in specific format for alignment and analysis. You could build a genome reference based on input fasta and gtf files. 

Command: 
```
build_reference.py manifest_file
```

Check `example.buildreference.txt` for manifest file format

## Run the pipeline

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
3) In order to speed up the process of NovaSeq data and NovaSeq S4 data, the pipeline splits each lane into a few slices, and runs the alignment steps on the slices parallelly and combines the alignment outputs together. 
4) See `user_doc.txt` for detailed usage of the Slide-seq pipeline. 

