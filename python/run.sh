#!/bin/bash

# This script is to call Python or C++ scripts

submission=$0

echo ${submission}

if [[ "$1" == "run_preparation" ]]; then
	python "$3"/run_preparation.py "$2"
fi

if [[ "$1" == "run_processbarcodes" ]]; then
	python "$4"/run_processbarcodes.py "$2" "$3"
fi

if [[ "$1" == "run_barcodes2sam" ]]; then
	python "$6"/run_barcodes2sam.py "$2" "$3" "$4" "$5"
fi

if [[ "$1" == "run_mergebarcodes" ]]; then
	python "$3"/run_mergebarcodes.py "$2"
fi

if [[ "$1" == "run_alignment" ]]; then
	python "$6"/run_alignment.py "$2" "$3" "$4" "$5"
fi

if [[ "$1" == "run_analysis" ]]; then
	python "$4"/run_analysis.py "$2" "$3"
fi

if [[ "$1" == "run_analysis_spec" ]]; then
	python "$4"/run_analysis_spec.py "$2" "$3" "$5"
fi

if [[ "$1" == "cmatcher" ]]; then
	"$2"/cmatcher "$3" "$4" "$5" "$6" "$7" "$8"
fi

if [[ "$1" == "run_cmatcher_combine" ]]; then
	python "$4"/run_cmatcher_combine.py "$2" "$3" "$5"
fi

if [[ "$1" == "tag_matched_bam" ]]; then
	python "$7"/tag_matched_bam.py "$2" "$3" "$4" "$5" "$6"
fi

if [[ "$1" == "filter_unmapped_bam" ]]; then
	python "$7"/filter_unmapped_bam.py "$2" "$3" "$4" "$5" "$6"
fi

if [[ "$1" == "generate_plots" ]]; then
	python "$4"/generate_plots.py "$2" "$3"
fi

if [[ "$1" == "generate_plots_cmatcher" ]]; then
	python "$4"/generate_plots_cmatcher.py "$2" "$3" "$5"
fi

