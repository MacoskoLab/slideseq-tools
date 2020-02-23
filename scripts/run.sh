#!/bin/bash

# This script is to call Python and C++ scripts

submission=$0

echo ${submission}

# submission_script, 'run_preparation', manifest_file, scripts_folder
if [[ "$1" == "run_preparation" ]]; then
	python "$3"/run_preparation.py "$2"
fi

# submission_script, 'run_processbarcodes', manifest_file, lane, scripts_folder
if [[ "$1" == "run_processbarcodes" ]]; then
	python "$4"/run_processbarcodes.py "$2" "$3"
fi

# submission_script, 'run_barcodes2sam', manifest_file, commandStr, lane, slice, scripts_folder
if [[ "$1" == "run_barcodes2sam" ]]; then
	python "$6"/run_barcodes2sam.py "$2" "$3" "$4" "$5"
fi

# submission_script, 'run_mergebarcodes', manifest_file, scripts_folder
if [[ "$1" == "run_mergebarcodes" ]]; then
	python "$3"/run_mergebarcodes.py "$2"
fi

# submission_script, 'run_alignment', manifest_file, library, lane, slice, barcode, scripts_folder
if [[ "$1" == "run_alignment" ]]; then
	python "$7"/run_alignment.py "$2" "$3" "$4" "$5" "$6"
fi

# submission_script, 'run_analysis', manifest_file, library, scripts_folder
if [[ "$1" == "run_analysis" ]]; then
	python "$4"/run_analysis.py "$2" "$3"
fi

# submission_script, 'run_analysis_spec', manifest_file, library, scripts_folder, locus_function_list
if [[ "$1" == "run_analysis_spec" ]]; then
	python "$4"/run_analysis_spec.py "$2" "$3" "$5"
fi

# submission_script, 'cmatcher', scripts_folder, bead_barcode_file, Illumina_barcode_file, output_distance_file, output_details_file, bead_type, hamming_distance
if [[ "$1" == "cmatcher" ]]; then
	"$2"/cmatcher "$3" "$4" "$5" "$6" "$7" "$8"
fi

# submission_script, 'run_cmatcher_combine', manifest_file, library, scripts_folder, locus_function_list
if [[ "$1" == "run_cmatcher_combine" ]]; then
	python "$4"/run_cmatcher_combine.py "$2" "$3" "$5"
fi

# submission_script, 'tag_matched_bam', manifest_file, library, lane, slice, barcode, locus_function_list, scripts_folder
if [[ "$1" == "tag_matched_bam" ]]; then
	python "$8"/tag_matched_bam.py "$2" "$3" "$4" "$5" "$6" "$7"
fi

# submission_script, 'filter_unmapped_bam', manifest_file, library, lane, slice, barcode, locus_function_list, scripts_folder
if [[ "$1" == "filter_unmapped_bam" ]]; then
	python "$8"/filter_unmapped_bam.py "$2" "$3" "$4" "$5" "$6" "$7"
fi

# submission_script, 'generate_plots', manifest_file, library, scripts_folder
if [[ "$1" == "generate_plots" ]]; then
	python "$4"/generate_plots.py "$2" "$3"
fi

# submission_script, 'generate_plots_cmatcher', manifest_file, library, scripts_folder, locus_function_list
if [[ "$1" == "generate_plots_cmatcher" ]]; then
	python "$4"/generate_plots_cmatcher.py "$2" "$3" "$5"
fi

