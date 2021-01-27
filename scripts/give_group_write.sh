#!/bin/bash

# This is a helper script to change permissions to give group 
# write access to all files and subdirs in provided directory

target_dir=$1

chmod -R g+wX $target_dir

