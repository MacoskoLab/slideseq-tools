import sys

# need to pass this to the qsub scripts so they activate the right thing
conda_env = sys.executable.split('/')[-3]
