#!/bin/sh
#$ -l h_vmem=40g
#$ -l h_rt=5:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# script for execution of deployed applications
#
# Sets up the MATLAB Runtime environment for the current $ARCH and executes 
# the specified command.
#

source /broad/software/scripts/useuse
reuse UGER

#/broad/software/nonfree/Linux/redhat_7_x86_64/pkgs/matlab_2019a

submission=$0
matlab_path=$1
scripts_folder=$2
bcb_file=$3
dge_file=$4
unique_bci_file=$5
genename_file=$6
location_file=$7
puckcaller_folder=$8
output_folder=$9

echo ${submission}
echo ${matlab_path}

exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"

if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
  export LD_LIBRARY_PATH;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=$1
      args="${args} \"${token}\"" 
      shift
  done
  eval "\"${scripts_folder}/WriteBijectiveMapping/WriteBijectiveMapping\"" ${bcb_file} ${dge_file} ${unique_bci_file} ${genename_file} ${location_file} ${puckcaller_folder}
  rm ${bcb_file}
  rm ${dge_file}
fi
exit

