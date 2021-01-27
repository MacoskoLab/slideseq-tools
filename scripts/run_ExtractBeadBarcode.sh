#!/bin/sh
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
puckcaller_folder=$2
scripts_folder=$3
output_folder=$4

echo ${submission}
echo ${matlab_path}
echo ${puckcaller_folder}
echo ${scripts_folder}
echo ${output_folder}

exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"

# in case of exit, set all permissions 
function finish {
    log_lib=${output_folder}/logs
    if [ -d "$log_lib" ]; then
echo "$log_lib" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi

}
trap finish SIGUSR2 EXIT

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
  eval "\"${scripts_folder}/puckcaller/ExtractBeadBarcode\"" ${puckcaller_folder}
fi
exit

