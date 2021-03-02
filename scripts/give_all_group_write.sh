#!/bin/bash

# Instead of give_group_write.sh which gave group permissions for individual
# files, this checks in the main output folder for any folders without group
# permissions. Instead of an exhaustive `find`, uses -maxdepth 2 to catch any
# sub-directories without g+w perms (/broad/macosko/data/workflows/flowcell has
# that) and then chmod -R to go deeper

find /broad/macosko/data/workflows/flowcell /broad/macosko/data/libraries -maxdepth 2 -user `whoami` ! -perm -g=rw -exec chmod -R g+rwx "{}" \;
