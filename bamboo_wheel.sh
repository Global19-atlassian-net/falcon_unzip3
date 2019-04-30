#!/bin/bash
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module purge
module load gcc
module load ccache
set -vex
ls -larth ..
ls -larth
pwd

export WHEELHOUSE=./wheelhouse

# Give everybody read/write access.
umask 0000


module load python/3.7.3
make wheel

# Select export dir based on Bamboo branch, but only for develop and master.
case "${bamboo_planRepository_branchName}" in
  develop|master)
    WHEELHOUSE="/mnt/software/p/python/wheelhouse/${bamboo_planRepository_branchName}/"
    rsync -av ./wheelhouse/ ${WHEELHOUSE}
    ;;
  *)
    ;;
esac
