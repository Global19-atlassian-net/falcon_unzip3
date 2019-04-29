type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module purge
module load gcc
module load ccache
module load python/3.7.3
which python
which pip3 # assumed

export PYTHONUSERBASE=$(pwd)/LOCAL3
mkdir -p LOCAL3/

export PATH=$PYTHONUSERBASE/bin:$PATH
