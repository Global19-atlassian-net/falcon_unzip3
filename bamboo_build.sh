#!/bin/bash
source bamboo_setup.sh
set -vex

WHEELHOUSE="/mnt/software/p/python/wheelhouse/develop/"
pip3 install -v --user --find-links=${WHEELHOUSE} pytest networkx pysam msgpack pylint future intervaltree pypeflow falcon_kit

#pushd ../pypeFLOW
#pip3 install -v --user --no-deps --edit .
#popd
#
#pushd ../FALCON
#pip3 install -v --user --no-deps --edit .
#popd

python3 -c 'import falcon_kit; print falcon_kit.falcon'

pip3 -v install --user --no-deps --use-wheel --find-links=${WHEELHOUSE} .

pip3 install --user pytest pytest-cov pylint
export MY_TEST_FLAGS="-v -s --durations=0 --cov=falcon_unzip --cov-report=term-missing --cov-report=xml:coverage.xml --cov-branch"
make test
#sed -i -e 's@filename="@filename="./falcon_unzip/@g' coverage.xml

pylint --errors-only falcon_unzip

ls -larth

bash bamboo_wheel.sh
