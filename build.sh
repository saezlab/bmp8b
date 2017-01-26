#!/bin/bash

rm -f src/*.pyc
rm -f src/pex100/*.pyc
rm -rf src/pex100.egg-info

sed -i 's/\([0-9]*\.[0-9]*\.\)\([0-9]*\)/echo \1$\(\(\2+1\)\)/ge' src/pex100/__version__

python setup.py sdist
