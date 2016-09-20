#!/bin/bash

rm -f src/*.pyc
rm -f src/bmp8/*.pyc
rm -rf src/bmp8.egg-info

sed -i 's/\([0-9]*\.[0-9]*\.\)\([0-9]*\)/echo \1$\(\(\2+1\)\)/ge' src/bmp8/__version__

python setup.py sdist
