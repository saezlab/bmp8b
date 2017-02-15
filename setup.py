#!/usr/bin/python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pex100` Python module
#
#  Copyright (c) 2016-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

__revision__ = "$Id$"
import sys
import os
from setuptools import setup, find_packages
import glob
import imp

with open(os.path.join('src', 'pex100', '__version__'), 'r') as v:
    __version__ = v.readline().strip()

metainfo = {
    'authors': {
    'Türei':('Dénes Türei','denes@ebi.ac.uk'),
    },
    'version': __version__,
    'license': 'GPLv3',
    'download_url': ['https://github.com/saezlab/pex100'],
    'url': ['https://github.com/saezlab/pex100'],
    'description': 'Phosphoassay data analysis',
    'platforms': ['Linux', 'Unix', 'MacOSX', 'Windows'],
    'keywords': ['data analysis', 'BMP8', 'noradrenaline', 'brown adipose tissue'],
    'classifiers': [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: GNU GPLv3',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Information Analysis',
    'Topic :: Scientific/Engineering :: Mathematics']
}

with open('README.rst') as f:
    readme = f.read()
with open('HISTORY.rst') as f:
    history = f.read()

deps = [
    'future'
]

extras = {
}

setup(
    name = 'pex100',
    version = __version__,
    maintainer = metainfo['authors']['Türei'][0],
    maintainer_email = metainfo['authors']['Türei'][1],
    author = metainfo['authors']['Türei'][0],
    author_email = metainfo['authors']['Türei'][1],
    long_description = readme + '\n\n' + history,
    keywords = metainfo['keywords'],
    description = metainfo['description'],
    license = metainfo['license'],
    platforms = metainfo['platforms'],
    url = metainfo['url'],
    download_url = metainfo['download_url'],
    classifiers = metainfo['classifiers'],
    # package installation
    package_dir = {'':'src'},
    packages = list(set(find_packages() + ['pex100'])),
    include_package_data = True,
    install_requires = deps,
    extras_require = extras,
    package_data = {'__version__': os.path.join('src', 'pex100', '__version__')}
)
