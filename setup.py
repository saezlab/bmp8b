# -*- coding: utf-8 -*-

__revision__ = "$Id$"
import sys
import os
from setuptools import setup, find_packages
import glob
import imp

with open(os.path.join('src', 'csim', '__version__'), 'r') as v:
    __version__ = v.readline().strip()

metainfo = {
    'authors': {
    'Türei':('Dénes Türei','denes@ebi.ac.uk'),
    },
    'version': __version__,
    'license': 'GPLv3',
    'download_url': ['http://bitbucket.org/deeenes/bmp8'],
    'url': ['http://bitbucket.org/deeenes/bmp8'],
    'description': 'Phosphoassay data analysis',
    'platforms': ['Linux', 'Unix', 'MacOSX', 'Windows'],
    'keywords': ['data analysis', 'BMP8', 'noradrenaline', 'brown adipose tissue'],
    'classifiers': [
    'Development Status :: 2 - Pre-Alpha',
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
    name = 'bmp8',
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
    packages = list(set(find_packages() + ['csim'])),
    include_package_data = True,
    install_requires = deps,
    extras_requre = extras,
    package_data = {'__version__': os.path.join('src', 'csim', '__version__')}
)
