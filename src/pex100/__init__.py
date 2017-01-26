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

import os as _os

from bmp8.main import *

_ROOT = _os.path.abspath(_os.path.dirname(__file__))

def _get_version():
    with open(_os.path.join(_ROOT, '__version__'), 'r') as v:
        return tuple([int(i) for i in v.read().strip().split('.')])

_MAJOR, _MINOR, _MICRO = _get_version()
__version__ = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
__release__ = '%d.%d' % (_MAJOR, _MINOR)
