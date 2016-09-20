#!/usr/bin/python
# -*- coding: utf-8 -*-

#
#  This file is part of the `bmp8` Python module
#
#  Copyright (c) 2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys

import xlrd
import openpyxl

def read_xls(xls_file, sheet = 0, csv_file = None,
    return_table = True):
    """
    Generic function to read MS Excel XLS file, and convert one sheet
    to CSV, or return as a list of lists
    """
    table = None
    try:
        book = xlrd.open_workbook(xls_file, on_demand = True)
        try:
            if type(sheet) is int:
                sheet = book.sheet_by_index(sheet)
            else:
                sheet = book.sheet_by_name(sheet)
        except xlrd.biffh.XLRDError:
            sheet = book.sheet_by_index(0)
        table = [[str(c.value) \
            for c in sheet.row(i)] \
            for i in xrange(sheet.nrows)]
    except IOError:
        sys.stdout.write('No such file: %s\n' % xls_file)
        sys.stdout.flush()
    except:
        book = openpyxl.load_workbook(filename = xls_file,
            read_only = True)
        try:
            if type(sheet) is int:
                sheet = book.worksheets[sheet]
            else:
                sheet = book.get_sheet_by_name(sheet)
        except:
            sheet = book.worksheets[0]
        cells = sheet.get_squared_range(1, 1,
            sheet.max_column, sheet.max_row)
        table = map(lambda row:
            map(lambda c:
                unicode(c.value),
                row
            ),
            cells
        )
    if csv_file:
        with open(csv_file, 'w') as csv:
            csv.write('\n'.join(['\t'.join(r) for r in table]))
    if not return_table:
        table = None
    if 'book' in locals() and hasattr(book, 'release_resources'):
        book.release_resources()
    return table
