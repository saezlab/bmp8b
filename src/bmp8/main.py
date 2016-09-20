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
import re

import xlrd
import openpyxl

import numpy as np

class Bmp8(object):
    
    def __init__(self):
        self.dirBase = '..'
        self.fnXlsAssayData = os.path.join(
            self.dirBase,
            'AssayData-Cambridge-Peirce-150515.xlsx'
        )
        self.fnIdMapping = os.path.join(
            self.dirBase,
            'PEX100_Layout.csv'
        )
        self.reAnnot = re.compile(r'([\-\s/\.,\(\)\+A-Za-z0-9]{2,}) '\
            r'\(([A-Z][a-z]+)-?([A-Za-z0-9/]*)\)')
        self.reRes = re.compile(r'([A-Z]?[a-z]*)([0-9]+)')
        self.dAaletters = {
            'Thr': 'T',
            'Tyr': 'Y',
            'Ser': 'S'
        }
    
    def read_xls(self, xls_file, sheet = 0, csv_file = None,
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
    
    def read_ids(self):
        pass
    
    def ll_table_slice(self, table, col, row, ncol = None, nrow = None):
        cend = None if ncol is None else col + ncol
        rend = None if nrow is None else row + nrow
        return list(map(lambda l: l[col:cend], table[row:rend]))
    
    def psites_diffs(self):
        def get_residues(res):
            ress = self.reRes.findall(res)
            aa = ''
            sres = []
            for r in ress:
                if r[0] in self.dAaletters:
                    aa = self.dAaletters[r[0]]
                sres.append((aa, int(r[1])))
            return tuple(sres)
        
        self.llAssayData = self.read_xls(self.fnXlsAssayData)
        self.llAssayData = self.ll_table_slice(self.llAssayData, 16, 9, 7, 582)
        self.aAnnots = \
            np.array(
                list(
                    map(
                        lambda m:
                            [m[0], m[1], get_residues(m[2])],
                        map(
                            lambda l:
                                list(
                                    self.reAnnot.match(l[0]).groups(0)
                                ),
                            self.llAssayData
                        )
                    )
                ),
                dtype = np.object
            )
        self.aAssayData = np.array(list(map(lambda l:
                            list(map(float, l[3:])),
                        self.llAssayData)))
