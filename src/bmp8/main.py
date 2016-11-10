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
import imp
import re
import itertools

import xlrd
import openpyxl

import numpy as np

import pypath

if 'unicode' not in __builtins__:
    unicode = str

class Bmp8(object):
    
    def __init__(self
                 dirBase = '..',
                 fnIdMapping = 'PEX100_Layout.csv',
                 fnXlsFullData = 'AssayData-Cambridge-Peirce-150515.xlsx'):
        self.dirBase = dirBase
        self.fnXlsFullData = os.path.join(
            self.dirBase,
            fnXlsFullData
        )
        self.fnIdMapping = os.path.join(
            self.dirBase,
            fnIdMapping
        )
        self.reAnnot = re.compile(r'([\-\s/\.,\(\)\+A-Za-z0-9]{2,}) '\
            r'\(([A-Z][a-z]+)-?([A-Za-z0-9/]*)\)')
        self.reRes = re.compile(r'([A-Z]?[a-z]*)([0-9]+)')
        self.dAaletters = {
            'Thr': 'T',
            'Tyr': 'Y',
            'Ser': 'S'
        }
        self.pa = pypath.PyPath(ncbi_tax_id = 10090)
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    #
    # Loading experimental data.
    #
    
    def read_tables(self):
        """
        Reads all data tables.
        """
        self.read_signals()
        self.read_coeffvar()
        self.read_normalized()
        self.read_psites_diffs()
    
    def read_psites_diffs(self):
        """
        Reads the measured values for each phosphosite at each condition.
        """
        
        self.read_data('Psite', 16, 9, 7, 582, 3)
    
    def read_signals(self):
        """
        Reads the mean of duplicates signal values.
        """
        
        self.read_data('Signal', 0, 9, 5, 1318, 1)
    
    def read_coeffvar(self):
        """
        Reads coefficient of variation values for each pair of spots.
        """
        
        self.read_data('CVar', 0, 9, 10, 1318, 6)
    
    def read_normalized(self):
        """
        Reads the values normalized to median.
        """
        
        self.read_data('Norm', 0, 9, 15, 1318, 11)
    
    def read_data(self, name, col, row, ncol, nrow, vcol):
        """
        Reads a data table with its annotations.
        """
        
        if not hasattr(self, 'llFullData') or self.llFullData is None:
            self.read_raw()
        
        name = name.capitalize()
        
        attr = 'll%s' % name
        
        setattr(self, attr,
                self.ll_table_slice(self.fullData, row, col, nrow, ncol))
        
        data, annot = self.get_arrays(getattr(self, attr), vcol)
        
        setattr(self, 'a%sAnnot', annot)
        setattr(self, 'a%sData', data)
    
    def read_raw(self):
        """
        Reads the full xls table.
        """
        self.llFullData = self.read_xls(self.fnXlsFullData)
    
    def read_xls(self, xls_file, sheet = 0, csv_file = None,
        return_table = True):
        """
        Generic function to read MS Excel XLS file. Converts one sheet
        to CSV, or returns as a list of lists
        """
        table = []
        try:
            book = xlrd.open_workbook(xls_file, on_demand = True)
            try:
                if type(sheet) is int:
                    sheet = book.sheet_by_index(sheet)
                else:
                    sheet = book.sheet_by_name(sheet)
            except xlrd.biffh.XLRDError:
                sheet = book.sheet_by_index(0)
            table = [[unicode(c.value) \
                for c in sheet.row(i)] \
                for i in xrange(sheet.nrows)]
        except IOError:
            sys.stdout.write('No such file: %s\n' % xls_file)
            sys.stdout.flush()
        except:
            try:
                book = openpyxl.load_workbook(filename = xls_file,
                    read_only = True)
            except:
                sys.stdout.write('\tCould not open xls: %s\n' % xls_file)
                if not os.path.exists(xls_file):
                    sys.stdout.write('\tFile does not exist.\n')
                sys.stdout.flush()
            try:
                if type(sheet) is int:
                    sheet = book.worksheets[sheet]
                else:
                    sheet = book[sheet]
            except:
                sheet = book.worksheets[0]
            cells = sheet.get_squared_range(1, 1,
                sheet.max_column, sheet.max_row)
            table = \
                list(
                    map(lambda row:
                        list(
                            map(lambda c:
                                unicode(c.value),
                                row
                            )
                        ),
                        cells
                    )
                )
        if csv_file:
            with open(csv_file, 'w') as csv:
                csv.write('\n'.join(['\t'.join(r) for r in table]))
        if not return_table:
            table = None
        if 'book' in locals() and hasattr(book, 'release_resources'):
            book.release_resources()
        return table
    
    def ll_table_slice(self, table, col, row, ncol = None, nrow = None):
        """
        Returns a slice from a list-of-lists table, starting at (row, col),
        height and width of (nrow, ncol).
        """
        cend = None if ncol is None else col + ncol
        rend = None if nrow is None else row + nrow
        return list(map(lambda l: list(l)[col:cend], table[row:rend]))
    
    def get_arrays(self, llTable, vcol):
        """
        Returns annotations and data as numpy arrays.
        Requires list-of-lists table and
        the column number where data values start.
        """
        
        def get_residues(res):
            """
            Matches residues in string,
            returns tuples of tuples with residue names and numbers.
            """
            ress = self.reRes.findall(res)
            aa = ''
            sres = []
            for r in ress:
                if r[0] in self.dAaletters:
                    aa = self.dAaletters[r[0]]
                sres.append((aa, int(r[1])))
            return tuple(sres)
        
        return \
            list(
                map(
                    lambda a:
                        np.array(a[0], a[1]),
                    # zip data and dtype
                    zip(
                        # zip annotations and values
                        zip(
                            *list(
                                # chain rows & residues
                                itertools.chain(
                                    map(
                                        # iterate rows
                                        lambda ll:
                                            list(
                                                map(
                                                    # iterate residues
                                                    lambda res:
                                                        (
                                                            # annotation
                                                            [
                                                                ll[0][0],
                                                                ll[0][1],
                                                                res[0],
                                                                res[1]
                                                            ],
                                                            # values
                                                            ll[1]
                                                        ),
                                                    get_residues(ll[0][2])
                                                )
                                            ),
                                            [m[0], m[1], get_residues(m[2])],
                                        map(
                                            # for each row, match annotation
                                            # and convert data to float
                                            lambda l:
                                                (
                                                    self.reAnnot.\
                                                        match(l[0]).groups(0),
                                                    list(map(float, l[3:]))
                                                ),
                                            # starting from list of lists
                                            llTable
                                        )
                                    )
                                )
                            )
                        ),
                        # dtypes for annot. and data
                        [np.object, np.float]
                    )
                )
            )
    
    #
    # Loading database data.
    #
    
    def load_ptms(self):
        """
        Obtains a list of all kinase-substrate interaction from
        Signor, phosphoELM, PhosphoSite and dbPTM.
        Creates a dict `dPhosDb` where (kinase, substrate) UniProts are keys
        and list of PTMs are values.
        """
        self.dPhosDb = {}
        
        lKinSub = itertools.chain(
            self.pa.load_signor_ptms(return_raw=True),
            self.pa.load_phosphoelm(return_raw=True),
            self.pa.load_dbptm(return_raw=True),
            self.pa.load_psite_phos(return_raw=True)
        )
        
        for ksub in lKinSub:
            tKs = (ksub.domain.protein, ksub.ptm.protein)
            if tKs not in self.dPhosDb:
                self.dPhosDb[tKs] = []
            self.dPhosDb[tKs].append(ksub)
