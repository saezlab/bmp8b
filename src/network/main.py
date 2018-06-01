#!/usr/bin/python
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2018 - EMBL
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

import imp
import sys
import openpyxl
import pandas as pd
import scipy as sp
import scipy.stats

import pypath

class NetworkAnalysis(object):
    
    basedir   = '/home/denes/documents/viv'
    qpcrdir   = 'qpcr'
    qpcrfiles = [
        'BTAX_4weeks_RNA_moretissues_3.xlsx',
        '5weeks_BTAX_Bmp8bKO_RNAanalysis_allBK.xlsx'
    ]
    
    def __init__(self):
        
        pass
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def main(self):
        
        self.read_qpcr()
        self.init_mapper()
        self.get_mouse_human_dict()
        self.map_qpcr_proteins()
        self.qpcr_to_human()
        self.qpcr_process()
    
    def read_qpcr(self):
        
        drop = {'36b4', '18s'}
        
        cols = [
            'Q', 'P', 'R', 'S', 'T', 'U', 'V',
            'W', 'X', 'Y', 'Z', 'AA', 'AB'
        ]
        qpcr = []
        
        for fn in self.qpcrfiles:
            
            path = '/'.join([self.basedir, self.qpcrdir, fn])
            
            book = openpyxl.load_workbook(
                filename = path,
                read_only = True,
                data_only = True
            )
            
            prg = pypath.progress.Progress(
                len(book.worksheets), 'Reading qPCR data', 1
            )
            
            weeks = '5w' if fn.startswith('5weeks') else '4w'
            
            for sheet in book.worksheets:
                
                prg.step(status = 'reading sheet `%s`' % sheet.title)
                
                if not sheet.title.endswith('analysis'): continue
                
                protein = sheet.title.split()[0].capitalize()
                
                if protein in drop: continue
                
                tbl = list(sheet.get_squared_range(16, 1, 27, sheet.max_row))
                
                for ci in xrange(12):
                    
                    tissue = 'unknown'
                    
                    for ri, r in enumerate(tbl):
                        
                        val = r[ci].value
                        
                        if not val or val == '#DIV/0!': continue
                        
                        if type(val) is str:
                            
                            tissue = val
                            
                        elif r[ci].fill.bgColor.rgb == 'FFCCFFFF':
                            
                            try:
                                
                                qpcr.append([
                                    protein,
                                    weeks,
                                    tissue,
                                    cols[ci],
                                    ri + 1,
                                    val
                                ])
                                
                            except:
                                
                                sys.stdout.write(
                                    (
                                        '\t > Not a number: %s (sheet %s,'
                                        'col %s row %u)\n' % (
                                            val, sheet.title, cols[ci], ri + 1
                                        )
                                    )
                                )
            
            prg.terminate()
            
        self.qpcr = qpcr
    
    def init_mapper(self):
        
        self.mapper = pypath.mapping.Mapper()
    
    def get_mouse_human_dict(self):
        
        self.mh = pypath.dataio.homologene_uniprot_dict(
            10090, 9606, mapper = self.mapper
        )
    
    def map_qpcr_proteins(self):
        
        nonstandard = {
            'Actr2b':      'P61161',
            'Adiponectin': 'Q60994',
            'Alpha1ar':    'P97718',
            'Alpha2ar':    'Q01338',
            'Beta2adr':    'P18762',
            'Beta3ar':     'P25962',
            'Betaact':     'P60710',
            'Betaactin':   'P60710',
            'Bmp8bendog':  'P55105',
            'Bmp8bko':     'P55105',
            'Bmp8btg':     'P55105',
            'Bmp8btotal':  'P55105',
            'Bmp8btotal2': 'P55105',
            'Cd206':       'Q61830',
            'Cd31':        'Q08481',
            'Collagen4a1': 'P02463',
            'Erbb42':      'Q61527',
            'F480':        'Q61549',
            'Irisin':      'Q8K4Z2',
            'Lr11':        'O88307',
            'Nrg4mix':     'Q9WTX4',
            'Nt3':         'P20181',
            'Pdgfrbeta':   'P05622',
            'Pparg1':      'P37238',
            'Pparg2':      'P37238',
            'Reelin':      'Q60841',
            'Srebp1c':     'Q9WTN3',
            'Trka':        'Q3UFB7',
            'Trkc':        'Q6VNS1',
            'Tsp':         'P35441',
            'Vegfa120':    'Q00731',
            'Vegfa164':    'Q00731',
            'Vegfa166':    'Q00731',
            'Vegfa188':    'Q00731',
            'Vegfr2':      'P35918'
        }
        
        self.qpcr_unmapped = set()
        
        for q in self.qpcr:
            
            u = self.mapper.map_name(q[0], 'genesymbol', 'uniprot', ncbi_tax_id = 10090)
            
            if not u:
                
                if q[0] in nonstandard:
                    
                    u = [nonstandard[q[0]]]
                    
                else:
                    
                    u = [None]
                    self.qpcr_unmapped.add(q[0])
            
            if len(q) == 7:
                
                q[6] = u[0]
                
            else:
                
                q.append(u[0])
    
    def qpcr_to_human(self):
        
        nonstandard = {
            'Vegfr2':      'P35968',
            'Ngfr':        'P08138',
            'Bmp8b':       'P34820',
            'Tsp':         'P07996',
            'Bmp8btg':     'P34820',
            'Fgfr3':       'P22607',
            'Bmp8bko':     'P34820',
            'Bmp8bendog':  'P34820',
            'Ang2':        'P03950',
            'Ttn':         'Q8WZ42',
            'Nrg2':        'O14511',
            'Bmp8btotal':  'P34820',
            'Bmp8btotal2': 'P34820',
            'Alpha2ar':    'P08913'
        }
        
        self.qpcr_human_unmapped = set()
        
        for q in self.qpcr:
            
            h = (
                self.mh[q[6]][0]
                if q[6] in self.mh and self.mh[q[6]] else
                nonstandard[q[0]] if q[0] in nonstandard else
                None
            )
            
            if not h:
                
                self.qpcr_human_unmapped.add(q[0])
            
            if len(q) == 8:
                
                q[7] = h
                
            else:
                
                q.append(h)
    
    def qpcr_process(self):
        
        self.qpcrdf = pd.DataFrame(
            self.qpcr,
            columns = [
                'gs_mouse', 'weeks', 'tissue', 'excelcol',
                'excelrow', 'raw', 'up_mouse', 'up_human'
            ]
        )
        
        self.qpcrdf.to_csv('qpcr.tsv', sep = '\t', index = False)
    
    
    def qpcr_tests(self):
        
        result = []
        
        conds = self.qpcrdf[['gs_mouse', 'tissue', 'weeks']].drop_duplicates()
        
        for i, c in conds.iterrows():
            
            if c.tissue.startswith('WT'): continue
            
            t_other = c.tissue.replace('KO', 'WT').replace('TG', 'WT')
            
            this = self.qpcrdf[
                (self.qpcrdf.tissue   == c.tissue) &
                (self.qpcrdf.weeks    == c.weeks)  &
                (self.qpcrdf.gs_mouse == c.gs_mouse)
            ]
            other = this = self.qpcrdf[
                (self.qpcrdf.tissue   == t_other) &
                (self.qpcrdf.weeks    == c.weeks)  &
                (self.qpcrdf.gs_mouse == c.gs_mouse)
            ]
            
            if not other.shape[0]:
                
                sys.stdout.write('\t > No WT found for %s %s %s' %
                    (c.gs_mouse, c.tissue, c.weeks)
                )
                continue
            
            t         = sp.stats(this.zscore, other.zscore)
            meanthis  = this.zscore.mean()
            meanother = other.zscore.mean()
            meandiff  = meanthis - meanother
            meanfc    = self.fold_change(meanthis, meanother)
            w         = sp.stats.mannwhitneyu(this.zscore, other.zscore)
            medthis   = this.zscore.median()
            medother  = other.zscore.median()
            meddiff   = medthis - medother
            medfc     = self.fold_change(medthis, medother)
            
            result.append([
                c.gs_mouse,
                list(this.up_mouse)[0],
                list(this.up_human)[0],
                t_other,
                c.tissue,
                c.weeks,
                meanthis,
                meanother,
                meandiff,
                meanfc,
                t.pvalue,
                medthis,
                medother,
                meddiff,
                medfc,
                w.pvalue
            ])
        
        self.qpcrcompdf = pd.DataFrame(
            result,
            columns = [
                'gs_mouse', 'up_mouse', 'up_human', 'tissue1', 'tissue2',
                'weeks', 'mean1', 'mean2', 'meandiff', 'meanfc', 't_pval',
                'med1', 'med2', 'meddiff', 'medfc', 'mw_pval'
            ]
        )
    
    def read_qpcr_tests(self):
        
        self.batexpdf = pd.DataFrame.from_csv(
            'bat.qpcr.tests.tsv', sep = '\t', index_col = None
        )
    
    def read_phospho(self):
        
        self.bmp8bphos = pd.DataFrame.from_csv(
            '../fc_BMP8b_none.csv', sep = '\t', index_col = None
        )
    
    def read_kinact(self):
        
        result = []
        
        with open('kinact_top.tsv', 'r') as fp:
            
            for l in fp:
                
                if l.startswith('='): continue
                
                if 'Ksea' in l or 'Networkin' in l:
                    
                    method = l.strip()
                
                l = l.split('\t')
                
                try:
                    
                    score = float(l[1])
                    pval  = float(l[2])
                    result.append([
                        method,
                        l[0].strip(),
                        score,
                        pval
                    ])
                    
                except:
                    
                    pass
        
        self.kinact = pd.DataFrame(
            result, columns = ['method', 'protein', 'score', 'kinact_p']
        )
        
        self.kinact['uniprot'] = [
            u[0] if u else None
            for u in
            [
                self.mapper.map_name(
                    p, 'genesymbol', 'uniprot', ncbi_tax_id = 9606
                )
                for p in
                self.kinact.protein
            ]
        ]
    
    def load_tfregulons(self):
        
        targets = set(self.batexpdf.up_human)
        tfofinterest = []
        
        tfreg = pypath.dataio.get_tfregulons(
            levels = {'A', 'B', 'C', 'D'}
        )
        
        for r in tfreg:
            
            tfup = self.mapper.map_name(
                r[0], 'genesymbol', 'uniprot',
                ncbi_tax_id = 9606
            )
            tgup = self.mapper.map_name(
                r[1], 'genesymbol', 'uniprot',
                ncbi_tax_id = 9606
            )
            
            if tgup and tfup and tgup[0] in targets:
                
                r.extend([tfup[0], tgup[0]])
                tfofinterest.append(r)
        
        self.tfreg = tfofinterest
    
    def load_network(self):
        
        self.pa = pypath.PyPath()
        self.pa.mapper = self.mapper
        self.pa.init_network(pypath.data_formats.pathway)
        #{
            #'signor': pypath.data_formats.pathway['signor'],
            #'slk3':   pypath.data_formats.pathway['signalink3'],
            #'spike':  pypath.data_formats.pathway['spike']
        #})
        #self.pa.load_resources(pypath.data_formats.ptm_misc)
        self.pa.third_source_directions()
        self.pa.get_directed()
    
    def collect_proteins(self):
        
        def net_lookup(ups):
            
            return dict(
                (v['name'], v.index)
                for v in 
                self.pa.dups(ups)
            )
        
        self.phos = set(self.bmp8bphos.uniprot)
        self.kin  = set(self.kinact.uniprot)
        self.tf   = set(t[-2] for t in self.tfreg)
        
        self.dphos = net_lookup(self.phos)
        self.dkin  = net_lookup(self.kin)
        self.dtf   = net_lookup(self.tf)
    
    def random_walks(self, sources, targets, **kwargs):
        
        result = []
        
        prg = pypath.progress.Progress(
            len(sources),
            'Random walks',
            1
        )
        
        i = 0
        
        for s, si in sources.items():
            
            i += 1
            
            prg.step()
            rwr = self.pa.random_walk_with_return(q = si, **kwargs)
            
            for t, ti in targets.items():
                
                result.append([s, t, rwr[ti]])
        
        prg.terminate()
        
        return result
    
    def shortest_paths(self, sources, targets):
        
        result = []
        
        srev = dict((v, k) for k, v in sources.items())
        trev = dict((v, k) for k, v in targets.items())
        slst = sorted(srev.keys())
        tlst = sorted(trev.keys())
        
        sp = self.pa.dgraph.shortest_paths(source = slst, target = tlst)
        
        for si, sid in enumerate(slst):
            
            for ti, tid in enumerate(tlst):
                
                result.append([
                    srev[sid],
                    trev[tid],
                    sp[si][ti]
                ])
        
        self.sp = result
    
    def export_sp(self):
        
        pd.DataFrame(self.sp, columns = ['phos', 'tf', 'splen']).to_csv(
            'sp_phos_tf.tsv', sep = '\t', index = False
        )
    
    def export_rw(self, data):
        
        pd.DataFrame(data, columns = ['phos', 'tf', 'dist']).to_csv(
            'rw_spike4.tsv', sep = '\t', index = False
        )
    
    def export_tf(self):
        
        pd.DataFrame(
                self.tfreg,
                columns = [
                    'gstf', 'gstg', 'sign', 'level',
                    'curated', 'chipseq', 'gex', 'tfbs',
                    'scurated', 'schipseq', 'sgex', 'stfbs',
                    'sall', 'uptf', 'uptg'
                ]
            ).to_csv(
                'tfreg.tsv', sep = '\t', index = False
            )
    
    def export_phos(self):
        
        self.bmp8bphos.to_csv('phos.tsv', sep = '\t', index = False)
    
    def export_exp(self):
        
        self.batexpdf.to_csv('exp.tsv', sep = '\t', index = False)
    
    def export_uniprot_gs(self):
        
        pd.DataFrame([
                [ups[0], gs]
                for gs, ups in
                self.mapper.tables[9606][
                    ('genesymbol', 'swissprot')
                ].mapping['to'].items()
                if ups
            ],
            columns = ['uniprot', 'genesymbol']
        ).to_csv('up_gs.tsv', sep = '\t', index = False)
    
    @staticmethod
    def fold_change(val, ctrl):
        """
        Calculates the fold change as {value}/{control} if 
        {control} < {value} and -1 * {control}/{value} otherwise.
        """
        return \
            0.0 if val == ctrl \
            else  val / ctrl if ctrl < val \
            else -1 * ctrl / val
