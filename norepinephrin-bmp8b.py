#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
    #After satisfying identification of ChEMBL IDs for 
    #compounds in GDSC screening, this script looks up
    #in the network their nominal tarets and other known 
    #targets form ChEMBL binding assays.
#

# generic modules #

import sys
import os
import re
import cPickle as pickle
import copy
import igraph
import louvain
import unicodedata

# stats and plotting modules #

import math
import ranking
import numpy as np
from numpy.random import randn
import pandas as pd
from scipy import stats
from scipy.misc import comb
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import hcluster as hc2
import matplotlib.patches as mpatches

# from bioigraph #

import pypath
from pypath import chembl
from pypath.common import *
import _sensitivity as sens
from pypath import progress, dataio, intera

# 0 ## Constants

assayfile = '/home/denes/Dokumentumok/viv/AssayData-Cambridge-Peirce-150515.csv'
pairedfile = '/home/denes/Dokumentumok/viv/AssayDataPaired.csv'
idfile = '/home/denes/Dokumentumok/viv/PEX100_Layout.csv'
ptmlog = '/home/denes/Dokumentumok/viv/matching_phosphosites.log'
chembl_mysql = (None, 'chembl_ebi')
mysql_gelati = (None,'mapping_gelati')
rebr = re.compile(r'\(([^\)]{3,})\)')
renum = re.compile(r'[0-9]+')
remod = re.compile(r'([A-Z][a-z]{2})([0-9]+)')
nores = re.compile(r'([A-Z][a-z]{2})([0-9]+/)([0-9]+)')
# reptm = re.compile(r'([A-Z][a-z]+)-([A-Z][a-z]{2})([0-9]+)')
reptm = re.compile(r'[A-Z][a-z]+-[A-Z][a-z]{2}[0-9A-Za-z/]+')
ptmtyp = {
    'Phospho': 'phosphorylation'
}
aarev = dict([(a3, a1) for a1, a3 in aacodes.iteritems()])
names_dict = {}
values_dict = {}
paired = {}
psitep = {}

all_sprot = set(dataio.all_uniprots(swissprot = 'yes'))

# 1 ## Initializing

net = pypath.PyPath(9606)
net.init_network(pfile = 'cache/plus_phospho.pickle')
net.genesymbol_labels()
sprot_seq = dataio.swissprot_seq(isoforms = True)
ps = dataio.get_psite_p(organism = 'human')
for pps in ps:
    if pps.protein not in psitep:
        psitep[pps.protein] = []
    psitep[pps.protein].append(pps)

with open(assayfile, 'r') as f:
    proteins = [[n[0].split('(')[0].strip(), 
                ''.join([nores.sub('\g<1>\g<2>\g<1>\g<3>', inbr) \
                    for inbr in reptm.findall(n[0])]), 
                n[0].strip()] + \
                [int(i) for i in n[1:5]] + \
                [float(i) for i in n[6:10]] + 
                [float(i) for i in n[12:15]] \
        for n in \
            [l.split('\t') \
                for l in f.read().split('\n') if len(l) > 0][9:-4]
        ]

with open(pairedfile, 'r') as f:
    paired = dict([(ll[0].strip(), tuple([float(v) for v in ll[3:7]])) \
        for ll in [l.split('\t') for l in f.read().split('\n')][:-1]])

mapper = net.mapper if 'net' in globals() else pypath.mapping.Mapper()

with open(idfile, 'r') as f:
    names_dict = dict([(ll[1].split('(')[0].strip(),
        flatList( \
        [[s] if s in all_sprot else mapper.map_name(s, 'uniprot', 'uniprot') for s in \
        [''.join([c for c in unicode(u.decode('utf-8')) \
                if not unicodedata.category(c).startswith('C')]).\
                    strip().replace(u'\xa0', '') \
            for u in ll[3].split('/')]])) \
        for ll in \
        [l.split('\t') for l in f.read().split('\n') if len(l) > 0]])

names = uniqList(map(lambda x: x[0], proteins))
uniprots = uniqList(flatList(names_dict.values()))

# missing from the network
# only one rat protein and Claudin-6
# and this is weird:
# Q13972 (RASGRF1), P07101 (Tyrosine 3-monooxygenase): 
# rat & mouse only, but human UniProt ID
set(uniprots) - set(net.graph.vs['name'])

# processing PTMs:
ptms_proteins = []
for p in proteins:
    ptmts = [] if len(p[1]) == 0 else remod.findall(p[1])
    if p[0] not in names_dict:
        print '\t#### Missing from ID file: %s' % p[0]
        continue
    uni = names_dict[p[0]]
    for ptmt in ptmts:
        residue_match = False
        mod = ptmtyp[p[1].split('-')[0]]
        aa = aarev[ptmt[0].upper()]
        num = int(ptmt[1])
        for u in uni:
            if u in sprot_seq:
                for plus in [0, 1, -1, 2, -2]:
                    for iso, sq in sprot_seq[u].isof.iteritems():
                        aau = sprot_seq[u].get(num + plus, isoform = iso)
                        if aau == aa:
                            residue_match = True
                            res = intera.Residue(num + plus, aa, u, isoform = iso)
                            ptm = intera.Ptm(u, residue = res, typ = mod, 
                                source = 'assay', isoform = iso)
                            values = None if p[2] not in paired else paired[p[2]]
                            ptms_proteins.append((u, ptm, values))
                    if residue_match:
                        break
        if not residue_match:
            print 'Residue does not match: %s, %u, %s, %s'%(u, num, aau, aa)
            print '\t', p[0], p[1]

# this was not the best way, although found almost all sites.
# maybe better to lookup all the sites in PhosphoSite

# processing PTMs with PhosphoSite:
ptms_proteins = []
with open(ptmlog, 'w') as logf:
    for p in proteins:
        ptmts = [] if len(p[1]) == 0 else remod.findall(p[1])
        if p[0] not in names_dict:
            msg = '<!!!> Missing from ID file: %s\n' % p[0]
            logf.write(msg)
            sys.stdout.write(msg)
            sys.stdout.flush()
            continue
        uni = names_dict[p[0]]
        if len(ptmts) > 0:
            for u in uni:
                residue_match = False
                if u in sprot_seq:
                    pdiff = {}
                    for ptmt in ptmts:
                        mod = ptmtyp[p[1].split('-')[0]]
                        aa = aarev[ptmt[0].upper()]
                        num = int(ptmt[1])
                        # Collecting all PTMs of all isoforms for this UniProt:
                        for iso, sq in sprot_seq[u].isof.iteritems():
                            if u in psitep:
                                for psptm in psitep[u]:
                                    psres = psptm.residue.number
                                    psaa = psptm.residue.name
                                    diff = psres - num
                                    uaa = sprot_seq[u].get(psres, isoform = iso)
                                    incons = uaa != psaa and psptm.isoform == iso
                                    if incons:
                                        msg = 'Inconsistency between '\
                                            'UniProt and PhosphoSite:\n'\
                                            '\t%s-%u:%u UP: %s PS: %s\n' % \
                                            (u, iso, psres, uaa, psaa)
                                        logf.write(msg)
                                        sys.stdout.write(msg)
                                        sys.stdout.flush()
                                    if abs(diff) not in pdiff:
                                        pdiff[abs(diff)] = []
                                    pdiff[abs(diff)].append((diff, uaa, incons, psptm))
                    for ptmt in ptmts:
                        mod = ptmtyp[p[1].split('-')[0]]
                        aa = aarev[ptmt[0].upper()]
                        num = int(ptmt[1])
                        # Checking first all isoforms assuming original PTM position:
                        for iso, sq in sprot_seq[u].isof.iteritems():
                            aau = sprot_seq[u].get(num, isoform = iso)
                            if aau == aa:
                                residue_match = True
                                srcs = ['assay']
                                # let's see if this PTM is in PhosphoSite:
                                if 0 in pdiff:
                                    for psptm in pdiff[0]:
                                        if not psptm[2] and psptm[3].isoform == iso \
                                            and psptm[1] == aa:
                                            srcs.append('PhosphoSite')
                                            break
                                res = intera.Residue(num, aa, u, isoform = iso)
                                ptm = intera.Ptm(u, residue = res, typ = mod, 
                                    source = srcs, isoform = iso)
                                values = None if p[2] not in paired else paired[p[2]]
                                ptms_proteins.append((u, ptm, values))
                        # OK, if there is no match 
                    if not residue_match and len(pdiff) > 0:
                        mindiff = min(pdiff.keys())
                        if mindiff <= 5:
                            for psptm in pdiff[mindiff]:
                                if not psptm[2]:
                                    for ptmt in ptmts:
                                        mod = ptmtyp[p[1].split('-')[0]]
                                        aa = aarev[ptmt[0].upper()]
                                        num = int(ptmt[1])
                                        if psptm[1] == aa:
                                            res = intera.Residue(psptm[3].residue.number, 
                                                aa, u, isoform = psptm[3].isoform)
                                            ptm = intera.Ptm(u, residue = res, typ = mod, 
                                                source = ['PhosphoSite'], 
                                                isoform = psptm[3].isoform)
                                            values = None if p[2] not in paired \
                                                else paired[p[2]]
                                            ptms_proteins.append((u, ptm, values))
                                            msg = 'Residue does not match:\n'\
                                                '\t%s:%s%u\n'%(u, aa, num)
                                            msg += '\tSubstituting with closest site'\
                                                ' from PhosphoSite (at offset %s):\n'\
                                                '\t%s-%u:%s%u\n'%(str(psptm[0]), u, 
                                                psptm[3].isoform, aa,
                                                psptm[3].residue.number)
                                            if mindiff > 2:
                                                msg += '\tWARNING: difference in '\
                                                    'offset larger than 2: %s\n'%\
                                                    (str(psptm[0]),)
                                            logf.write(msg)
                                            sys.stdout.write(msg)
                                            sys.stdout.flush()
                                            residue_match = True
                    if not residue_match:
                        msg = '<!!!> No match found in PhosphoSite within range +-5:\n'
                        msg += '\t%s:%s\n'%(u, \
                            ':'.join(['%s%u'%(aarev[ptmt[0].upper()], int(ptmt[1])) \
                            for ptmt in ptmts]))
                        logf.write(msg)
                        sys.stdout.write(msg)
                        sys.stdout.flush()
                        # if not found in PhosphoSite, take a look at the range +-2:
                        for plus in [1, -1, 2, -2]:
                            for iso, sq in sprot_seq[u].isof.iteritems():
                                aau = sprot_seq[u].get(num + plus, isoform = iso)
                                if aau == aa:
                                    residue_match = True
                                    res = intera.Residue(num + plus, aa, u, isoform = iso)
                                    ptm = intera.Ptm(u, residue = res, typ = mod, 
                                        source = 'assay', isoform = iso)
                                    values = None if p[2] not in paired else paired[p[2]]
                                    ptms_proteins.append((u, ptm, values))
                                    msg = '\tUsing residue at offset %s, '\
                                        'although this site '\
                                        'is not in PhosphoSite:\n'%str(plus)
                                    msg += '\t%s-%u:%s%u\n'%(u, iso, aa, num + plus)
                                    logf.write(msg)
                                    sys.stdout.write(msg)
                                    sys.stdout.flush()
                            if residue_match:
                                break
                        if not residue_match:
                            msg = '<!!!> No match found within range +-2:\n'
                            msg += '\t%s:%s\n'%(u, \
                                ':'.join(['%s%u'%(aarev[ptmt[0].upper()], int(ptmt[1])) \
                                for ptmt in ptmts]))
                            logf.write(msg)
                            sys.stdout.write(msg)
                            sys.stdout.flush()

# OMG, this was a mess!

# creating network attributes:
net.graph.vs['sites'] = [[] for _ in net.graph.vs]
for site in ptms_proteins:
    if site[0] in net.graph.vs['name']:
        net.graph.vs[net.graph.vs['name'].index(site[0])]['sites'].append(site)

net.load_ptms()
haskinase = []
nokinase = []
for v in net.graph.vs:
    for ptm in v['sites']:
        this_site_has_kinase = False
        for k in v.neighbors():
            e = net.graph.get_eid(v.index, k.index)
            for kinsub in net.graph.es[e]['ptm']:
                if kinsub.ptm == ptm[1]:
                    this_site_has_kinase = True
        if this_site_has_kinase:
            haskinase.append((v['name'], ptm[1].residue.number))
        else:
            nokinase.append((v['name'], ptm[1].residue.number))

net.graph.vs['assay'] = [v['name'] in uniprots for v in net.graph.vs]
sum(net.graph.vs['assay'])
len([v for v in net.graph.vs if v['assay'] and 'Signor' in v['sources']])

set(uniqList(haskinase)) & set(uniqList(nokinase))
len(uniqList(haskinase))
len(set(uniqList(nokinase)) - set(uniqList(haskinase)))
len([(p[0], p[1].residue.number) for p in  ptms_proteins])
