#!/usr/bin/python
# -*- coding: utf-8 -*-

import bmp8
import os

b = bmp8.Bmp8(ncbi_tax_id = 9606)
b.init()
b.export_table(to_file = True)
b.export_fc_table(to_file = True)
b.fc_top_table()
b.kinact_analysis()
b.kinact_top(threshold = 1.0)
b.kinact_top(threshold = 1.0, fname = 'kinact_top.tab')

b.top_fc_venn(fname = 'FC_1.5_both.pdf', threshold = 1.5, title = '> 1.5 :: up & down')
b.top_fc_venn(fname = 'FC_1.5_down.pdf', sign = lambda s: s < 0.0, threshold = 1.5, title = '> 1.5 :: down')
b.top_fc_venn(fname = 'FC_1.5_up.pdf', sign = lambda s: s > 0.0, threshold = 1.5, title = '> 1.5 :: up')

b.top_fc_venn(fname = 'FC_2.0_both.pdf', threshold = 2.0, title = '> 2.0 :: up & down')
b.top_fc_venn(fname = 'FC_2.0_down.pdf', sign = lambda s: s < 0.0, threshold = 2.0, title = '> 2.0 :: down')
b.top_fc_venn(fname = 'FC_2.0_up.pdf', sign = lambda s: s > 0.0, threshold = 1.678, title = '> 1.678 :: up')
