#!/usr/bin/python
# -*- coding: utf-8 -*-

import bmp8
import os

b = bmp8.Bmp8(ncbi_tax_id = 9606)
b.init()
b.export_table(to_file = True)
b.export_fc_table(to_file = True)
b.kinact_analysis()
