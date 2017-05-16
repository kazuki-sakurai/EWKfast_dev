#! /usr/bin/env python

import os, sys
import numpy as np
from numpy import linalg as la
from numpy import random as rn
import iminuit as im
from math import *
import subprocess as sbproc
from collections import OrderedDict

try:
    order=sys.argv[1]
    rs=sys.argv[2]
    pdf=sys.argv[3]    
except:
    print '[order] [rs] [pdf]'
    exit()

data = OrderedDict()
i = 0
for line in open('validation_points.dat'):
    i += 1
    slhaname, i1, i2, mode = line.split()
    tagid = int(slhaname.split('-')[1].split('.dat')[0])
    if tagid not in data.keys(): data[tagid] = []
    data[tagid].append([mode, i1, i2, i])

id_list = sorted(data.keys())

header = '''
ScaleVariation: OFF
Unit: pb
'''

for tagid in id_list:
    f = open('input_{order}{rs}/input_{id}.dat'.format(id=tagid, order=order, rs=rs), 'w')
    f.write(header)
    for elems in data[tagid]:
        mode, i1, i2, iline = elems
        f.write('{rs}  {order}  {pdf}  {mode} \n'.format(order=order, rs=rs, pdf=pdf, mode=mode))
