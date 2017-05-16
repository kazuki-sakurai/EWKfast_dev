#! /usr/bin/env python

import os, sys
import numpy as np
import json
from math import *
import subprocess as sbproc
from collections import OrderedDict
from functions import *
from StringIO import StringIO

def get_mass(i, mode):
    infile = 'slha_files/slha-{i}.dat'.format(i=i)
    p = get_params(infile)
    if mode.count('C') == 2:
        i1 = int(mode.split('C')[1][:1])
        i2 = int(mode.split('C')[2][:1])        
        ma, mb = p['mC'][int(i1)-1], p['mC'][int(i2)-1]
    if mode.count('N') == 2:
        dm, i1, i2 = mode.split('N')
        ma, mb = p['mN'][int(i1)-1], p['mN'][int(i2)-1]
    if mode.count('N') == 1 and mode.count('C') == 1:
        nn = int(mode.split('N')[1][:1])
        cc = int(mode.split('N')[1].split('C')[1][:1])
        ma, mb = p['mN'][int(nn)-1], p['mC'][int(cc)-1]

    mQ = p['mQ']

    for line in open(infile):
        elems = line.split()
        if len(elems) > 3:
            if 'M_1(Q)' in elems[3] and elems[0] == '1': 
                M1 = float(line.split()[1])
            if 'M_2(Q)' in elems[3] and elems[0] == '2': 
                M2 = float(line.split()[1])
            if 'mu(Q)MSSM' in elems[3] and elems[0] == '1': 
                mu = float(line.split()[1])

    return M1, M2, mu, ma, mb, mQ

rs = '13'
order = 'NLO'
#pdf='cteq6'
#pdf='CT14'
#pdf='MMHT14'
pdf='NNPDF'

#method = 'nearest'
method = ''

if method == 'nearest': 
    fout = open('comparison/comparison_{order}{rs}_{pdf}_nearest.dat'.format(order=order, rs=rs, pdf=pdf), 'w')
else:
    fout = open('comparison/comparison_{order}{rs}_{pdf}.dat'.format(order=order, rs=rs, pdf=pdf), 'w')

for i in xrange(1, 601):
    if method == 'nearest': 
        infile = 'EWKFresult_{order}{rs}_nearest/result-{i}.dat'.format(i=i, order=order, rs=rs)    
    else:
        infile = 'EWKFresult_{order}{rs}/result-{i}.dat'.format(i=i, order=order, rs=rs)        
    for line in open(infile):
        elems = line.split()
        xs = ''
        if len(elems) == 5:  # version 0
            rs2, order2, pdfdm, mode, xs = elems                            
            if rs != rs2: continue
            if order2 != order:
                print 'ERROR 1'
                exit()
        if len(elems) == 0: 
            break
        if elems[0] == 'Energy': continue                            
        #print infile
        #print elems
        xs = float(xs)
        if xs < 10**-6: continue
        M1, M2, mu, ma, mb, mQ = get_mass(i, mode)

        ####
        resum_file = 'data_resum/E{rs}_{pdf}/point-{i}/{mode}.json'.format(i=i, rs=rs, pdf=pdf, mode=mode)
        if not os.path.exists(resum_file):
            print resum_file, 'does not exist.'
            continue 
        f = open(resum_file).read()
        f = f.replace('\n', ' ')
        f = f.replace('.,', ',')
        data = json.loads(f)
        if order == 'LO': xsv = float(data['lo'])
        if order == 'NLO': xsv = float(data['nlo'])        
        # try:
        #     f = open(pros_file)        
        #     xsp05  = float( f.readline().split()[column])
        #     elems = f.readline().split()
        #     xsp1     = float( elems[column]  )
        #     xsp_dm  = float( elems[column2] )            
        #     xsp2  = float( f.readline().split()[column])
        #     unc1 = abs( xsp1 - xsp05 )
        #     unc2 = abs( xsp2 - xsp1 )            
        #     sc_unc = (unc1 + unc2)/2.
        # except:
        #     print 'error in ', pros_file
        #     continue

        outline = '{i}  {mode}  {M1} {M2} {mu}  {ma} {mb} {mQ}  {xs}  {xsv}'.format(i=i, mode=mode, M1=M1, M2=M2, mu=mu, ma=ma, mb=mb, mQ=mQ, xs=xs, xsv=xsv)
        #print outline
        fout.write(outline + '\n')



