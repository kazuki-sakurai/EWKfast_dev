#!/usr/bin/env python

import sys, os
import numpy as np

def get_neighbour(i1, m2, m3, m1ar, Fdic, key_orig):    
    for i in xrange(i1, len(m1ar) ):
        m1 = m1ar[i]
        key = (m1,m2,m3)
        if key in key_orig:
            Fs = Fdic[key]
            return Fs

def get_line(List):
    line = str(List[0])
    for i in range(1, len(List)):
        line += '  ' + str(List[i])
    return line

infiles = sys.argv[1:]

pref = 'LO-13'
pdf = 'CT14lo'
#pdf = 'NNPDF30'

#modes = ['NN', 'NC-', 'NC+', 'C2-C1+', 'C2+C1-']
#modes = ['NN', 'C2-C1+', 'C2+C1-']
#modes = ['C2-C1+', 'C2+C1-',]
modes = ['NN']
#modes = ['NC-', 'NC+']

for mode in modes:

    if mode in ['C2-C1+', 'C2+C1-']:
        signs = ['_pp']
    else:
        signs = ['_pp', '_pm']

    outfile = '{pdf}/{pref}_{mode}_{pdf}.table'.format(pref=pref, mode=mode, pdf=pdf)

    m1ar, m2ar, m3ar = [], [], []
    Fdic = {}
    for sign in signs:
        infile = '{pdf}/pp_pm/{pref}_{mode}{sign}_{pdf}.table'.format(pref=pref, mode=mode, sign=sign, pdf=pdf)
        for line in open(infile):
            elems = line.split()
            m1, m2, m3 = int(elems[0]), int(elems[1]), int(elems[2])
            # Fs = elems[3:-1]
            # chi2 = elems[-1]
            Fs = elems[3:]
            key = (m1, m2, m3)
            m1ar.append(m1)
            m2ar.append(m2)
            m3ar.append(m3)
            Fdic[key] = Fs
    m1ar = sorted(list(set(m1ar)))
    m2ar = sorted(list(set(m2ar)))
    m3ar = sorted(list(set(m3ar)))
    key_orig = Fdic.keys()

    for m3 in m3ar:
        for m2 in m2ar:
            for i1 in xrange(len(m1ar)):
                m1 = m1ar[i1]
                key = (m1, m2, m3)
                #print Fdic.keys()
                if key not in key_orig:                 
                    Fs = get_neighbour(i1, m2, m3, m1ar, Fdic, key_orig) 
                    Fdic[key] = Fs

    f = open(outfile,'w')
    for m3 in m3ar:
        for m2 in m2ar:
            for m1 in m1ar:                    
                key = (m1, m2, m3)
                line = '{m1}  {m2}  {m3}   '.format(m1=m1, m2=m2, m3=m3)
                line += get_line(Fdic[key]) 
                f.write(line + '\n')

    print outfile