#! /usr/bin/env python

import os, sys
import numpy as np
from math import *
from scipy import interpolate
from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator
from collections import OrderedDict
from re import findall 

pref = 'LO-13'
#pdfs = ['NNPDF30']
pdfs = ['CT14lo']

#grids = ['C2+C1-', 'C2-C1+']
#grids = ['NN']
#grids = ['NC-']
#grids = ['CCsame', 'NNsame']
#grids = ['NNsame']
#grids = ['CCsame']
#grids = ['C2-C1+','C2+C1-']
#grids = ['NC+','NC-']
grids = ['NN']

for grid in grids:
    for pdf in pdfs:

        print grid, pdf

        if 'same' in grid:
            dim = 2
        else:
            dim = 3

        outdir = '{pdf}/lookups'.format(pdf=pdf)
        tag = '{pdf}/{pref}_{grid}_{pdf}'.format(pref=pref, grid=grid, pdf=pdf)        
        filename = '{tag}.table'.format(tag=tag)

        if dim == 2:
            dir_path = os.path.dirname(os.path.realpath(__file__))
            fpath = os.path.join(dir_path, filename)
            data = np.loadtxt(fpath)
            m1_ar = data[:,0]
            m2_ar = data[:,1]
            m1 = np.array(sorted(list(set(m1_ar))))
            m2 = np.array(sorted(list(set(m2_ar))))            
            m1.dump('{outdir}/{tag}.m1'.format(outdir=outdir, tag=tag))
            m2.dump('{outdir}/{tag}.m2'.format(outdir=outdir, tag=tag))
            F_ar = []
            sizeF = np.shape(data)[1] - dim           
            for ii in xrange(sizeF): 
                Fdm = data[:,ii+dim]
                Far = Fdm.reshape(len(m1), len(m2))
                Far.dump('{outdir}/{tag}.F{i}'.format(outdir=outdir, tag=tag, i=ii+1))

        if dim == 3:
            dir_path = os.path.dirname(os.path.realpath(__file__))
            fpath = os.path.join(dir_path, filename)
            data = np.loadtxt(fpath)
            m1_ar = data[:,0]
            m2_ar = data[:,1]
            m3_ar = data[:,2]            
            m1 = np.array(map(int, sorted(list(set(m1_ar)))))
            m2 = np.array(map(int, sorted(list(set(m2_ar)))))
            m3 = np.array(map(int, sorted(list(set(m3_ar)))))              
            m1.dump('{outdir}/{tag}.m1'.format(outdir=outdir, tag=tag))
            m2.dump('{outdir}/{tag}.m2'.format(outdir=outdir, tag=tag))            
            m3.dump('{outdir}/{tag}.m3'.format(outdir=outdir, tag=tag))
            sizeF = np.shape(data)[1] - dim
            Fs = []
            for ii in xrange(sizeF):
                Fdic = {}
                for line in open(fpath):
                    elems = line.split()
                    d1, d2, d3 = int(elems[0]), int(elems[1]), int(elems[2])
                    Fdic[(d1,d2,d3)] = float(elems[ii+dim])
                Fs.append(Fdic)

            for ii in xrange(sizeF): 
                Far = []
                for x1 in m1:
                    Fdm2 = []
                    for x2 in m2:
                        Fdm1 = []
                        for x3 in m3:
                            key = (x1,x2,x3)
                            Fdm1.append(Fs[ii][key])
                        Fdm2.append(Fdm1)
                    Far.append(Fdm2)                
                Far = np.array(Far)
                print np.shape(Far)
                Far.dump('{outdir}/{tag}.F{i}'.format(outdir=outdir, tag=tag, i=ii+1))

        # if dim == 3 and False:
        #     dir_path = os.path.dirname(os.path.realpath(__file__))
        #     fpath = os.path.join(dir_path, filename)
        #     data = np.loadtxt(fpath)
        #     m1_ar = data[:,0]
        #     m2_ar = data[:,1]
        #     m3_ar = data[:,2]            
        #     m1 = np.array(sorted(list(set(m1_ar))))
        #     m2 = np.array(sorted(list(set(m2_ar))))            
        #     m3 = np.array(sorted(list(set(m3_ar))))                        
        #     m1.dump('{outdir}/{tag}.m1'.format(outdir=outdir, tag=tag))
        #     m2.dump('{outdir}/{tag}.m2'.format(outdir=outdir, tag=tag))            
        #     m3.dump('{outdir}/{tag}.m3'.format(outdir=outdir, tag=tag))                        
        #     F_ar = []
        #     sizeF = np.shape(data)[1] - dim - 1            
        #     # print filename
        #     # print len(m1), len(m2), len(m3), np.shape(data[:,1+dim])
        #     for ii in xrange(sizeF): 
        #         Fdm = data[:,ii+dim]
        #         Far = Fdm.reshape(len(m1), len(m2), len(m3))
        #         F_ar.append( RegularGridInterpolator( (m1, m2, m3), Far) )
        #         Far.dump('{outdir}/{tag}.F{i}'.format(outdir=outdir, tag=tag, i=ii+1))

