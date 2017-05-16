#!/usr/bin/env python

import sys, os, json
import numpy as np
from func_coeff import * 

def get_line(List):
    line = str(List[0])
    for i in range(1, len(List)):
        line += '  ' + str(List[i])
    return line

order = 'LO'
rs = '13'
OE = '{order}-{rs}'.format(order=order, rs=rs)
OEdm = 'NLO-13'
#modes = ['NC+_pp', 'NC-_pp', 'NC+_pm', 'NC-_pm']
#modes = ['CCsame', 'NNsame']
modes = ['NN_pp', 'NN_pm']
#modes = ['C2-C1+', 'C2+C1-']

for mode in modes:

    if order == 'LO': order_tag = 'lo'
    if order == 'NLO': order_tag = 'nlo'

    if mode == 'NNsame':    nsamp, i1, i2, grid_path = 15, 0, 0, '../grids/grid_same.dat'
    if mode == 'CCsame':    nsamp, i1, i2, grid_path = 20, 0, 0, '../grids/grid_same.dat'
    if mode == 'C2+C1-':    nsamp, i1, i2, grid_path = 10, 0, 1, '../grids/grid_XX++.dat'
    if mode == 'C2-C1+':    nsamp, i1, i2, grid_path = 10, 0, 1, '../grids/grid_XX++.dat'
    if mode == 'NN_pp':     nsamp, i1, i2, grid_path = 15, 0, 1, '../grids/grid_XX++.dat'
    if mode == 'NN_pm':     nsamp, i1, i2, grid_path = 15, 0, 1, '../grids/grid_XX+-.dat'
    if mode == 'NC+_pp':    nsamp, i1, i2, grid_path = 20, 0, 0, '../grids/backup_3/grid_NC++.dat'
    if mode == 'NC+_pm':    nsamp, i1, i2, grid_path = 20, 0, 0, '../grids/backup_3/grid_NC+-.dat'
    if mode == 'NC-_pp':    nsamp, i1, i2, grid_path = 20, 0, 0, '../grids/backup_3/grid_NC++.dat'
    if mode == 'NC-_pm':    nsamp, i1, i2, grid_path = 20, 0, 0, '../grids/backup_3/grid_NC+-.dat'

    if 'same' in grid_path: 
        m1, m2 = np.loadtxt(grid_path).transpose()
    else: 
        m1, m2, m3 = np.loadtxt(grid_path).transpose()
    ptot = len(m1)

    pdf = ''

    Coeffs = []
    Chisq = []
    for ip in xrange(1, ptot+1):

        xsec = []
        vecs = []
        for isamp in xrange(1, nsamp+1):
            dirname = 'data/{OEdm}_{mode}/isamp-{isamp}/point-{ip}'.format(OEdm=OEdm, mode=mode, isamp=isamp, ip=ip)
            jsonfile = dirname + '/result.json'
            f = open(jsonfile).read()
            f = f.replace('\n', ' ')
            f = f.replace('.,', ',')
            data = json.loads(f)
            xsec.append(data[order_tag])
            if order == 'LO' and pdf =='': pdf = data['pdflo'].split('_')[0]
            if order == 'NLO' and pdf =='': pdf = data['pdfnlo'].split('_')[0]

            slhafile = dirname + '/slha.in'
            #print slhafile
            #blocks, decays = readSLHAFile(slhafile, ignorenobr=True)
            #N, V, U = get_mixings(blocks)
            N, V, U = get_mixings(slhafile)

            vec = get_vector(mode, N, V, U, i1, i2)
            vecs.append(vec)


        Fini = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        chi2_result, coeff = get_coeff(mode, vecs, xsec, Fini)

        #print chi2_result
        Coeffs.append(coeff)
        Chisq.append(chi2_result[0])

    stag = ''
    if 'pp' in mode: stag = '/pp_pm'
    if 'pm' in mode: stag = '/pp_pm'
    fout = open('tables/{pdf}{stag}/{OE}_{mode}_{pdf}.table'.format(OE=OE, mode=mode, pdf=pdf, stag=stag), 'w')

    for i in xrange(len(m1)):
        if 'same' in mode:
            line = '{m1}  {m2}   '.format(m1=int(m1[i]), m2=int(m2[i])) + get_line(Coeffs[i]) + '\n'
        else:
            line = '{m1}  {m2}  {m3}   '.format(m1=int(m1[i]), m2=int(m2[i]), m3=int(m3[i])) + get_line(Coeffs[i]) + '\n'
        fout.write(line)

    print mode, max(Chisq), min(Chisq), np.array(Chisq).mean()

