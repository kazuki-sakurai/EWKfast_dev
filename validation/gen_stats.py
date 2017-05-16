#!/usr/bin/env python

import sys, os
from math import *
from subprocess import Popen, PIPE, call

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cls

order = 'NLO'
rs = '13'
# pdf='cteq6'
# pdf='CT14'
# pdf = 'MMHT14'
pdf = 'NNPDF'

#method = 'nearest'
method = ''

thres, tag, extra = 10**-4, '_01fb', r'$\,\,\,\, \sigma > 0.1\,\rm fb$'
#thres, tag, extra = 10**-10, '', ''

xs_d, xsp_d, unc_d = {},{},{}

min_runc = 1000000.

if method == 'nearest':
    compfile = 'comparison/comparison_{order}{rs}_{pdf}_nearest.dat'.format(order=order, rs=rs, pdf=pdf)
else:
    compfile = 'comparison/comparison_{order}{rs}_{pdf}.dat'.format(order=order, rs=rs, pdf=pdf)

for line in open(compfile):
    elems = line.split()
    i, mode = elems[:2]    
    m1, m2, mu, ma, mb, mQ, xs, xsp = map(float, elems[2:])

    #if min(xsp, xs) < 0.5*10**-5: continue
    if min(xsp, xs) < thres: continue

    if mode not in xs_d.keys():
        xs_d[mode], xsp_d[mode], unc_d[mode] = [],[],[]
    xs_d[mode].append(xs)
    xsp_d[mode].append(xsp)
    #unc_d[mode].append(unc)

    #if unc/xsp < min_runc:
    #    min_runc = unc/xsp
    #    print mode, unc/xsp

modes1 = xs_d.keys()

cols = ['r']*2 + ['orange']*2 + ['g']*4 + ['b']*6 + ['k']*16
#cols = ['r']*2 + ['orange']*2 + ['g']*4 + ['b']*6

modes = ['C1C1', 'C2C2']
modes += ['C1+C2-', 'C1-C2+']
modes += ['N1N1', 'N2N2', 'N3N3', 'N4N4']
modes += ['N1N2', 'N1N3', 'N1N4']
modes += ['N2N3', 'N2N4']
modes += ['N3N4']
#cols = ['orange']*2 + ['g']*4 + ['b']*6 + ['k']*16
modes += ['N1C1+', 'N1C2+', 'N2C1+', 'N2C2+']
modes += ['N3C1+', 'N3C2+', 'N4C1+', 'N4C2+']
modes += ['N1C1-', 'N1C2-', 'N2C1-', 'N2C2-']
modes += ['N3C1-', 'N3C2-', 'N4C1-', 'N4C2-']


for mode in modes:
    xs_d[mode] = np.array(xs_d[mode])
    xsp_d[mode] = np.array(xsp_d[mode])
    unc_d[mode] = np.array(unc_d[mode])

#modes = ['N1N1']
#for val in ['ratio', 'scale']:
for val in ['ratio']:

    v_list, v_mean, v_std = [], [], []

    for mode in modes:

        ratio = 0.
        scale = 0.

        xs = xs_d[mode] 
        xsp = xsp_d[mode] 
        #unc = unc_d[mode]
        if val == 'ratio': v = (xs - xsp)/np.maximum(xsp, xs)
        #if val == 'scale': v = (xs - xsp)/unc

        v_list.append(v)
        v_mean.append(v.mean())
        v_std.append(v.std())

        #print val, mode, '  ', np.mean(v), np.std(v)



    fig = plt.figure()
    ax = fig.add_subplot(111) 

    x_list = np.arange(len(v_mean))
    # ax.errorbar(x_list, v_mean, v_std, color=['red']*30, fmt='ok', lw=1.5)

    for i in xrange(len(v_mean)):
        s = v_std[i]
        m = v_mean[i]
        c = cols[i]
        ax.errorbar(i, m, s, color=c, fmt='ok', lw=1.5)

    xmin, xmax = -1, len(v_mean)
    ax.plot([xmin, xmax], [0, 0], lw=1, ls='-', color='gray')
    if val == 'ratio':
        for y in [0.1, 0.2, 0.3, 0.4]:
            ax.plot([xmin, xmax], [y, y],   lw=0.3, ls='-', color='gray')
            ax.plot([xmin, xmax], [-y, -y], lw=0.3, ls='-', color='gray')


    ax.set_xlim([xmin, xmax])

    if val == 'ratio': ymin, ymax = -0.5, 0.5
    if val == 'scale': 
        ymin, ymax = -50, 50
        ymin, ymax = -100, 100
        if tag != '':
            ymin, ymax = -20, 20
    if [val, order] == ['scale', 'LO']:  
        ymin, ymax = -20, 20

    ax.set_ylim([ymin, ymax])
    plt.xticks(x_list, modes, rotation=90, fontsize=9)

    if val == 'ratio': ax.set_title('Means and STDs of (EWKF - Resum)/Resum @{order} {rs}TeV {pdf} {method}'.format(order=order, rs=rs, pdf=pdf, method=method) + extra, fontsize=11)
    if val == 'scale': ax.set_title('Means and STDs of (EWKF - Resum)/(Scale Variation) @{order} {rs}TeV {pdf} {method}'.format(order=order, rs=rs, pdf=pdf, method=method) + extra, fontsize=11)

    #ax.boxplot(r_list)
    #ax.violinplot(r_list)

    # outdir = 'plots/{zval}'.format(zval=zval)
    # if not os.path.exists(outdir): call(['mkdir', outdir]) 
    # figname = os.path.join(outdir, '{mode}.pdf'.format(mode=mode))

    #if not os.path.exists('plots/' + mode): call(['mkdir', 'plots/' + mode]) 
    #figname = 'plots/{mode}/{zval}_{cut}_{p_mode}.pdf'.format(cut=cut, mode=mode, p_mode=p_mode, zval=zval)
    #fig.savefig(figname)
    #print figname
    if method == 'nearest':
        figname = 'plots/{order}{rs}_nearest/statistics_{val}{tag}_{pdf}.pdf'.format(val=val, order=order, rs=rs, tag=tag, pdf=pdf)
    else:
        figname = 'plots/{order}{rs}/statistics_{val}{tag}_{pdf}.pdf'.format(val=val, order=order, rs=rs, tag=tag, pdf=pdf)
    fig.savefig(figname)
    print figname


exit()
