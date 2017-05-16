#! /usr/bin/env python

import os, sys
import numpy as np
from numpy import linalg as la
from numpy import random as rn
import iminuit as im
from math import *
import subprocess as sbproc

FNULL = open(os.devnull, 'w')

sw = sqrt(1. - 80.410003662109375**2 / 91.186996459960938**2)
cw = sqrt(1 - sw**2)
r2 = sqrt(2.)

deno = sw * cw
Lu = ( 1./2. - 2./3. * sw**2 ) / deno
Ld = (-1./2 - (-1./3.) * sw**2 ) / deno
Ru = (-2./3. * sw**2 ) / deno
Rd = ( - (-1./3.) * sw**2 ) / deno

def get_line(List):
    line = str(List[0])
    for i in range(1, len(List)):
        line += '  ' + str(List[i])
    return line

def Lut(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (2./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Lu

def Ldt(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (-1./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Ld

def Rut(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (2./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Ru

def Rdt(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (-1./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Rd

def gNC( i, j, n, v, u ):
    elem1 = n[i, 1] * v[j, 0] - n[i, 3] * v[j, 1]/sqrt(2.) 
    elem2 = v[j, 0] * Lut(i, n)
    elem3 = n[i, 1] * u[j, 0] + n[i, 2] * u[j, 1]/sqrt(2.) 
    elem4 = u[j, 0] * Ldt(i, n)
    return [elem1, elem2, elem3, elem4]

def cNC( v, c ):
    if c == 0: return v[0]**2 + v[2]**2
    if c == 1: return v[1]**2 + v[3]**2
    if c == 2: return 2.*( v[0]*v[1] - v[2]*v[3] )
    if c == 3: return 2.*( v[1]*v[2] - v[0]*v[3] )
    if c == 4: return 2.* v[0]*v[2]
    if c == 5: return 2.* v[1]*v[3]

def gNN( i, j, n ):
    elem1 = n[i, 3] * n[j, 3] - n[i, 2] * n[j, 2]
    elem2 = Lut(i, n) * Lut(j, n)
    elem3 = Rut(i, n) * Rut(j, n)
    elem4 = Ldt(i, n) * Ldt(j, n)
    elem5 = Rdt(i, n) * Rdt(j, n)
    return [elem1, elem2, elem3, elem4, elem5]

def cNN( v, c ):
    if c == 0: return v[0]**2
    if c == 1: return v[1]**2 + v[2]**2
    if c == 2: return 2.*( Lu * v[0] * v[1] - Ru * v[0] * v[2] )
    if c == 3: return v[3]**2 + v[4]**2
    if c == 4: return 2.*( Ld * v[0] * v[3] - Rd * v[0] * v[4] )    


def gCC( i, j, u, v ):
    elem1 = u[i, 0] * u[j, 0]
    elem2 = v[i, 0] * v[j, 0]
    elem3 = 0
    if i == j: elem3 = 1
    return [elem1, elem2, elem3]

def cCC( v, c ):
    if c == 0: return v[0]**2
    if c == 1: return 2. * v[0] * v[1]
    if c == 2: return v[1]**2
    if c == 3: return 2. * v[0] * v[2]
    if c == 4: return 2. * v[1] + v[2]
    if c == 5: return v[2]**2

def get_vector(mode, N, V, U, i1, i2):

    if 'NC' in mode:
        gnc = gNC(i1, i2, N, V, U)
        vec = [ cNC(gnc, i) for i in xrange(6) ]
    if 'NN' in mode:
        gnn = gNN(i1, i2, N)
        vec = [ cNN(gnn, i) for i in xrange(5) ]
    if mode in ['CCsame', 'C2+C1-', 'C2-C1+']:
        gcc = gCC(i1, i2, U, V)
        n = 3
        if mode == 'CCsame': n=6
        vec = [ cCC(gcc, i) for i in xrange(n) ]
    return vec

######################################
######################################
######################################

def get_coeff(run_mode, vecs, xsec, Fini):

    #xsec_repro = []
    chi2_result = []
    #Coeff = []

    ######################################
    if 'NC' in run_mode:

        def chi2(F1, F2, F3, F4, F5, F6):
            __coeff = np.array([F1, F2, F3, F4, F5, F6])
            MR = (np.dot( vecs, __coeff ) - xsec) / xsec
            ML = np.transpose(MR)
            return np.dot(ML, MR)

        m = im.Minuit(chi2, 
            F1=Fini[0], F2=Fini[1], F3=Fini[2], F4=Fini[3], F5=Fini[4], F6=Fini[5],            
            #F1=0, F2=0, F3=0, F4=0, F5=0, F6=0,
            #F1=1, F2=2, F3=3, F4=4, F5=5, F6=6,
            #F1=6, F2=5, F3=4, F4=3, F5=2, F6=1,
            error_F1=0.001, 
            error_F2=0.001, 
            error_F3=0.001, 
            error_F4=0.001, 
            error_F5=0.001, 
            error_F6=0.001, 
            print_level=1) 

        m.migrad()
        vals = m.values

        F1, F2, F3, F4, F5, F6 = vals['F1'], vals['F2'], vals['F3'], vals['F4'], vals['F5'], vals['F6']
        chi2_result.append(chi2(F1, F2, F3, F4, F5, F6)) 
        __coeff = np.array([F1, F2, F3, F4, F5, F6])

    ######################################
    if 'NN' in run_mode:

        def chi2(F1, F2, F3, F4, F5):
            __coeff = np.array([F1, F2, F3, F4, F5])
            MR = (np.dot( vecs, __coeff ) - xsec) / xsec
            ML = np.transpose(MR)
            return np.dot(ML, MR)

        m = im.Minuit(chi2, 
            #F1=0, F2=0, F3=0, F4=0, F5=0,
            F1=Fini[0], F2=Fini[1], F3=Fini[2], F4=Fini[3], F5=Fini[4],            
            #F1=1, F2=2, F3=3, F4=4, F5=5,
            #F1=5, F2=4, F3=3, F4=2, F5=1,
            error_F1=0.001, 
            error_F2=0.001, 
            error_F3=0.001, 
            error_F4=0.001, 
            error_F5=0.001, 
            print_level=1) 

        m.migrad()
        vals = m.values

        F1, F2, F3, F4, F5 = vals['F1'], vals['F2'], vals['F3'], vals['F4'], vals['F5']
        chi2_result.append(chi2(F1, F2, F3, F4, F5)) 
        __coeff = np.array(    [F1, F2, F3, F4, F5])

    ######################################
    if run_mode == 'CCsame':

        def chi2(F1, F2, F3, F4, F5, F6):
            __coeff = np.array([F1, F2, F3, F4, F5, F6])
            MR = (np.dot( vecs, __coeff ) - xsec) / xsec
            ML = np.transpose(MR)
            return np.dot(ML, MR)

        m = im.Minuit(chi2, 
            #F1=0, F2=0, F3=0, F4=0, F5=0, F6=0,
            F1=Fini[0], F2=Fini[1], F3=Fini[2], F4=Fini[3], F5=Fini[4], F6=Fini[5],                           
            #F1=fi1, F2=fi2, F3=fi3, F4=fi4, F5=fi5, F6=fi6, 
            #F1=1, F2=2, F3=3, F4=4, F5=5, F6=6,                
            #F1=6, F2=5, F3=4, F4=3, F5=2, F6=1,                                
            error_F1=0.001, 
            error_F2=0.001, 
            error_F3=0.001, 
            error_F4=0.001, 
            error_F5=0.001, 
            error_F6=0.001,                 
            print_level=1) 

        m.migrad()
        vals = m.values

        F1, F2, F3, F4, F5, F6 = vals['F1'], vals['F2'], vals['F3'], vals['F4'], vals['F5'], vals['F6']
        chi2_result.append(chi2(F1, F2, F3, F4, F5, F6)) 
        __coeff = np.array(    [F1, F2, F3, F4, F5, F6])

    ######################################
    if run_mode in ['C2+C1-', 'C2-C1+']:

        def chi2(F1, F2, F3):
            __coeff = np.array([F1, F2, F3])
            MR = (np.dot( vecs, __coeff ) - xsec) / xsec
            ML = np.transpose(MR)
            return np.dot(ML, MR)

        m = im.Minuit(chi2, 
            #F1=0, F2=0, F3=0,
            F1=Fini[0], F2=Fini[1], F3=Fini[2],             
            error_F1=0.001, 
            error_F2=0.001, 
            error_F3=0.001, 
            print_level=1) 

        m.migrad()
        vals = m.values

        F1, F2, F3 = vals['F1'], vals['F2'], vals['F3']
        chi2_result.append(chi2(F1, F2, F3)) 
        __coeff = np.array(    [F1, F2, F3])

    ######################################
    ######################################

    return chi2_result, __coeff



def get_mixings_2(blocks):
    N11 = blocks['NMIX'].entries[1][1]
    N12 = blocks['NMIX'].entries[1][2]
    N13 = blocks['NMIX'].entries[1][3]
    N14 = blocks['NMIX'].entries[1][4]
    N21 = blocks['NMIX'].entries[2][1]
    N22 = blocks['NMIX'].entries[2][2]
    N23 = blocks['NMIX'].entries[2][3]
    N24 = blocks['NMIX'].entries[2][4]
    N31 = blocks['NMIX'].entries[3][1]
    N32 = blocks['NMIX'].entries[3][2]
    N33 = blocks['NMIX'].entries[3][3]
    N34 = blocks['NMIX'].entries[3][4]
    N41 = blocks['NMIX'].entries[4][1]
    N42 = blocks['NMIX'].entries[4][2]
    N43 = blocks['NMIX'].entries[4][3]
    N44 = blocks['NMIX'].entries[4][4]
    N = [[N11, N12, N13, N14],
        [N21, N22, N23, N24],
        [N31, N32, N33, N34],
        [N41, N42, N43, N44]]
    N = np.array(N)

    V11 = blocks['VMIX'].entries[1][1];
    V12 = blocks['VMIX'].entries[1][2];
    V21 = blocks['VMIX'].entries[2][1];
    V22 = blocks['VMIX'].entries[2][2];
    V = [[V11, V12],[V21, V22]]
    V = np.array(V)

    U11 = blocks['UMIX'].entries[1][1];
    U12 = blocks['UMIX'].entries[1][2];
    U21 = blocks['UMIX'].entries[2][1];
    U22 = blocks['UMIX'].entries[2][2];
    U = [[U11, U12],[U21, U22]]
    U = np.array(U)
    return N, V, U


def get_mixings(slhafile):

    for line in open(slhafile):
        elems = line.split()
        if len(elems) < 3: continue
        if elems[-1] == 'N_11': N11 = float(elems[2])
        if elems[-1] == 'N_12': N12 = float(elems[2])
        if elems[-1] == 'N_13': N13 = float(elems[2])
        if elems[-1] == 'N_14': N14 = float(elems[2])
        if elems[-1] == 'N_21': N21 = float(elems[2])
        if elems[-1] == 'N_22': N22 = float(elems[2])
        if elems[-1] == 'N_23': N23 = float(elems[2])
        if elems[-1] == 'N_24': N24 = float(elems[2])
        if elems[-1] == 'N_31': N31 = float(elems[2])
        if elems[-1] == 'N_32': N32 = float(elems[2])
        if elems[-1] == 'N_33': N33 = float(elems[2])
        if elems[-1] == 'N_34': N34 = float(elems[2])
        if elems[-1] == 'N_41': N41 = float(elems[2])
        if elems[-1] == 'N_42': N42 = float(elems[2])
        if elems[-1] == 'N_43': N43 = float(elems[2])
        if elems[-1] == 'N_44': N44 = float(elems[2])

        if elems[-1] == 'V_11': V11 = float(elems[2])
        if elems[-1] == 'V_12': V12 = float(elems[2])
        if elems[-1] == 'V_21': V21 = float(elems[2])
        if elems[-1] == 'V_22': V22 = float(elems[2])

        if elems[-1] == 'U_11': U11 = float(elems[2])
        if elems[-1] == 'U_12': U12 = float(elems[2])
        if elems[-1] == 'U_21': U21 = float(elems[2])
        if elems[-1] == 'U_22': U22 = float(elems[2])

    N = [[N11, N12, N13, N14],
        [N21, N22, N23, N24],
        [N31, N32, N33, N34],
        [N41, N42, N43, N44]]
    N = np.array(N)

    V = [[V11, V12],[V21, V22]]
    V = np.array(V)

    U = [[U11, U12],[U21, U22]]
    U = np.array(U)
    return N, V, U











