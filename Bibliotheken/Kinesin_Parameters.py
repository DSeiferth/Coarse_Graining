import numpy as np

from matplotlib import pyplot as plt
import math

import time
import random
import networkx as nx
import os

#import sys
#sys.path.append("/home/david.seiferth/Desktop/MA_Arbeit/JupNotebook/Bibliotheken/")

import Steady_State_Calculation_Spanning_Trees as auto



def Kinesin(force):
    #Parametrisierung von Kinesin vgl. Liepelt
    #force 
    F=force#1.e-12 #pico Newtons Groessenordnung
    #Konzentrationen P, ADP, ATP
    cP=1.e-6 #M=Mol?
    cADP=cP
    cATP=cP
    #stepping size l=8nm
    l=8.e-9
    #BoltzmannFaktor Ruamtemperatur [Joule]
    kT=4.e-21

    #rate constants (see Table1 Liepelt) Expt. [9]
    k25=3.0e5 # dimension 1/s
    k52=0.24
    k21=100.
    k56=200.
    k61=k56
    k65=0.09e6 # dimension 1/(micro M s) = 1e6 1/Mol*s
    k16=0.02e6
    k12=1.8e6


    #equlibirum konstancts
    Keq=4.9e11*1e-6 #M=mol
    Keq1=k25*k56*k61*k12/(k52*k65*k16*k21) # corresponds to forward cycle

    k23=k56 # from ADP+ATP to empty + ATP
    k32=k65 # from empty+ATP to ADP +ATP
    k34=k61 #from ATP to ADP (rest empty)
    k43=k16 #from ADP to ATP (rest empyt)
    k45=k12 # from empty +ADP +to ATP +ADP
    #k54=   k52*k23*k34*k45/(k25*k32*k43*Keq1)
    k54=(k52/k25)**2*k21
    #print([k54, k21])

    #equlibirum konstancts
    Keq2=k52*k23*k34*k45/(k25*k32*k43*k54) # corresponds to backward cycle
    #print('equilibrium constant')
    #print([Keq, Keq1, Keq2])
    #chemical transitions

    #mechanical tranistions
    theta=0.3 #load distribution factor aacording to Ref[9]
    Fdim=l*F/kT
    #print('Fdim '+str(Fdim))
    #transition matrix
    w=np.zeros((7,7))
    w[2][5]=k25*math.exp(-theta*Fdim)
    w[5][2]=k52*math.exp((1-theta)*Fdim)
    #chemical transitions
    XI61=0.05 #they are symmetric
    XI56=XI61
    XI23=XI61
    XI34=XI61
    XI12=0.25 #according to Ref[9]
    XI45=XI12

    w[1][2]=k12*cATP*2./(1+math.exp(XI12*Fdim))
    w[2][1]=k21*     2./(1+math.exp(XI12*Fdim))

    w[4][5]=k45*cATP*2./(1+math.exp(XI45*Fdim))
    w[5][4]=k54*     2./(1+math.exp(XI45*Fdim))

    w[6][5]=k65*cADP*2./(1+math.exp(XI56*Fdim))
    w[5][6]=k56*     2./(1+math.exp(XI56*Fdim))

    w[3][2]=k32*cADP*2./(1+math.exp(XI23*Fdim))
    w[2][3]=k23*     2./(1+math.exp(XI23*Fdim))

    w[1][6]=k16*cP  *2./(1+math.exp(XI61*Fdim))
    w[6][1]=k61*     2./(1+math.exp(XI61*Fdim))

    w[4][3]=k43*cP  *2./(1+math.exp(XI34*Fdim))
    w[3][4]=k34*     2./(1+math.exp(XI34*Fdim))

    #print [w[1][2] , w[2][1]]
    #print [w[4][5] , w[5][4]]
    #print [w[6][5] , w[5][6]]
    #print [w[3][2] , w[2][3]]
    #print [w[1][6] , w[6][1]]
    #print [w[4][3] , w[3][4]]
    #print [w[5][2] , w[2][5]]
    
    wsliced= w[:,1:]
    wsliced=np.delete(wsliced, 0, 0)
    #print wsliced
    
    return [wsliced , w]


def Coarse_Grain_Kinesin(w):
    #calculate steady-state
    G=auto.Matrix2Graph(w)
    p=auto.steady_state(G)
    # coarse graining 
    # I start counting from 0
    # first merge 0 and 5

    # wsliced is the original matrix with steady state p
    wcg1=np.zeros((5,5))
    #copying
    for i in range(5):
        for j in range(5):
            wcg1[i][j]=w[i][j]
    # new rates
    wcg1[0][1]=w[0][1]*p[0]/(p[0]+p[5])
    wcg1[0][4]=w[5][4]*p[5]/(p[0]+p[5])
    wcg1[4][0]=w[4][5]
    #print wcg1
    Gcg1=auto.Matrix2Graph(wcg1)
    pcg1=auto.steady_state(Gcg1)
    #print pcg1
    #print p

    # merge states 3 and 4 in the paper notation
    # merge 2 and 3 in notation starting with 0

    #copying
    wcg2=wcg1
    #new rates
    wcg2[4][2]=wcg1[4][3]
    wcg2[2][1]=wcg1[2][1]*pcg1[2]/(pcg1[2]+pcg1[3])
    wcg2[2][4]=wcg1[3][4]*pcg1[3]/(pcg1[2]+pcg1[3])

    #print wcg2
    #delete row and column 3
    #print wcg2[:,3] #column
    #print wcg2[3] #row
    wcg2=np.delete(wcg2, 3, 0)
    #print wcg2

    wcg2= np.delete(wcg2, 3, 1)
    #print wcg2
    Gcg2=auto.Matrix2Graph(wcg2)
    pcg2=auto.steady_state(Gcg2)
    #print 
    #print 'probabilites'
    #print pcg2
    #print pcg1
    #print p
    return wcg2

def Kinesin_CycleFluxes(w, w1):
    G=auto.Matrix2Graph(w)
    normierung=auto.Normfactor(G)
    #print normierung

    #print 'forward cycle F'
    zuflussF=(w1[4][5]*w1[3][2]+w1[3][4]*w1[4][5]+w1[4][3]*w1[3][2])
    JF_plus=w1[1][2]*w1[2][5]*w1[5][6]*w1[6][1]*zuflussF/normierung
    JF_minus=w1[1][6]*w1[6][5]*w1[5][2]*w1[2][1]*zuflussF/normierung


    #print 'backward cycle B'
    zuflussB=(w1[1][2]*w1[6][5]+w1[6][1]*w1[1][2]+w1[1][6]*w1[6][5])
    JB_plus=w1[2][3]*w1[3][4]*w1[4][5]*w1[5][2]*zuflussB/normierung
    JB_minus=w1[2][5]*w1[5][4]*w1[4][3]*w1[3][2]*zuflussB/normierung

    #print 'dissipative cycle D'
    zuflussD=1.0
    JD_plus=w1[1][2]*w1[2][3]*w1[3][4]*w1[4][5]*w1[5][6]*w1[6][1]*zuflussD/normierung
    JD_minus=w1[1][6]*w1[6][5]*w1[5][4]*w1[4][3]*w1[3][2]*w1[2][1]*zuflussD/normierung


    JF=JF_plus-JF_minus
    JB=JB_plus-JB_minus
    JD=JD_plus-JD_minus
    
    mean=[JF, JB, JD]
    var=[JF_plus+JF_minus , JB_plus+JB_minus, JD_plus+JD_minus]
    print('mean and variance of Kinesin cycles fluxes')
    print(mean)
    print(var)
    return [mean, var]


def Kinsein_CG_CycleFluxes(wcg2):
    Gcg2=auto.Matrix2Graph(wcg2)
    normierungCG2=auto.Normfactor(Gcg2)
    #print normierungCG2
    #steady state cycle fluxes
    #print 'CG2 forward cycle F'
    zuflussFCG2=(wcg2[2][1]+wcg2[2][3])
    JF_plusCG2=wcg2[0][1]*wcg2[1][3]*wcg2[3][0]*zuflussFCG2/normierungCG2
    JF_minusCG2=wcg2[0][3]*wcg2[3][1]*wcg2[1][0]*zuflussFCG2/normierungCG2


    #print 'CG2 backward cycle B'
    zuflussBCG2=(wcg2[0][3]+wcg2[0][1])
    JB_plusCG2=wcg2[1][2]*wcg2[2][3]*wcg2[3][1]*zuflussBCG2/normierungCG2
    JB_minusCG2=wcg2[1][3]*wcg2[3][2]*wcg2[2][1]*zuflussBCG2/normierungCG2

    #print 'dissipative cycle D'
    zuflussDCG2=1.0
    JD_plusCG2=wcg2[0][1]*wcg2[1][2]*wcg2[2][3]*wcg2[3][0]*zuflussDCG2/normierungCG2
    JD_minusCG2=wcg2[0][3]*wcg2[3][2]*wcg2[2][1]*wcg2[1][0]*zuflussDCG2/normierungCG2
    
    mean=[JF_plusCG2-JF_minusCG2 , JB_plusCG2-JB_minusCG2, JD_plusCG2-JD_minusCG2]
    var=[JF_plusCG2+JF_minusCG2 , JB_plusCG2+JB_minusCG2, JD_plusCG2+JD_minusCG2]
    print('mean and variance of coarse-grained cycles fluxes')
    print(mean)
    print(var)
    return [mean, var]

def Kinesin_OneWayCycleFluxes(w, w1):
    G=auto.Matrix2Graph(w)
    normierung=auto.Normfactor(G)
    #print normierung

    #print 'forward cycle F'
    zuflussF=(w1[4][5]*w1[3][2]+w1[3][4]*w1[4][5]+w1[4][3]*w1[3][2])
    JF_plus=w1[1][2]*w1[2][5]*w1[5][6]*w1[6][1]*zuflussF/normierung
    JF_minus=w1[1][6]*w1[6][5]*w1[5][2]*w1[2][1]*zuflussF/normierung


    #print 'backward cycle B'
    zuflussB=(w1[1][2]*w1[6][5]+w1[6][1]*w1[1][2]+w1[1][6]*w1[6][5])
    JB_plus=w1[2][3]*w1[3][4]*w1[4][5]*w1[5][2]*zuflussB/normierung
    JB_minus=w1[2][5]*w1[5][4]*w1[4][3]*w1[3][2]*zuflussB/normierung

    #print 'dissipative cycle D'
    zuflussD=1.0
    JD_plus=w1[1][2]*w1[2][3]*w1[3][4]*w1[4][5]*w1[5][6]*w1[6][1]*zuflussD/normierung
    JD_minus=w1[1][6]*w1[6][5]*w1[5][4]*w1[4][3]*w1[3][2]*w1[2][1]*zuflussD/normierung


    JF=JF_plus-JF_minus
    JB=JB_plus-JB_minus
    JD=JD_plus-JD_minus
    
    return [[JF_plus,JF_minus],[JB_plus,JB_minus],[JD_plus,JD_minus]]


def Kinsein_CG_OneWayCycleFluxes(wcg2):
    Gcg2=auto.Matrix2Graph(wcg2)
    normierungCG2=auto.Normfactor(Gcg2)
    #print normierungCG2
    #steady state cycle fluxes
    #print 'CG2 forward cycle F'
    zuflussFCG2=(wcg2[2][1]+wcg2[2][3])
    JF_plusCG2=wcg2[0][1]*wcg2[1][3]*wcg2[3][0]*zuflussFCG2/normierungCG2
    JF_minusCG2=wcg2[0][3]*wcg2[3][1]*wcg2[1][0]*zuflussFCG2/normierungCG2


    #print 'CG2 backward cycle B'
    zuflussBCG2=(wcg2[0][3]+wcg2[0][1])
    JB_plusCG2=wcg2[1][2]*wcg2[2][3]*wcg2[3][1]*zuflussBCG2/normierungCG2
    JB_minusCG2=wcg2[1][3]*wcg2[3][2]*wcg2[2][1]*zuflussBCG2/normierungCG2

    #print 'dissipative cycle D'
    zuflussDCG2=1.0
    JD_plusCG2=wcg2[0][1]*wcg2[1][2]*wcg2[2][3]*wcg2[3][0]*zuflussDCG2/normierungCG2
    JD_minusCG2=wcg2[0][3]*wcg2[3][2]*wcg2[2][1]*wcg2[1][0]*zuflussDCG2/normierungCG2
    
    return [[JF_plusCG2,JF_minusCG2] , [JB_plusCG2,JB_minusCG2], [JD_plusCG2,JD_minusCG2]]

