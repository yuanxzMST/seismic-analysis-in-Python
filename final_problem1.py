# -*- coding: utf-8 -*-
"""
Created on Fri May 11 12:54:09 2018

@author: xyvm4

6205 structural dynamics & earthquake engineering

Final Project Problem-1

Xinzhe Yuan
"""

import numpy as np
import math
import matplotlib.pyplot as plt
################################################################################

def Newmark_coefficients(dt,garma,belta):

    a0=1/(belta*(dt**2))
    a1=garma/(belta*dt)
    a2=1/(belta*dt)
    a3=(1/(2*belta))-1
    a4=(garma/belta)-1
    a5=(dt/2)*((garma/belta)-2)
    a6=dt*(1-garma)
    a7=garma*dt

    return a0,a1,a2,a3,a4,a5,a6,a7

################################################################################
# input values

garma = 0.5
belta = 0.25

OmegaP = 0.5
TP = 2*math.pi/OmegaP
TN = TP*np.arange(10,0.1,-0.05)
TP_over_TN = TP/TN

stiffness = 1
mass = TN**2*stiffness/(2*math.pi)**2
massT = zip(mass,TN)
vel_s = 30 #in/sec
damping_ratio = 0.05

loaddata = np.loadtxt('excitation.txt')
inD = 0
dT = 0.02

ratio_force = np.array([0.1,0.2,0.3,0.5,0.7,0.9,1.0,1.2,1.5,2.0])
epslon = 1e-3
 
elasticS = stiffness
plasticS = 0
a0,a1,a2,a3,a4,a5,a6,a7 = Newmark_coefficients(dT,garma,belta)
umax = np.zeros((len(massT),len(ratio_force)))
DuctileFactor = np.zeros((len(massT),len(ratio_force)))

################################################################################
# computation

jj = 0
for ratio in ratio_force:
    ii = 0
    for item in massT:
        yieldforce = ratio*item[0]*vel_s*OmegaP
        yielddisp = yieldforce/stiffness
        eq_force = -item[0]*loaddata
        damping = 2*damping_ratio*item[0]*2*math.pi/item[1]
        steps = len(eq_force)-1
        R = np.zeros(steps)
        kT = np.zeros(steps+1); kT[0] = stiffness
        kTilter = np.zeros(steps)
        deltU = np.zeros(steps)
        u = np.zeros(steps+1); u[0] = inD
        fs = np.zeros(steps+1); fs[0] = u[0]*kT[0]
        ud = np.zeros(steps+1); ud[0] = 0
        udd = np.zeros(steps+1); udd[0] = loaddata[0]
        pTilter = np.zeros(steps)
        pStep = eq_force
        cA = a0*item[0] + a1*damping
        cB = a2*item[0] + a4*damping
        cC = a3*item[0] + a5*damping
        
        for i in range(steps):
            u[i+1] = u[i]
            kT[i+1] = kT[i]
            fs[i+1] = fs[i]
            pTilter[i] = pStep[i+1]+cA*u[i]+cB*ud[i]+cC*udd[i]
            R[i] = pTilter[i]-fs[i+1]-cA*u[i+1]
            
            while abs(R[i]) > epslon:
                 kTilter[i] = kT[i+1]+cA
                 deltU[i] = R[i]/kTilter[i]
                 u[i+1] = u[i+1] + deltU[i]
                 fs[i+1] = fs[i+1]+kT[i+1]*deltU[i]
                 
                 if fs[i+1] >= yieldforce:
                     kT[i+1] = plasticS
                     fs[i+1] = yieldforce
                 elif fs[i+1] <= -yieldforce:
                     kT[i+1] = plasticS
                     fs[i+1] = -yieldforce
                 else:
                     kT[i+1] = elasticS
                 R[i] = pTilter[i]-fs[i+1]-cA*u[i+1]    
                
        
            ud[i+1] = a1*(u[i+1]-u[i])-a4*ud[i]-a5*udd[i]
            udd[i+1] = a0*(u[i+1]-u[i])-a2*ud[i]-a3*udd[i]
        
            if ud[i+1]*ud[i] <= 0:
                kT[i+1] = elasticS
         
        umax[ii][jj] = max(abs(u))
        DuctileFactor[ii][jj] = umax[ii][jj]/yielddisp
        ii += 1
    jj += 1

################################################################################
# Chart plot

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)
for i in range(len(ratio_force)):
    text = r'R_u/ms$\omega_p$={}'.format(ratio_force[i])
    line, = ax.plot(TP_over_TN,DuctileFactor[:,i])
    ax.annotate(text, xy=(TP_over_TN[-1],DuctileFactor[:,i][-1]))
    ax.set_yscale('log')
    ax.set_xscale('log')

plt.xlabel('$T_p/T_n$')
plt.ylabel(r'$\mu$=y_m/y_e')
plt.title(r'Graphical Solution Chart for $\omega_p$=0.4, s=50')
fig.savefig('finalpl1.png')

        
        
        
                     
            
        
        
        
