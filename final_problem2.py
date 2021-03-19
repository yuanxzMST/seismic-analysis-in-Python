# -*- coding: utf-8 -*-
"""
Created on Sat May 12 16:10:49 2018

@author: xyvm4

6205 structural dynamics & earthquake engineering

Final Project Problem-2

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
TN = TP*0.5
stiffness = 1
mass = TN**2*stiffness/(2*math.pi)**2
vel_s = 30 #in/sec
damping_ratio = 0.05
loaddata = np.loadtxt('excitation.txt')
eq_force = -mass*loaddata
pStep = eq_force
inD = 0
dT = 0.02
ratio_force = np.array([0.2,0.5,1.0])
epslon = 1e-3
elasticS = stiffness
a0,a1,a2,a3,a4,a5,a6,a7 = Newmark_coefficients(dT,garma,belta)
damping = 2*damping_ratio*mass*2*math.pi/TN
steps = len(eq_force)-1
cA = a0*mass + a1*damping
cB = a2*mass + a4*damping
cC = a3*mass + a5*damping

gamma = np.array([0.1,0.2,0.3,0.4])
col_num = np.array([2,3,4,5])

umax = np.zeros((len(col_num),len(ratio_force)))
DuctileFactor = np.zeros((len(col_num),len(ratio_force)))
solu_23 = np.zeros((len(gamma),len(col_num),len(ratio_force)))

for ii in range(len(gamma)):
#    kk = 0
    gar = gamma[ii]
    for kk in range(len(ratio_force)):
        ratio = ratio_force[kk]
#        yieldforce = ratio*mass*OmegaP*vel_s
#        yielddisp = yieldforce/stiffness
        
        for jj in range(len(col_num)):
            cnum = col_num[jj]
            R = np.zeros(steps)
            kT = np.zeros(steps+1); kT[0] = stiffness
            kTilter = np.zeros(steps)
            deltU = np.zeros(steps)
            u = np.zeros(steps+1); u[0] = inD
            fs = np.zeros(steps+1); fs[0] = u[0]*kT[0]
            ud = np.zeros(steps+1); ud[0] = 0
            udd = np.zeros(steps+1); udd[0] = loaddata[0]
            pTilter = np.zeros(steps)
            yieldforce = ratio*mass*OmegaP*vel_s
            
            if cnum == 2:
                plasticS = 0
                yielddisp = yieldforce/stiffness
                
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
                
                umax[jj][kk] = max(abs(u))
                DuctileFactor[jj][kk] = umax[jj][kk]/yielddisp
            
            if cnum == 3:
                
                plasticS1 = (2-gamma[ii])*stiffness/3; force1 = yieldforce/(1+gamma[ii])
                plasticS2 = (1-gamma[ii])*stiffness/3; force2 = (3-gamma[ii])*yieldforce/3
                plasticS3 = 0
                yielddisp = force1/stiffness + (force2-force1)/plasticS1 + (yieldforce-force2)/plasticS2
                
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
                        
                        if fs[i+1] >= -force1 and fs[i+1] <= force1:
                            kT[i+1] = elasticS
                        elif fs[i+1] >= force1 and fs[i+1] <= force2:
                            kT[i+1] = plasticS1
                        elif fs[i+1] >= force2 and fs[i+1] <= yieldforce:
                            kT[i+1] = plasticS2
                        elif fs[i+1] >= yieldforce:
                            kT[i+1] = plasticS3
                            fs[i+1] = yieldforce
                        elif fs[i+1] >= -force2 and fs[i+1] <= -force1:
                            kT[i+1] = plasticS1
                        elif fs[i+1] >= -yieldforce and fs[i+1] <= -force2:
                            kT[i+1] = plasticS2
                        elif fs[i+1] <= -yieldforce:
                            kT[i+1] = plasticS3
                            fs[i+1] = -yieldforce
                        R[i] = pTilter[i]-fs[i+1]-cA*u[i+1]
                    
                    ud[i+1] = a1*(u[i+1]-u[i])-a4*ud[i]-a5*udd[i]
                    udd[i+1] = a0*(u[i+1]-u[i])-a2*ud[i]-a3*udd[i]
                    
                    if ud[i+1]*ud[i] <= 0:
                        kT[i+1] = elasticS
                        
                umax[jj][kk] = max(abs(u))
                DuctileFactor[jj][kk] = umax[jj][kk]/yielddisp
          
            if cnum == 4:
                
                plasticS1 = (3-gamma[ii])/4; force1 = yieldforce/(1+gamma[ii])
                plasticS2 = (1-gamma[ii])/4; force2 = (4-gamma[ii])*yieldforce/4
                plasticS3 = 0;
                yielddisp = force1/stiffness + (force2-force1)/plasticS1 + (yieldforce-force2)/plasticS2
                
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
                        
                        if fs[i+1] >= -force1 and fs[i+1] <= force1:
                            kT[i+1] = elasticS
                        elif fs[i+1] >= force1 and fs[i+1] <= force2:
                            kT[i+1] = plasticS1
                        elif fs[i+1] >= force2 and fs[i+1] <= yieldforce:
                            kT[i+1] = plasticS2
                        elif fs[i+1] >= yieldforce:
                            kT[i+1] = plasticS3
                            fs[i+1] = yieldforce
                        elif fs[i+1] >= -force2 and fs[i+1] <= -force1:
                            kT[i+1] = plasticS1
                        elif fs[i+1] >= -yieldforce and fs[i+1] <= -force2:
                            kT[i+1] = plasticS2
                        elif fs[i+1] <= -yieldforce:
                            kT[i+1] = plasticS3
                            fs[i+1] = -yieldforce
                        R[i] = pTilter[i]-fs[i+1]-cA*u[i+1]
                    
                    ud[i+1] = a1*(u[i+1]-u[i])-a4*ud[i]-a5*udd[i]
                    udd[i+1] = a0*(u[i+1]-u[i])-a2*ud[i]-a3*udd[i]
                    
                    if ud[i+1]*ud[i] <= 0:
                        kT[i+1] = elasticS
                        
                umax[jj][kk] = max(abs(u))
                DuctileFactor[jj][kk] = umax[jj][kk]/yielddisp
                
            if cnum == 5:
                
                plasticS1 = (4-2*gamma[ii])*stiffness/5; force1 = yieldforce/(1+2*gamma[ii])
                plasticS2 = (3-3*gamma[ii])*stiffness/5; force2 = (5-gamma[ii])*yieldforce/(1+gamma[ii])/5
                plasticS3 = (2-3*gamma[ii])*stiffness/5; force3 = (5-3*gamma[ii])*yieldforce/5
                plasticS4 = (1-2*gamma[ii])*stiffness/5; force4 = (5-6*gamma[ii])/5/(1-gamma[ii])*yieldforce
                plasticS5 = 0
                yielddisp = force1/stiffness + (force2-force1)/plasticS1 + (force3-force2)/plasticS2 + (force4-force3)/plasticS3 + (yieldforce-force4)/plasticS4
                
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
                        
                        if fs[i+1] >= -force1 and fs[i+1] <= force1:
                            kT[i+1] = elasticS
                        elif fs[i+1] >= force1 and fs[i+1] <= force2:
                            kT[i+1] = plasticS1
                        elif fs[i+1] >= force2 and fs[i+1] <= force3:
                            kT[i+1] = plasticS2
                        elif fs[i+1] >= force3 and fs[i+1] <= force4:
                            kT[i+1] = plasticS3
                        elif fs[i+1] >= force4 and fs[i+1] <= yieldforce:
                            kT[i+1] = plasticS4
                        elif fs[i+1] >= yieldforce:
                            kT[i+1] = plasticS5
                            fs[i+1] = yieldforce
                        elif fs[i+1] >= -force2 and fs[i+1] <= -force1:
                            kT[i+1] = plasticS1
                        elif fs[i+1] >= -force3 and fs[i+1] <= -force2:
                            kT[i+1] = plasticS2
                        elif fs[i+1] >= -force4 and fs[i+1] <= -force3:
                            kT[i+1] = plasticS3
                        elif fs[i+1] >= -yieldforce and fs[i+1] <= -force4:
                            kT[i+1] = plasticS4
                        elif fs[i+1] <= -yieldforce:
                            kT[i+1] = plasticS5
                            fs[i+1] = -yieldforce
                        R[i] = pTilter[i]-fs[i+1]-cA*u[i+1]
                    
                    ud[i+1] = a1*(u[i+1]-u[i])-a4*ud[i]-a5*udd[i]
                    udd[i+1] = a0*(u[i+1]-u[i])-a2*ud[i]-a3*udd[i]
                    
                    if ud[i+1]*ud[i] <= 0:
                        kT[i+1] = elasticS
                        
                umax[jj][kk] = max(abs(u))
                DuctileFactor[jj][kk] = umax[jj][kk]/yielddisp
                    
            jj += 1
        kk += 1
        
    solu_23[ii] = DuctileFactor
    
    ii += 1

################################################################################  
# plot the results
# problem 2 for each gamma value
# key = 0
for key in range(len(solu_23)):
    item = solu_23[key]
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)
    for pp in range(len(ratio_force)):
        text = r'$R_u/ms\omega_p={}$'.format(ratio_force[pp])
        line, = ax.plot(col_num,item[:,pp])
        ax.annotate(text, xy=(col_num[-1],item[:,pp][-1]))
        ax.set_yscale('log')
    plt.xlabel('column number')
    plt.ylabel(r'$\mu=u_m/u_y$')
    plt.title(r'$\gamma={}$'.format(gamma[key]))
    fig.savefig('fig{}.png'.format(key+1))
    
# problem 3 for each column number

for kkey in range(len(col_num)):
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)
    itemm = solu_23[:,kkey,:]
    for ppp in range(len(ratio_force)):
        text = r'$R_u/ms\omega_p={}$'.format(ratio_force[ppp])
        line, = ax.plot(gamma,itemm[:,ppp])
        ax.annotate(text, xy = (gamma[-1],itemm[:,ppp][-1]))
        ax.set_yscale('log')
    plt.xlabel(r'$\gamma$')
    plt.ylabel(r'$\mu=u_m/u_y$')
    plt.title('column number={}'.format(col_num[kkey]))
    fig.savefig('figg{}.png'.format(kkey+1))