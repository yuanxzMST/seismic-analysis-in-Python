# -*- coding: utf-8 -*-
"""
Created on Fri May 11 12:54:09 2018

@author: xyvm4

6205 structural dynamics & earthquake engineering

Final Project

Xinzhe Yuan
"""

import numpy as np
import math

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

#   input valuse
garma = 0.5
belta = 0.25

OmegaP = 0.4
TP = 2*math.pi/OmegaP
TN = TP*np.arange(10,0.0,-0.03)
TP_over_TN = TP/TN

stiffness = 1
mass = TN**2*stiffness/(2*math.pi)**2

damping = 0.1592
inD = 0
inV = 0
inF = 0
dT = 0.1
yieldforce = 7.5
elasticS = 10
plasticS = 0