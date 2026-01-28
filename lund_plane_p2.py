import numpy as np
import matplotlib.pyplot as plt
import random as random
import sys

from shower import *
from particle import CheckEvent
from qcd import AlphaS
from matrix import eetojj

alphas = AlphaS(91.1876,0.118)
hardxs = eetojj(alphas)
shower = Shower(alphas,t0=1.)

ecom = 91.1876

def Magnitude3(mom):
    return m.sqrt(mom[1]**2 + mom[2]**2 + mom[3]**2)

def EmissionAngle(mom1, mom2):
    dot = mom1[1]*mom2[1] + mom1[2]*mom2[2] + mom1[3]*mom2[3]
    return m.acos(dot / (Magnitude3(mom1 * Magnitude3(mom2))))

def emissionsN(evno, emno):
    
    for ev in evno:
    
        moms_list = []
        em_angle = []
        em_energy = []
        
        event, weight = hardxs.GenerateLOPoint()
        t = (event[0].mom+event[1].mom).M2()
        shower.c = 1
        shower.t = t
        
        
        while shower.t > shower.t0:
            moms = shower.GeneratePoint(event)
            if moms == None: continue
            moms_list.append(moms)
            
        if not CheckEvent(event):
            print("Something went wrong:")
            for p in event:
                print(p)
        
        for em in range(0, len(moms_list)):
            if em >= emno: break
            pi = moms_list[em][0]
            pj = moms_list[em][1]
            pk = moms_list[em][2]
            
            angle = EmissionAngle(pi,pj) / np.pi
            energy = 2*pj[0]/ecom
            
            if energy > 1- angle or energy <= 0 or angle <= 0: break
            
            em_angle.append(angle)
            em_energy.append(energy)

        sys.stdout.write('\rEvent {0}'.format(ev))
        sys.stdout.flush()
            
            
            
    
            
        
        
        
    
    