import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scipy.optimize import curve_fit
import pandas as pd

from shower import *
from particle import CheckEvent
from qcd import AlphaS
from matrix import eetojj

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

ecom = 91.1876

alphas = AlphaS(ecom,0.118)
hardxs = eetojj(alphas)
shower = Shower(alphas,t0=1.)

def Magnitude3(mom):
    return m.sqrt(mom[1]**2 + mom[2]**2 + mom[3]**2)
    
def EmissionAngle(mom1, mom2):
    dotprod = mom1[1]*mom2[1]+mom1[2]*mom2[2]+mom1[3]*mom2[3]
    return m.acos(dotprod / (Magnitude3(mom1) * Magnitude3(mom2)))

def emissionsy(evno, log = False):

    em_angle = np.zeros(evno)
    em_energy = np.zeros(evno)

    for ev in range(evno):
        event, weight = hardxs.GenerateLOPoint()
        t = (event[0].mom+event[1].mom).M2()
        while len(event) < 5:
            shower.c = 1
            shower.t = t
            moms = shower.GeneratePoint(event)
        if not CheckEvent(event):
            print("Something went wrong:")
            for p in event:
                print(p)
        pi = moms[0]
        pj = moms[1]
        em_angle[ev] = EmissionAngle(pi, pj) / np.pi
        em_energy[ev] = 2*pj[0] / ecom
        if em_energy[ev]*(em_angle[ev])**2>0.0004:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan
        elif em_energy[ev] < 0 or em_angle[ev] < 0:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan

        progress_bar(ev + 1, evno)
    
    return em_angle, em_energy


def emissionsx(evno, log = False):

    em_angle = np.zeros(evno)
    em_energy = np.zeros(evno)

    for ev in range(evno):
        event, weight = hardxs.GenerateLOPoint()
        t = (event[0].mom+event[1].mom).M2()
        while len(event) < 5:
            shower.c = 1
            shower.t = t
            moms = shower.GeneratePoint(event)
        if not CheckEvent(event):
            print("Something went wrong:")
            for p in event:
                print(p)
        pi = moms[0]
        pj = moms[1]
        em_angle[ev] = EmissionAngle(pi, pj) / np.pi
        em_energy[ev] = 2*pj[0] / ecom
        if em_energy[ev]*em_angle[ev]>0.009:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan
        elif em_energy[ev] < 0 or em_angle[ev] < 0:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan

        progress_bar(ev + 1, evno)
    
    return em_angle, em_energy

def psboundary(angle, scale_factor):
    curve = scale_factor*pow(angle, -2)
    return curve

em_angle1, em_energy1 = emissionsy(10000, False)
em_angle2, em_energy2 = emissionsx(10000, True)

em_angle = np.append(em_angle1, em_angle2)
em_energy = np.append(em_energy1, em_energy2)

df = pd.DataFrame({"Angle" : em_angle, "Energy" : em_energy})
cleaned_df = df.dropna()
cleaned_angle = cleaned_df["Angle"].to_numpy()
cleaned_energy = cleaned_df["Energy"].to_numpy()

popt, pcov = curve_fit(psboundary, cleaned_angle, cleaned_energy)
print(f'{popt[0]:2f}')
#emissionsN(10, 5)
offset = 0.00000001

x1, y1 = np.linspace(offset, 1, 100), np.full(100, offset)
x2, y2 = np.linspace(offset, 1, 100), 1 - np.linspace(offset, 1, 100)
x3, y3 = np.full(100, offset), np.linspace(offset, 1, 100)
angle = np.linspace(offset, 1, 1000)
plt.ylim(-0.01, 1)
plt.plot(x1, y1, color="k", linestyle="--", label="Phase space boundary")
plt.plot(x2, y2, color="k", linestyle="--")
plt.plot(x3, y3, color="k", linestyle="--")
plt.plot(angle, psboundary(angle, *popt))

plt.scatter(em_angle, em_energy, s=1, label = "Emissions")

plt.xlabel(r"$\theta/\pi$")
plt.ylabel(r"$2E/Q$")
plt.title("Lund jet plane for first emission")
plt.legend()
plt.savefig("lund_scipy_fit.pdf")
plt.show()