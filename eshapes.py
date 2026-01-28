"""
Making Your Own Parton Shower in Python - Thrust
------------------------------------------------
This programme generates a Thrust Distribution. This programme is based on
S. Höche's Tutorial on Parton Showers and Matching [1, 2].

[1] - Höche, S. (2014). Introduction to parton-shower event generators.
[2] - The MCnet Collaboration, “MCnet-CTEQ Summer School 2021 VIRTUAL,” Sept. 
2021.
"""

import math as m
from vector import Vec4
from particle import Particle
from histogram import Histo1D


# ------------------------------------------------------------------------------
# Utitlity Functions

def mag(v):
    return m.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])


def mag2(v):
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]


def add(a, b):
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]


def sub(a, b):
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]


def mul(v, s):
    return [v[0]*s, v[1]*s, v[2]*s]


def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def crs(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def determinant(m):
    a, b, c = m[0]
    d, e, f = m[1]
    g, h, i = m[2]
    det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    return det


# ------------------------------------------------------------------------------
class EAlgorithm:

    # --------------------------------------------------------------------------

    def GetMomenta(self, event):

        moms = []
        for p in event[2:]:
            moms.append([p.mom.px, p.mom.py, p.mom.pz])

        moms.sort(key=mag2, reverse=True)
        return moms

    # --------------------------------------------------------------------------

    def CalculateThrust(self, moms):

        # FORTRAN ARRAYS START AT 1
        # PYTHON IGNORES N IN RANGE(N)

        thrust = 0.
        t_axis = [0., 0., 0.]

        for k in range(1, len(moms)):

            for j in range(k):

                tmp_axis = crs(moms[j], moms[k])
                p_thrust = [0., 0., 0.]
                p_combin = []

                for i in range(len(moms)):

                    if (i != j) and (i != k):

                        if (dot(moms[i], tmp_axis) >= 0):

                            p_thrust = add(p_thrust, moms[i])

                        else:

                            p_thrust = sub(p_thrust, moms[i])

                p_combin.append(add(add(p_thrust, moms[j]), moms[k]))
                p_combin.append(sub(add(p_thrust, moms[j]), moms[k]))
                p_combin.append(add(sub(p_thrust, moms[j]), moms[k]))
                p_combin.append(sub(sub(p_thrust, moms[j]), moms[k]))

                for p in p_combin:

                    temp = mag(p)
                    if temp > thrust:
                        thrust = temp
                        t_axis = mul(p, 1/mag(p))

        momsum = sum(mag(mo) for mo in moms)
        thrust /= momsum
        t_axis = mul(t_axis, 1./mag(t_axis))

        if thrust < 0.5:
            print(thrust)

        # Make sure t_axis in +z direction
        t_axis = mul(t_axis, -1.) if t_axis[2] < 0 else t_axis

        return thrust, t_axis

    # --------------------------------------------------------------------------

    def CalculateJetMBr(self, moms, t_axis):

        # Momentum Sum (Should be 91.2 GeV)
        momsum = sum(mag(mo) for mo in moms)

        # Get Thrust Axis
        t_axis = mul(t_axis, 1./mag(t_axis))
        t_axis = mul(t_axis, -1.) if t_axis[2] < 0 else t_axis  # t_axis in +z

        p_with = Vec4()
        p_against = Vec4()

        n_with = 0
        n_against = 0

        e_vis = 0.
        broad_with = 0.
        broad_against = 0.
        broad_denominator = 0.

        for i in range(len(moms)):

            mo_para = dot(moms[i], t_axis)
            mo_perp = mag(sub(moms[i], mul(t_axis, mo_para)))

            e_vis += mag(moms[i])
            broad_denominator += 2. * mag(moms[i])

            if (mo_para > 0.):

                p_with += Vec4(mag(moms[i]), moms[i][0],
                               moms[i][1], moms[i][2])
                broad_with += mo_perp
                n_with += 1

            elif (mo_para < 0.):

                p_against += Vec4(mag(moms[i]), moms[i]
                                  [0], moms[i][1], moms[i][2])
                broad_against += mo_perp
                n_against += 1

            else:

                p_with += 0.5 * \
                    Vec4(mag(moms[i]), moms[i][0], moms[i][1], moms[i][2])
                p_against += 0.5 * \
                    Vec4(mag(moms[i]), moms[i][0], moms[i][1], moms[i][2])
                broad_with += 0.5 * mo_perp
                broad_against += 0.5 * mo_perp
                n_with += 1
                n_against += 1

        e2_vis = m.pow(e_vis, 2.)

        mass2_with = m.fabs(p_with.M2() / e2_vis)
        mass2_against = m.fabs(p_against.M2() / e2_vis)

        mass_with = m.sqrt(mass2_with)
        mass_against = m.sqrt(mass2_against)

        m2H = max(mass_with, mass_against)
        m2L = min(mass_with, mass_against)

        broad_with /= broad_denominator
        broad_against /= broad_denominator

        bMax = max(broad_with, broad_against)
        bMin = min(broad_with, broad_against)

        if (n_with == 1) or (n_against == 1):

            return [m2H, -5., bMax, -5.]

        return [m2H, m2L, bMax, bMin]

    # --------------------------------------------------------------------------

    def CalculateSphero(self, moms):

        # Get Transverse components of the momenta and their sum, and unit vecs
        trmoms = []
        trsum = 0.
        units = []

        for mo in moms:
            tr = [mo[0], mo[1], 0]
            trmoms.append(tr)
            trsum += mag(tr)
            units.append(mul(tr, 1./mag(tr)))

        val = 99999.
        for u in units:

            s = 0.

            for tr in trmoms:

                s += m.fabs(mag(crs(tr, u)))

            if (s < val):
                val = s

        sphero = 0.25 * m.pi * m.pi * val * val / (trsum * trsum)
        return sphero

    # --------------------------------------------------------------------------

    def CalculateEigenvalues(self, t):

        q = (t[0][0] + t[1][1] + t[2][2])/3.
        p1 = t[0][1]*t[0][1] + t[0][2]*t[0][2] + t[1][2]*t[1][2]
        p2 = (t[0][0] - q)**2 + (t[1][1] - q)**2 + (t[2][2] - q)**2 + 2.*p1
        p = m.sqrt(p2/6.)

        t2 = t
        for i in range(3):
            t2[i][i] -= q

        for i in range(3):
            for j in range(3):
                t2[i][j] /= p

        r = determinant(t2) / 2.

        phi = 0
        if (r <= -1):
            phi = m.pi / 3.
        elif (r >= 1):
            phi = 0
        else:
            phi = m.acos(r) / 3.

        l1 = q + 2 * p * m.cos(phi)
        l3 = q + 2 * p * m.cos(phi + (2*m.pi/3.))
        l2 = 3 * q - l1 - l3

        return l1, l2, l3

    def CalculateCDPrms(self, moms):

        momsum = sum(mag(mo) for mo in moms)

        t = [[0., 0., 0.],
             [0., 0., 0.],
             [0., 0., 0.]]

        for a in range(3):
            for b in range(3):

                t[a][b] = sum((mo[a]*mo[b]/mag(mo)) for mo in moms)
                t[a][b] /= momsum

        l1, l2, l3 = self.CalculateEigenvalues(t)

        c = 3.*(l1*l2 + l2*l3 + l3*l1)

        if len(moms) < 4:
            d = -5.
        else:
            d = 27.*l1*l2*l3

        return [c, d]

    # --------------------------------------------------------------------------

    def Calculate(self, event):

        # No 0-emission events allowed
        if (len(event) <= 4):
            return [-5.] * 8

        # Get Momenta, sorted by largest magnitude2 (leading vectors)
        moms = self.GetMomenta(event)

        # Thrust
        thrust, t_axis = self.CalculateThrust(moms)

        # Jet Mass and Broadening
        hjm, ljm, wjb, njb = self.CalculateJetMBr(moms, t_axis)

        # Spherocity
        sph = self.CalculateSphero(moms)

        # C and D Parameters
        cpr, dpr = self.CalculateCDPrms(moms)

        # Same Coplanarity Check for D Param as njb and wjb
        # Instead of redoing check, just remove if njb/ljm = -5.
        if (ljm == -5.) or (njb == -5.):
            dpr = -5.

        return [1. - thrust, hjm, ljm, wjb, njb, sph, cpr, dpr]


# Histogram Generation

class EAnalysis:

    def __init__(self, n=1):

        self.wtot = 0.
        self.thr = Histo1D(50, 0., .5, '/EShapes/tvalue\n')
        self.tzm = Histo1D(50, 0., .5, '/EShapes/tzoomd\n')
        self.hjm = Histo1D(100, 0., 1., '/EShapes/hjm\n')
        self.ljm = Histo1D(100, 0., .5, '/EShapes/ljm\n')
        self.wjb = Histo1D(100, 0., .5, '/EShapes/wjb\n')
        self.njb = Histo1D(100, 0., .2, '/EShapes/njb\n')
        self.sph = Histo1D(100, 0., 1., '/EShapes/sph\n')
        self.cpr = Histo1D(100, 0., 1., '/EShapes/cpr\n')
        self.dpr = Histo1D(100, 0., 1., '/EShapes/dpr\n')
        self.eshape = EAlgorithm()

    def Analyze(self, event, w, weights=[]):

        arr = self.eshape.Calculate(event)

        self.wtot += w
        self.thr.Fill(arr[0], w)
        self.tzm.Fill(arr[0], w)
        self.hjm.Fill(arr[1], w)
        self.ljm.Fill(arr[2], w)
        self.wjb.Fill(arr[3], w)
        self.njb.Fill(arr[4], w)
        self.sph.Fill(arr[5], w)
        self.cpr.Fill(arr[6], w)
        self.dpr.Fill(arr[7], w)

    def Finalize(self, name):

        self.thr.ScaleW(1./self.wtot)
        self.tzm.ScaleW(1./self.wtot)
        self.hjm.ScaleW(1./self.wtot)
        self.ljm.ScaleW(1./self.wtot)
        self.wjb.ScaleW(1./self.wtot)
        self.njb.ScaleW(1./self.wtot)
        self.sph.ScaleW(1./self.wtot)
        self.cpr.ScaleW(1./self.wtot)
        self.dpr.ScaleW(1./self.wtot)

        file = open(name+".yoda", "a")

        for object in [self.thr, self.tzm, self.hjm, self.ljm, self.wjb,
                       self.njb, self.sph, self.cpr, self.dpr]:
            file.write("\n\n")
            file.write(str(object))

        file.close()
