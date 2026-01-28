import math as m
import random as r

from vector import Vec4
from particle import Particle, CheckEvent
from qcd import AlphaS, NC, TR, CA, CF

class Kernel:
    """
    This class represents the parton splitting functions.

    Attributes
    ----------
    flavs (arr) : The possible quark flavours
    """

    # Constructor
    def __init__(self, flavs):
        self.flavs = flavs

class Pqq (Kernel):
    '''
    This class represents the splitting function for the process q->qg

    Methods
    -------
    Value : Calculates the actual value of the splitting function
    Estimate : Overestimates the value of the splitting function
    Integral : Defines the integral of the overestimate
    GenerateZ : Generates a value for the kinematic quantity z
    '''

    def Value(self,z,y,t):
        '''
        This function calculates the value of the splitting function

        Parameters
        ----------
        z (float) : The energy fraction for the emission itself
        y (float) : The energy fraction for the recoil of partons
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The value of the splitting function
        '''
        return CF*(2./(1.-z*(1.-y))-(1.+z))
    
    def Estimate(self,z,t):
        '''
        This function overestimates the true value of the splitting function

        Parameters
        ----------
        z (float) : The energy fraction for the emission itself
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The overestimated value of the splitting function
        '''
        return CF*2./(1.-z)
    
    def Integral(self,zm,zp,t):
        '''
        This function defines the integral of the overestimate

        Parameters
        ----------
        zm (float) : The lower bound on z
        zp (float) : The upper bound on z
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The integral of the overestimate
        '''
        return CF*2.*m.log((1.-zm)/(1.-zp))
    
    def GenerateZ(self,zm,zp):
        '''
        This function generates a value of z between the upper and lower bounds

        Parameters
        ----------
        zm (float) : The lower bound on z
        zp (float) : The upper bound on z

        Returns
        -------
        (float) : The generated value of z        
        '''
        return 1.+(zp-1.)*m.pow((1.-zm)/(1.-zp),r.random())

class Pgg (Kernel):
    '''
    This class represents the splitting function for the process g->gg

    Methods
    -------
    Value : Calculates the actual value of the splitting function
    Estimate : Overestimates the value of the splitting function
    Integral : Defines the integral of the overestimate
    GenerateZ : Generates a value for the kinematic quantity z
    '''

    def Value(self,z,y,t):
        '''
        This function calculates the value of the splitting function

        Parameters
        ----------
        z (float) : The energy fraction for the emission itself
        y (float) : The energy fraction for the recoil of partons
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The value of the splitting function
        '''
        return CA/2.*(2./(1.-z*(1.-y))-2.+z*(1.-z))

    def Estimate(self,z,t):
        '''
        This function overestimates the true value of the splitting function

        Parameters
        ----------
        z (float) : The energy fraction for the emission itself
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The overestimated value of the splitting function
        '''
        return CA/(1.-z)

    def Integral(self,zm,zp,t):
        '''
        This function defines the integral of the overestimate

        Parameters
        ----------
        zm (float) : The lower bound on z
        zp (float) : The upper bound on z
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The integral of the overestimate
        '''
        return CA*m.log((1.-zm)/(1.-zp))

    def GenerateZ(self,zm,zp):
        '''
        This function generates a value of z between the upper and lower bounds

        Parameters
        ----------
        zm (float) : The lower bound on z
        zp (float) : The upper bound on z

        Returns
        -------
        (float) : The generated value of z        
        '''
        return 1.+(zp-1.)*m.pow((1.-zm)/(1.-zp),r.random())

class Pgq (Kernel):
    '''
    This class represents the splitting function for the process g->qqbar

    Methods
    -------
    Value : Calculates the actual value of the splitting function
    Estimate : Overestimates the value of the splitting function
    Integral : Defines the integral of the overestimate
    GenerateZ : Generates a value for the kinematic quantity z
    '''

    def Value(self,z,y,t):
        '''
        This function calculates the value of the splitting function

        Parameters
        ----------
        z (float) : The energy fraction for the emission itself
        y (float) : The energy fraction for the recoil of partons
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The value of the splitting function
        '''
        return TR/2.*(1.-2.*z*(1.-z))

    def Estimate(self,z,t):
        '''
        This function overestimates the true value of the splitting function

        Parameters
        ----------
        z (float) : The energy fraction for the emission itself
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The overestimated value of the splitting function
        '''
        return TR/2.

    def Integral(self,zm,zp,t):
        '''
        This function defines the integral of the overestimate

        Parameters
        ----------
        zm (float) : The lower bound on z
        zp (float) : The upper bound on z
        t (float) : The ordering variable for the process

        Returns
        -------
        (float) : The integral of the overestimate
        '''
        return TR/2.*(zp-zm)

    def GenerateZ(self,zm,zp):
        '''
        This function generates a value of z between the upper and lower bounds

        Parameters
        ----------
        zm (float) : The lower bound on z
        zp (float) : The upper bound on z

        Returns
        -------
        (float) : The generated value of z        
        '''
        return zm+(zp-zm)*r.random()
    
class Shower:
    '''
    This class represents the parton shower itself

    Attributes
    ----------
    t0 (float) : The cutoff scale of the ordering parameter
    alpha (float) : The strong running coupling
    alphamax (float) : The maximum value of the strong running coupling
    kernels (arr) : The list of possible parton splitting functions

    Methods
    -------
    MakeKinematics : Constructs the kinematics for a 2->3 parton branching
    MakeColours : Ensures that colour charge is conserved in the 2->3 parton branching
    GeneratePoint : Implements the veto algorithm for generating parton emissions
    Run : Responsible for running the shower
    Run_Return : Runs the shower while also returning all intermediate momenta and ordering variables
    '''

    # Constructor
    def __init__(self,alpha,t0):
        self.t0 = t0
        self.alpha = alpha
        self.alphamax = alpha(self.t0)
        self.kernels = [ Pqq([fl,fl,21]) for fl in [-5,-4,-3,-2,-1,1,2,3,4,5] ]
        self.kernels += [ Pgg([21,21,21]) ]
        self.kernels += [ Pgq([21,fl,-fl]) for fl in [1,2,3,4,5] ]

    def MakeKinematics(self,z,y,phi,pijt,pkt):
        '''
        This function constructs the kinematics of the outgoing partons in a 2->3 branching

        Parameters
        ----------
        z (float) : The energy fraction associated with the actual emission
        y (float) : The energy fraction associated with the parton recoil
        phi (float) : The azimuthal angle of the emission
        pijt (arr) : The momentum of the emitter before the emission
        pkt (arr) : The momentum of the spectator before the emission

        Returns
        -------
        pi (arr) : The momentum of the emitter after the emission
        pj (arr) : The momentum of the emitted parton
        pk (arr) : The momentum of the spectator after the emission
        '''

        # Define the centre of mass energy
        Q = pijt+pkt

        # Define the transverse momentum
        rkt = m.sqrt(Q.M2()*y*z*(1.-z))
        kt1 = pijt.Cross(pkt)
        if kt1.P() < 1.e-6: kt1 = pijt.Cross(Vec4(0.,1.,0.,0.))
        kt1 *= rkt*m.cos(phi)/kt1.P()
        kt2cms = Q.Boost(pijt).Cross(kt1)
        kt2cms *= rkt*m.sin(phi)/kt2cms.P()
        kt2 = Q.BoostBack(kt2cms)

        # Compute the outgoing kinematics
        pi = z*pijt + (1.-z)*y*pkt + kt1 + kt2
        pj = (1.-z)*pijt + z*y*pkt - kt1 - kt2
        pk = (1.-y)*pkt        
        return [pi,pj,pk]

    def MakeColors(self,flavs,colij,colk):
        '''
        This function ensures that colour charge is conserved in the 2->3 branching

        Parameters
        ----------
        colij (arr) : The colour and anticolour of the emitter before the emission
        colk (arr) : The colour and anticolour of the spectator before the emission

        Returns
        -------
        (arr) : The colour and anticolour of the emitter and emitted parton after the emission
        '''

        # Set an arbitrary new colour
        self.c += 1

        # Check emitter is not a gluon
        if flavs[0] != 21:
            # Check emitter is a quark or antiquark
            if flavs[0] > 0:
                return [ [self.c,0], [colij[0],self.c] ]
            else:
                return [ [0,self.c], [self.c,colij[1]] ]
        else:
            # Check emitted parton is a gluon
            if flavs[1] == 21:
                if colij[0] == colk[1]:
                    if colij[1] == colk[0] and r.random()>0.5:
                        return [ [colij[0],self.c], [self.c,colij[1]] ]
                    return [ [self.c,colij[1]], [colij[0],self.c] ]
                else:
                    return [ [colij[0],self.c], [self.c,colij[1]] ]
            else:
                if flavs[1] > 0:
                    return [ [colij[0],0], [0,colij[1]] ]
                else:
                    return [ [0,colij[1]], [colij[0],0] ]
                
    def GeneratePoint(self,event):
        '''
        This function is responsible for implementing the veto algorithm, and generating all of the
        emissions in the shower

        Parameters
        ----------
        event (arr) : The list of parton momenta in the event

        Returns
        -------
        moms (arr) : The list of momenta for three partons at each emission
        tlist (arr) : The list of ordering variables for which emissions occurred
        '''

        tlist = []

        # Loop until the hadronisation scale is reached
        while self.t > self.t0:
            t = self.t0
            # Iterate through all emitter/spectator pairs
            for split in event[2:]:
                for spect in event[2:]:
                    # Skip pairs where the emitter is the spectator
                    if spect == split: continue

                    # Skip pairs where there will be no colour interaction
                    if not split.ColorConnected(spect): continue

                    # Iterate through all splitting functions
                    for sf in self.kernels:
                        # Skip emitters that are not partons
                        if sf.flavs[0] != split.pid: continue
                        
                        # Skip emissions with an empty z range
                        m2 = (split.mom+spect.mom).M2()
                        if m2 < 4.*self.t0: continue

                        # Compute z boundaries and calculate ordering variable for emission
                        zp = .5*(1.+m.sqrt(1.-4.*self.t0/m2))
                        g = self.alphamax/(2.*m.pi)*sf.Integral(1.-zp,zp,t)
                        tt = self.t*m.pow(r.random(),1./g)

                        # Remember only the highest accepted ordering variable
                        if tt > t:
                            t = tt
                            s = [ split, spect, sf, m2, zp ]
            self.t = t
            
            # Check hadronisation scale has not been reached
            if t > self.t0:
                # Compute energy fractions
                z = s[2].GenerateZ(1.-s[4],s[4])
                y = t/s[3]/z/(1.-z)

                # Check that y is allowed
                if y < 1.:
                    # Define weighting of emission
                    f = (1.-y)*self.alpha(t)*s[2].Value(z,y,t)
                    g = self.alphamax*s[2].Estimate(z,t)

                    # Perform acceptance test
                    if f/g > r.random():
                        # Extract t value
                        tlist.append(t)

                        # Generate azimuthal angle at random
                        phi = 2.*m.pi*r.random()

                        # Construct kinematics and colour conservation
                        moms = self.MakeKinematics(z,y,phi,s[0].mom,s[1].mom)
                        cols = self.MakeColors(s[2].flavs,s[0].col,s[1].col)

                        # Add new parton to list of partons, and update momenta and colours
                        event.append(Particle(s[2].flavs[2],moms[1],cols[1]))
                        s[0].Set(s[2].flavs[1],moms[0],cols[0])
                        s[1].mom = moms[2]
                        return moms, tlist
                    
        # Assuming no emissions are possible, return zeros
        return (0,0,0,0), [0]
                    
    def Run(self,event,t):
        '''
        This function runs the full shower

        Parameters
        ----------
        event (arr) : The list of particles in the event
        t (float) : The ordering parameters
        '''

        # Set up the initial colour and ordering variable
        self.c = 1
        self.t = t

        # Generate emissions until hadronisation scale is reached
        while self.t > self.t0:
            self.GeneratePoint(event)
        
    def Run_Return(self,event,t):
        '''
        This function runs the full shower, and returns the all of the intermediate momenta and
        ordering variables for which emisisons occurred

        Parameters
        ----------
        event (arr) : The list of particles in the event
        t (float) : The ordering parameters

        Returns
        -------
        moms (arr) : The list of momenta for three partons at each emission
        tlist (arr) : The list of ordering variables for which emissions occurred
        '''

        # Set up the initial colour and ordering variable
        self.c = 1
        self.t = t

        # Generates emissions until hadronisation scale is reached, and returns intermediate
        # momenta and ordering variables
        while self.t > self.t0:
            moms, tlist = self.GeneratePoint(event)
            yield moms, tlist
            

def progress_bar(iteration, total, length=30):
    '''
    This function prints the current stage of progress and is used within other functions as a visual
    aid for how close to completion the program is

    Parameters
    ----------
    iteration (int) : A value representing the current stage of the program
    total (int) : The maximum value of iteration, at which point the program should be finished
    length (int) : Physical length of the progress bar string
    '''

    # Normalise current stage of program as a percentage
    percent = 100 * (iteration / float(total))

    # Calculate how many characters should be filled up in the string so far
    filled_length = int(length * iteration // total)

    # Print the current stage of progress
    bar = '█' * filled_length + '-' * (length - filled_length - 1)
    print(f'\r|{bar}| {iteration} Events', end='\r')
    if iteration == total:
        print()
