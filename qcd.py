import math as m

NC = 3.
TR = 1./2.
CA = NC
CF = (NC*NC-1.)/(2.*NC)

class AlphaS:

    def __init__(self,mz,asmz,order=1,cmw=False,mb=4.75,mc=1.3):
        self.order = order
        self.cmw = cmw
        self.mc2 = mc*mc
        self.mb2 = mb*mb
        self.mz2 = mz*mz
        self.asmz = asmz
        self.asmb = self(self.mb2)
        self.asmc = self(self.mc2)
        print("\\alpha_s({0}) = {1}".format(mz,self(self.mz2)))

    def beta0(self,nf):
        return 11./6.*CA-2./3.*TR*nf

    def beta1(self,nf):
        return 17./6.*CA*CA-(5./3.*CA+CF)*TR*nf
    
    def kcmw(self,nf):
        '''
        This function computes the K_CMW constant.

        Parameters
        ----------
        nf (int) : The number of quark flavours active at the given energy scale

        Returns
        -------
        (float) : The value of the K_CMW constant
        '''
        return 1./(2.*m.pi)*((67./18.-m.pi**2/6.)*CA-5./9.*nf)

    def as0(self,t):
        if t >= self.mb2:
            tref = self.mz2
            asref = self.asmz
            b0 = self.beta0(5)/(2.*m.pi)
        elif t >= self.mc2:
            tref = self.mb2
            asref = self.asmb
            b0 = self.beta0(4)/(2.*m.pi)
        else:
            tref = self.mc2
            asref = self.asmc
            b0 = self.beta0(3)/(2.*m.pi)
        return 1./(1./asref+b0*m.log(t/tref))

    def as1(self,t):
        if t >= self.mb2:
            tref = self.mz2
            asref = self.asmz
            b0 = self.beta0(5)/(2.*m.pi)
            b1 = self.beta1(5)/pow(2.*m.pi,2)
        elif t >= self.mc2:
            tref = self.mb2
            asref = self.asmb
            b0 = self.beta0(4)/(2.*m.pi)
            b1 = self.beta1(4)/pow(2.*m.pi,2)
        else:
            tref = self.mc2
            asref = self.asmc
            b0 = self.beta0(3)/(2.*m.pi)
            b1 = self.beta1(3)/pow(2.*m.pi,2)
        w = 1.+b0*asref*m.log(t/tref)
        return asref/w*(1.-b1/b0*asref*m.log(w)/w)
    
    def ascmw(self, t):
        '''
        This function calculates the CMW correction to the one-loop strong running coupling

        Parameters
        ----------
        t (float) : The ordering variable

        Returns
        -------
        (float) : The value of the CMW correction to the one-loop strong running coupling
        '''

        # Calculate the one-loop strong running coupling
        alpha_s_t = self.as1(t)

        # Check the energy scale against the quark masses
        if t >= self.mb2: nf = 5
        elif t >= self.mc2: nf = 4
        else: nf = 3

        # Calculate the K_CMW constant
        kcmw = self.kcmw(nf)
        
        return alpha_s_t * (1 + kcmw * alpha_s_t)

    def __call__(self, t):
        if self.order == 1:
            if self.cmw: return self.ascmw(t)
            return self.as1(t)
        return self.as0(t)
