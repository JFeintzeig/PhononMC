import numpy, sys
from enum import Enum

#material
from base import *

#TODO: error checking!

rand=numpy.random.rand

class Proc(Enum):
    SCATTERING=0
    DECAY=1
    REFLECT=2
    DETECTED=3

class Mode(Enum):
    L=0 # longitudinal
    FT=1 # fast transverse
    ST=2 # slow transverse

def SphToCart(a):
    x=numpy.sin(a[0])*numpy.cos(a[1])
    y=numpy.sin(a[0])*numpy.sin(a[1])
    z=numpy.cos(a[0])
    return numpy.array([x,y,z])

class Phonon:
    def __init__(self,mode,x,k,w,t):
        self.mode=mode
        self.x=x
        self.k=k
        self.w=w
        self.t=t

class Phonons:
    '''
    container class to hold variables describing a group of phonons:
            -mode
            -x position
            -k wavevector
            -w frequency
            -t time
    '''
    def __init__(self):
        self.count=0
        self.count=0
        self.mode={}
        self.x={}
        self.k={}
        self.w={}
        self.t={}

    def Get(self,count=-1):
        if count==-1:
            '''return phonon with smallest time'''
            num=min(self.t,key=self.t.get)
        else:
            num=self.count
        data=[]
        for item in [self.mode,self.x,self.k,self.w,self.t]:
            data+=[item.pop(num)]
        return data

    def Add(self,mode,x,k,w,t):
        self.mode.update({self.count:mode})
        self.x.update({self.count:x})
        self.k.update({self.count:k})
        self.w.update({self.count:w})
        self.t.update({self.count:t})
        self.count+=1

class Propagator:
    def __init__(self,lattice=None,orient=0.):
        self.lattice=lattice
        #TODO: make sure orient angle units are correct
        self.orient=orient
        if lattice!=None:
            self.C11=lattice.elasticity[0,0]
            self.C33=lattice.elasticity[2,2]
            self.C44=lattice.elasticity[3,3]
            self.C66=lattice.elasticity[5,5]
            self.C12=lattice.elasticity[0,1]
            self.C13=lattice.elasticity[0,2]
            self.SC12=self.C12+self.C66
            self.SC13=self.C13+self.C44
        #orientation of the x-axis in terms of crystal structure
        #k-to-v mapping coords are according to crystal structure
        #but we need k and v according to physical orientation
        # so rotate v back so that [1,0,0] is the orientation vector
        # of the crystal

    def Rotate(self,v):
        # orient lets the user set the crystal orientation along the [1,0,0]
        # x-axis of th simulated crystal.  orient is the angle between [1,0,0] and
        # the desired crystal orientation vector for the simulated crystal
        # ie, if i want crystal direction [110] to be the x-axis of the crystal,
        # i set orient to 45 deg (but in rad)...similarly, if i want phonons
        # propagating along [010] to be moving down the x-axis of the simulated
        # crystal, i set orient to 90 deg
        #TODO: only x-rotations (around z axis) supported right now
        c=numpy.cos(-self.orient)
        s=numpy.sin(-self.orient)
        rot_mat=numpy.matrix([c,-s,0,s,c,0,0,0,1]).reshape(3,3)
        v=numpy.matrix(v).T
        v_det=(rot_mat*v).T
        return v_det.view(numpy.ndarray).flatten()

    def MapKToV(self,k,mode=None):
        if self.lattice==None:
            return self.MapKToV_Isotropic(k)
        else:
            return self.KToV(k,mode)

    def MapKToV_Isotropic(self,k):
        return 4200.*m/s*SphToCart(k)

    def KToV(self,k,mode):
        if mode==None:
            sys.exit('gimme a mode!!')
        # from wolfe appendix 2, eqn 5, 10, 11, 12, but for D4 crystal
        evals,evecs=self.SolveM(k)
        # choose mode by velocity magnitude order
        #TODO: is this really correct!?!
        inds=numpy.argsort(evals)
        evals=evals[inds]
        evecs=evecs[:,inds]
        if mode==Mode.ST:
            ind=0
        elif mode==Mode.FT:
            ind=1
        elif mode==Mode.L:
            ind=2
        else:
            sys.exit('unknown mode')
        w=numpy.sqrt(evals[ind]/self.lattice.density)
        ex,ey,ez=evecs[:,ind]
        e=numpy.zeros((3,6))
        e[0,0]=ex
        e[0,4]=ez
        e[0,5]=ey
        e[1,1]=ey
        e[1,3]=ez
        e[1,5]=ex
        e[2,2]=ez
        e[2,3]=ey
        e[2,4]=ex
        e=numpy.matrix(e)
        k=numpy.matrix(SphToCart(k))
        ke=k*e
        v_g=e*(self.lattice.elasticity*ke.T)*1/(self.lattice.density*w)
        return self.Rotate(v_g.view(numpy.ndarray).flatten())

    def SolveM(self,k):
        # from Tamura Sangu Maris PHYSICAL REVIEW B 68, 014302 (2003)
        # if i wanted to clean this up i could figure out how to just calculate M
        # from elasticity and k
        k=SphToCart(k)
        kx=k[0]
        ky=k[1]
        kz=k[2]
        M=numpy.matrix(numpy.zeros((3,3)))
        M[0,0]=self.C11*kx*kx+self.C66*ky*ky+self.C44*kz*kz
        M[0,1]=self.SC12*kx*ky
        M[0,2]=self.SC13*kz*kx
        M[1,0]=self.SC12*kx*ky
        M[1,1]=self.C66*kx*kx+self.C11*ky*ky+self.C44*kz*kz
        M[1,2]=self.SC13*ky*kz
        M[2,0]=self.SC13*kz*kx
        M[2,1]=self.SC13*ky*kz
        M[2,2]=self.C44*(kx*kx+ky*ky)+self.C33*kz*kz
        eigenvalues,eigenvectors=numpy.linalg.eig(M)
        return eigenvalues,eigenvectors

class ProcessGenerator:
    def __init__(self,lattice):
        self.C_iso=lattice.C_iso
        self.C_dec=lattice.C_dec
        # mode_populations is like (0.15,0.6,0.25)
        # mode_mixing should be [0.15,0.75]
        self.mode_mixing=[lattice.mode_populations[0],lattice.mode_populations[0]+\
            lattice.mode_populations[1]]
        self.lattice=lattice

    def GetTimeStep(self,w,mode):
        #TODO: convert w to nu???
        tau_iso=1./(self.C_iso*w**4)
        tau_dec=1./(self.C_dec*w**5)
        if mode==Mode.L:
            tau=1/(1/tau_iso+1/tau_dec)
        else:
            tau=tau_iso
        timestep=-1.*tau*numpy.log(rand())
        if mode==Mode.L and rand()<(1/tau_dec)/(1/tau_iso+1/tau_dec):
            interaction_type=Proc.DECAY
        else:
            interaction_type=Proc.SCATTERING
        return timestep,interaction_type

    def Scatter(self):
        # need to implement non-isotropic scattering
        u=rand()
        # choose new mode
        if u<self.mode_mixing[0]:
            mode=Mode.L
        elif u>self.mode_mixing[1]:
            mode=Mode.ST
        else:
            mode=Mode.FT
        # get random direction
        d=rand(2)
        d[0]=numpy.arccos(2*d[0]-1)
        d[1]*=2*numpy.pi
        return mode, d

    def Decay(self,w):
        u=-1
        while u>1 or u<0:
            # approximate as a one-branch process with gaussian energy split
            # approximate that both phonons end with same k as initial phonon
            # and that both are FT modes
            #TODO: make this realistic, what is the real split between L/FT/ST?
            u=numpy.random.normal(0.5,0.21)
        mode1=Mode.FT
        mode2=Mode.FT
        if rand()<0.2:
            mode2=Mode.L
        return mode1, u*w, mode2, (1.-u)*w

class Controller:
    #TODO: implement other stopping condition?
    #TODO: implement randomizing initial phonon direction?
    def __init__(self,crystal,process_generator,propagator,
            randomize_k=False,randomize_mode=False):
        self.Crystal=crystal
        self.mode_mixing=[crystal.lattice.mode_populations[0],\
            crystal.lattice.mode_populations[0]+crystal.lattice.mode_populations[1]]
        self.Process=process_generator
        self.Prop=propagator
        self.Tracking=Phonons()
        self.Detected=Phonons()
        self.RandomizeK=randomize_k
        self.RandomizeMode=randomize_mode

    def DefineInitialPhonon(self,mode,x,k,w,t):
        self.InitialPhonon=Phonon(mode,x,k,w,t)

    def Execute(self,repeat=1):
        n_complete=0
        while n_complete<repeat:
            ip=self.InitialPhonon
            if self.RandomizeK:
                ip.k=numpy.array([numpy.arccos(2*rand()-1),2*numpy.pi*rand()])
            if self.RandomizeMode:
                u=rand()
                if u<self.mode_mixing[0]:
                    ip.mode=Mode.L
                elif u>self.mode_mixing[1]:
                    ip.mode=Mode.ST
                else:
                    ip.mode=Mode.FT
            self.Tracking.Add(ip.mode,ip.x,ip.k,ip.w,ip.t)
            # stop when all phonons have been tracked to completion...
            while len(self.Tracking.x)>0:
                mode,x,k,w,t=self.Tracking.Get()
                # get velocity
                v=self.Prop.MapKToV(k,mode)
                # get time for next process
                t_step,proc=self.Process.GetTimeStep(w,mode)

                if self.Crystal.Detected(x+v*t_step):
                    # track to edge, then save phonon as detected
                    t_surf=self.Crystal.DetectorCrossingTime(x,v)
                    x=x+v*t_surf
                    t+=t_surf
                    proc=Proc.DETECTED
                    self.Detected.Add(mode,x,k,w,t)

                elif self.Crystal.OutsideCrystal(x):
                    #TODO: implement reflection
                    proc==Proc.REFLECT
                    print "reflect"
                    sys.exit('not implemented yet')

                else:
                    # move to that location, apply process
                    x=x+v*t_step
                    t+=t_step
                    if proc==Proc.SCATTERING:
                        mode,k=self.Process.Scatter()
                        self.Tracking.Add(mode,x,k,w,t)
                    elif proc==Proc.DECAY:
                        mode1,w1,mode2,w2=self.Process.Decay(w)
                        self.Tracking.Add(mode1,x,k,w1,t)
                        self.Tracking.Add(mode2,x,k,w2,t)
                    else:
                        print "unknown process with id:", proc
                        sys.exit('not implemented yet')
            n_complete+=1
            if n_complete%1000==0:
                print "Completed:", n_complete
        return True

    def GetOutput(self):
        dtype=[]
        for item in ['mode','x0','x1','x2','k0','k1','w','t']:
            if item=='mode':
                dtype+=[(item,int)]
            else:
                dtype+=[(item,float)]
        data=numpy.zeros(self.Detected.count,dtype=dtype)
        for num in range(self.Detected.count):
            data[num]=(int(self.Detected.mode[num].value),
            self.Detected.x[num][0],
            self.Detected.x[num][1],
            self.Detected.x[num][2],
            self.Detected.k[num][0],
            self.Detected.k[num][1],
            self.Detected.w[num],
            self.Detected.t[num])

        data=data.view(numpy.recarray)
        return data
