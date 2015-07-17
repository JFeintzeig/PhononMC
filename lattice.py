import numpy
from base import *

norm=numpy.linalg.norm

#material
#TODO: get values for TeO2!!!
#TODO: in this file i state all material constants and create TeO2 intance of Lattice
# i could put this stuff in a different file, away from the class definitions
C_dec=1e-50 # totally random just to test decays
C_iso=5.0*10**(-41)*s**3
MODE_POPULATIONS=(0.1,0.35,0.55) # (L,FT,ST)

# from Tamura Sangu Maris PHYSICAL REVIEW B 68, 014302 (2003)
C11=55.9*DC
C33=105.5*DC
C44=26.6*DC
C66=66.3*DC
C12=51.6*DC
C13=23.9*DC
C111=-160.*DC
C112=-600.*DC
C113=-140.*DC
C123=-110.*DC
C133=180.*DC
C144=-41.*DC
C155=36.*DC
C166=-640.*DC
C333=-2110.*DC
C344=-54.*DC
C366=-260.*DC
C456=-250.*DC

density = 5.99 # g/cm^3

# C squiggle for matrix, from Tamura Sangu Maris
SC12=C12+C66
SC13=C13+C44

# inelasticity tensor, pg 76. of Symmetry, Group Theory,
# and the Physical Properties of Crystals by Powell, for D4
elasticity=numpy.zeros((6,6))
elasticity[0,0]=C11
elasticity[1,1]=C11
elasticity[2,2]=C33
elasticity[3,3]=C44
elasticity[4,4]=C44
elasticity[5,5]=C66
elasticity[0,1]=C12
elasticity[1,0]=C12
elasticity[0,2]=C13
elasticity[1,2]=C13
elasticity[2,0]=C13
elasticity[2,1]=C13
elasticity=numpy.matrix(elasticity)

from phonons import Mode

#TODO: move this to phonons.py, reformulate with self.'s, etc, accordingly
def KToV(k,mode):
    # from wolfe appendix 2, eqn 5, 10, 11, 12, but for D4 crystal
    evals,evecs=SolveM(k)
    if mode==Mode.ST:
        ind=0
    elif mode==Mode.FT:
        ind=1
    elif mode==Mode.L:
        ind=2
    else:
        sys.exit('unknown mode')
    w=numpy.sqrt(evals[ind]/density)
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
    k=numpy.matrix(k)
    ke=k*e
    v_g=e*(elasticity*ke.T)*1/(density*w)
    return v_g

def SolveM(k):
    # from Tamura Sangu Maris PHYSICAL REVIEW B 68, 014302 (2003)
    # if i wanted to clean this up i could figure out how to just calculate M
    # from elasticity and k
    kx=k[0]
    ky=k[1]
    kz=k[2]
    M=numpy.matrix(numpy.zeros((3,3)))
    M[0,0]=C11*kx*kx+C66*ky*ky+C44*kz*kz
    M[0,1]=SC12*kx*ky
    M[0,2]=SC13*kz*kx
    M[1,0]=SC12*kx*ky
    M[1,1]=C66*kx*kx+C11*ky*ky+C44*kz*kz
    M[1,2]=SC13*ky*kz
    M[2,0]=SC13*kz*kx
    M[2,1]=SC13*ky*kz
    M[2,2]=C44*(kx*kx+ky*ky)+C33*kz*kz
    eigenvalues,eigenvectors=numpy.linalg.eig(M)
    return eigenvalues,eigenvectors

# phase velocities, should be 2.11 and 4.20 km/s
# first is ST, then FT, then L
evals,evecs=SolveM([0,0,1])
print numpy.sqrt(evals/density)

# volumes
class Shape(Enum):
    SPHERE=0
    CUBE=1

class Volume:
    def __init__(self,shape,size,origin):
        self.shape=shape
        self.size=size
        self.origin=origin

class Lattice:
    def __init__(self,C_iso,C_dec,modes,density,elasticity):
        self.C_iso=C_iso
        self.C_dec=C_dec
        self.mode_populations=modes
        self.density=density
        self.elasticity=elasticity

class Crystal:
    def __init__(self,volume,detector_volume,lattice):
        self.volume=volume
        self.detector_volume=detector_volume
        self.lattice=lattice

    def Outside(self,volume,pos):
        if volume.shape==Shape.SPHERE:
            if numpy.sqrt(sum((pos-volume.origin)**2))>=volume.size:
                return True
            else:
                return False
        if volume.shape==Shape.CUBE:
            if max(abs(pos-volume.origin))>=volume.size/2.:
                return True
            else:
                return False

    def OutsideCrystal(self,pos):
        return self.Outside(self.volume,pos)

    def Detected(self,pos):
        # TODO: for a spherical detection surface surrounding
        # initial phonon, i want to know when it crosses outside
        # but for detecting on a small square on edge of surface,
        # that's really more of an Inside() than Outside() sort of
        # problem...I'll need to edit this to do that...
        return self.Outside(self.detector_volume,pos)

    def SurfaceCrossingTime(self,volume,pos,v):
        if volume.shape==Shape.SPHERE:
            max_coord=volume.origin+volume.size
            min_coord=volume.origin-volume.size
            #http://en.wikipedia.org/wiki/Line-sphere_intersection
            l=v/norm(v)
            a=numpy.dot(l,pos-volume.origin)
            b=numpy.sqrt(a**2-norm(pos-volume.origin)**2+volume.size**2)
            t1=(-a+b)/norm(v)
            t2=(-a-b)/norm(v)
            t_surf=max(t1,t2)
            assert(t_surf>0)
        elif volume.shape==Shape.CUBE:
            max_coord=volume.origin+volume.size/2.
            min_coord=volume.origin-volume.size/2.
            # x=vt+x0, so t=(x-x0)/v
            t_max=(max_coord-pos)/v
            t_min=(min_coord-pos)/v
            # find smallest positive time (neg times are in opposite dir)
            t_pos=numpy.concatenate((t_min,t_max))
            t_surf=min(t_pos[t_pos>=0])
        else:
            sys.exit('what kinda volume is this!?')
        return t_surf

    def DetectorCrossingTime(self,pos,v):
        return self.SurfaceCrossingTime(self.detector_volume,pos,v)

TeO2=Lattice(C_iso,C_dec,MODE_POPULATIONS,density,elasticity)
