import numpy, os
from lattice import Shape, Volume, Crystal, TeO2
from phonons import Mode, ProcessGenerator, Propagator, Controller
from base import *

OutDir='/Users/jfeintzeig/CUORE/Scripts/2015/MCInteractions/plots/Phonons/'
os.system('mkdir -p %s' % (OutDir))

# initial conditions
mode=Mode.L
x=numpy.ones(3)*0.*cm
w=1e12*Hz #TODO: or nu????
k=numpy.array([0,0])
t=0 
print "Pos:", x
print "k:", k
print "w:", w

# crystal is 10cm sphere, detector is 5cm sphere inside crystal
CrystalShape=Volume(Shape.SPHERE,10.*cm,numpy.array([0,0,0]))
DetectorShape=Volume(Shape.SPHERE,5.*cm,numpy.array([0,0,0]))

C_iso=TeO2.C_iso
C_dec=TeO2.C_dec

TeO2.C_dec=5.5e-53

# kill decay

for r in [0.1,0.5,1.,2.5,5.]:
    print '%e' % w
    DetectorShape=Volume(Shape.SPHERE,r*cm,numpy.array([0,0,0]))
    TeO2Sphere=Crystal(CrystalShape,DetectorShape,TeO2)
    Process=ProcessGenerator(TeO2)
    Prop=Propagator(TeO2)
    MissionControl=Controller(TeO2Sphere,Process,Prop)

    # create container
    MissionControl.DefineInitialPhonon(mode,x,k,w,t)

    n=5000

    print "n:", n
    import time
    t_start=time.time()
    MissionControl.Execute(n)
    t_end=time.time()
    data=MissionControl.GetOutput()
    print (t_end-t_start)/n

import pylab
bins=numpy.linspace(0,3e2,101)
for item in data.keys():
    pylab.hist(data[item].t*1e6,bins,histtype='step',log=True,normed=True,
        linewidth=2,label='r = %3.2f cm' % (item))

pylab.legend(loc='lower right')
pylab.xlabel(r'Time ($\mu$s)')
pylab.ylabel('Counts/bin')
pylab.grid()
pylab.title('Phonon time distribution at different radii')
pylab.savefig(OutDir+'TimeVsRadius.png')
pylab.show()

bins=numpy.linspace(0,3e2,101)
for item in data.keys():
    pylab.hist(data[item].w/1e9,bins,histtype='step',log=False,normed=True,
        linewidth=2,label='r = %3.2f cm' % (item))

pylab.legend(loc='upper right')
pylab.xlabel(r'$\nu$ (GHz)')
pylab.ylabel('Counts/bin')
pylab.grid()
pylab.title('Phonon freq distribution at different radii')
pylab.savefig(OutDir+'FreqVsRadius.png')
pylab.show()
