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

# kill decay

data={}
dec_rates=numpy.logspace(-50,-54,5)
for dr in dec_rates:
    print dr
    TeO2.C_dec=dr
    TeO2Sphere=Crystal(CrystalShape,DetectorShape,TeO2)
    Process=ProcessGenerator(TeO2)
    Prop=Propagator()
    MissionControl=Controller(TeO2Sphere,Process,Prop)

    # create container
    MissionControl.DefineInitialPhonon(mode,x,k,w,t)

    n=1000

    print "n:", n
    import time
    t_start=time.time()
    MissionControl.Execute(n)
    t_end=time.time()
    data[dr]=MissionControl.GetOutput()
    print (t_end-t_start)/n

import pylab
bins=numpy.linspace(0,3e2,101)
for item in numpy.sort(data.keys()):
    pylab.hist(data[item].t*1e6,bins,histtype='step',log=True,normed=True,
        linewidth=2,label=r'log$_{10}$($C_{dec}$) = %3.2f' % (numpy.log10(item)))

# ballistic time
t_ball=DetectorShape.size/(4200.*m/s)
pylab.axvline(t_ball*1e6,linestyle='dashed',linewidth=4,
    color='red',label='Ballistic time')

pylab.legend(loc='lower right')
pylab.xlabel(r'Time ($\mu$s)')
pylab.ylabel('Counts/bin')
pylab.grid()
pylab.title('1 THz Phonon time distribution for different decay rates')
pylab.savefig(OutDir+'TimeVsC_dec.png')
pylab.show()
