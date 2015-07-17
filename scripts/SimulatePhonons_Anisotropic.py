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

names=['Isotropic','Anisotropic']
items=[None,TeO2]

data={}
for name,item in zip(names,items):
    print name
    TeO2Sphere=Crystal(CrystalShape,DetectorShape,TeO2)
    Process=ProcessGenerator(TeO2)
    Prop=Propagator(item)
    MissionControl=Controller(TeO2Sphere,Process,Prop)

    # create container
    MissionControl.DefineInitialPhonon(mode,x,k,w,t)

    n=5000

    print "n:", n
    import time
    t_start=time.time()
    MissionControl.Execute(n)
    t_end=time.time()
    data[name]=MissionControl.GetOutput()
    print (t_end-t_start)/n

import pylab
bins=numpy.linspace(0,2e2,51)
cols={0:'black',1:'red',2:'blue'}
ls={'Isotropic':'dashed','Anisotropic':'solid'}
for name in data.keys():
    for a in [0,1,2]:
        pylab.hist(data[name].t[data[name].mode==a]*1e6,bins,histtype='step',
            log=True,linewidth=2,label='%s, Mode=%i' % (name,a),color=cols[a],
            linestyle=ls[name])

pylab.legend()
pylab.xlabel(r'Time ($\mu$s)')
pylab.ylabel('Counts/bin')
pylab.grid()
pylab.title('Phonon time distribution')
pylab.savefig(OutDir+'Anisotropic_Time.png')
pylab.show()
