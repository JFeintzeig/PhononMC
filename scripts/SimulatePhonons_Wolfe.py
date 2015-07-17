import numpy, os
from lattice import Shape, Volume, Crystal, TeO2
from phonons import Mode, ProcessGenerator, Propagator, Controller
from base import *

OutDir='/Users/jfeintzeig/CUORE/Scripts/2015/MCInteractions/plots/Phonons/'
os.system('mkdir -p %s' % (OutDir))

# initial conditions
mode=Mode.L
x=numpy.ones(3)*2.5*cm
w=1e11*Hz #TODO: or nu????
k=numpy.array([0,0])
t=0 
print "Pos:", x
print "k:", k
print "w:", w

# crystal is 10cm sphere, detector is 5cm sphere inside crystal
CrystalShape=Volume(Shape.CUBE,5.*cm,numpy.array([2.5,2.5,2.5]))
DetectorShape=Volume(Shape.CUBE,2*0.145*cm,numpy.array([2.5,2.5,2.5]))

C_iso=TeO2.C_iso
C_dec=TeO2.C_dec

TeO2.C_dec=5.5e-53
#TeO2.C_dec=5.5e-100
#TeO2.C_iso=5.5e-100

names=['Anisotropic']
items=[TeO2]

data={}
for name,item in zip(names,items):
    print name
    TeO2Sphere=Crystal(CrystalShape,DetectorShape,TeO2)
    Process=ProcessGenerator(TeO2)
    Prop=Propagator(item,orient=numpy.radians(45))
    MissionControl=Controller(TeO2Sphere,Process,Prop,randomize_k=True,
        randomize_mode=True)

    # create container
    MissionControl.DefineInitialPhonon(mode,x,k,w,t)

    n=int(5e5)

    print "n:", n
    import time
    t_start=time.time()
    MissionControl.Execute(n)
    t_end=time.time()
    data[name]=MissionControl.GetOutput()
    print (t_end-t_start)/n

import pylab
import matplotlib.gridspec as gridspec
bins=numpy.linspace(0,1.5,51)
cols={0:'black',1:'red',2:'blue'}
ls={'Isotropic':'dashed','Anisotropic':'solid'}
#for disp in [0.0,0.05,0.06,0.07,0.09,0.11,0.12]: # scan z
i=6
fig=pylab.figure(figsize=(8,10))
gs1 = gridspec.GridSpec(7,1)
gs1.update(wspace=0,hspace=0)
for disp in [0.06,0.07,0.08,0.09,0.10,0.14,0.15]: # scan y
    ax=pylab.subplot(gs1[i])
    pylab.axis('on')
    for name in data.keys():
        for a in [0,1,2]:
            ind=(data[name].mode==a)&(data[name].x0==2.645)\
                &(data[name].x1<(2.5-disp+0.01))&(data[name].x1>(2.5-disp-0.01))\
                &(data[name].x2<2.51)&(data[name].x2>2.49)
            if numpy.sum(ind)<2:
                continue
            pylab.hist(data[name].t[ind]*1e6,bins,histtype='step',
                log=False,linewidth=2,linestyle=ls[name],
                label=Mode(a).name,color=cols[a])
            pylab.figtext(0.15,0.08+0.114*(7-i),'Offset = %3.2f cm' % (disp))

        pylab.xlim(0,1.5)
        if i==6:
            pylab.xlabel(r'Time ($\mu$s)',size='16')
            pylab.legend()
        else:
            ax.set_xticklabels([])
        pylab.ylabel('Counts/bin')
        pylab.grid()
    i-=1

pylab.title('Phonon time distribution for y-axis offsets')
pylab.savefig(OutDir+'Wolfe_Time.png')
pylab.show()

#f=open('wolfe_3e5.npy','w')
#numpy.save(f,data)
#f.close()
