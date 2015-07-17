# tests to verify parts of simulation
import numpy, pylab
from phonons import ProcessGenerator, SphToCart
from base import *
import lattice

#Process=ProcessGenerator(lattice.C_iso,lattice.C_dec,
#    lattice.MODE_POPULATIONS)

test=4

# test scattering direction is isotropic
if test==0:
    x=[]
    y=[]
    z=[]

    for i in range(10000):
        a=SphToCart(Process.Scatter()[1])
        x+=[a[0]]
        y+=[a[1]]
        z+=[a[2]]
        i+=1

    pylab.hist(x,30,histtype='step',label='x')
    pylab.hist(y,30,histtype='step',label='y')
    pylab.hist(z,30,histtype='step',label='z')
    pylab.title('Random directions (cartesian) from scattering')
    pylab.legend(loc='lower right')
    pylab.show()

# test time step
if test==1:
    #TODO: implement time step test
    print "Time step not implemented yet!!!"

# test decay energy distribution

if test==2:
    print "Decay energy distribution test"
    def pv(v0,v):
        C=30/v0**4 # normed
        return C*v**2*(v0-v)**2

    u=numpy.linspace(0,1,100)

    for v0 in [1e11,5e11,1e12,5e12]:
        pylab.plot(u,pv(v0,u*v0),label=v0)
        print v0, pv(v0,0.5*v0), pv(v0,0.25*v0),pv(v0,0.85*v0)

    pylab.plot(u,stats.norm.pdf(u,0.5,0.21),label='norm')
    #TODO: decay energy test: plot results from actually sampling distribution...

    pylab.legend()
    pylab.show()

if test==3:
    #TODO: phonon position at detector test
    print "Phonon position at detector test"
    # plot x-coords of phonons when they are detected
    # for volumes with different centers and shapes
    #for item in Detected.x:
    #    print norm(Detected.x[item]-crystal.detector_volume.origin)

if test==4:
    #TODO: slowness surface test
    print "Slowness surface test"
    from lattice import KToV, SolveM
    from phonons import SphToCart
    from lattice import density
    theta=numpy.linspace(0,numpy.pi,100)
    phi=numpy.linspace(0,2*numpy.pi,100)
    # form array of isotropic k vectors
    k=[]
    for a in theta:
        for b in phi:
            k+=[SphToCart([a,b])]
    #for each, find phase velocity of each mode
    inv_vp_L=[]
    inv_vp_FT=[]
    inv_vp_ST=[]
    for kvec in k:
        evals,evecs=SolveM(kvec)
        dum=numpy.sqrt(evals/density)/1e5
        inv_vp_ST+=[1./dum[0]]
        inv_vp_FT+=[1./dum[1]]
        inv_vp_L+=[1./dum[2]]

    inv_vp_ST=numpy.array(inv_vp_ST)
    inv_vp_FT=numpy.array(inv_vp_FT)
    inv_vp_L=numpy.array(inv_vp_L)
    kx=numpy.array([item[0] for item in k])
    ky=numpy.array([item[1] for item in k])
    kz=numpy.array([item[2] for item in k])
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(inv_vp_L*kx,inv_vp_L*ky,inv_vp_L*kz)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    pylab.title('L')
    pylab.show()

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(inv_vp_ST*kx,inv_vp_ST*ky,inv_vp_ST*kz)
    pylab.title('ST')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    pylab.show()

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(inv_vp_FT*kx,inv_vp_FT*ky,inv_vp_FT*kz)
    pylab.title('FT')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    pylab.show()

    x=numpy.cos(phi)
    y=numpy.sin(phi)

    slow_xy_ST=[]
    slow_xy_FT=[]
    slow_xy_L=[]
    eptk_ST=[]
    eptk_FT=[]
    eptk_L=[]
    for a,b in zip(x,y):
        evals,evecs=SolveM([a,b,0])
        dum=numpy.sqrt(evals/density)/1e5
        dum=numpy.sort(dum)
        slow_xy_ST+=[1./dum[0]]
        slow_xy_FT+=[1./dum[1]]
        slow_xy_L+=[1./dum[2]]
        eptk_ST+=[numpy.dot([a,b,0],numpy.array(evecs[:,0]).flatten())]
        eptk_FT+=[numpy.dot([a,b,0],numpy.array(evecs[:,1]).flatten())]
        eptk_L+=[numpy.dot([a,b,0],numpy.array(evecs[:,2]).flatten())]
        
    slow_xy_ST=numpy.array(slow_xy_ST)
    slow_xy_FT=numpy.array(slow_xy_FT)
    slow_xy_L=numpy.array(slow_xy_L)
    eptk_ST=numpy.array(eptk_ST)
    eptk_FT=numpy.array(eptk_FT)
    eptk_L=numpy.array(eptk_L)
    fig=pylab.figure(figsize=(6,12))
    pylab.subplot(211)
    pylab.plot(slow_xy_ST*x,slow_xy_ST*y,label='ST')
    pylab.plot(slow_xy_FT*x,slow_xy_FT*y,label='FT')
    pylab.plot(slow_xy_L*x,slow_xy_L*y,label='L')
    pylab.legend()

    pylab.subplot(212)
    pylab.plot(phi,eptk_ST,label='ST')
    pylab.plot(phi,eptk_FT,label='FT')
    pylab.plot(phi,eptk_L,label='L')
    pylab.legend()
    pylab.show()

    slow_xy_ST=[]
    slow_xy_FT=[]
    slow_xy_L=[]
    eptk_ST=[]
    eptk_FT=[]
    eptk_L=[]
    for a,b in zip(x,y):
        evals,evecs=SolveM([0,a,b])
        dum=numpy.sqrt(evals/density)/1e5
        slow_xy_ST+=[1./dum[0]]
        slow_xy_FT+=[1./dum[1]]
        slow_xy_L+=[1./dum[2]]
        eptk_ST+=[numpy.dot([0,a,b],numpy.array(evecs[:,0]).flatten())]
        eptk_FT+=[numpy.dot([0,a,b],numpy.array(evecs[:,1]).flatten())]
        eptk_L+=[numpy.dot([0,a,b],numpy.array(evecs[:,2]).flatten())]
            
    slow_xy_ST=numpy.array(slow_xy_ST)
    slow_xy_FT=numpy.array(slow_xy_FT)
    slow_xy_L=numpy.array(slow_xy_L)
    eptk_ST=numpy.array(eptk_ST)
    eptk_FT=numpy.array(eptk_FT)
    eptk_L=numpy.array(eptk_L)

    fig=pylab.figure(figsize=(6,12))
    pylab.subplot(211)
    pylab.plot(slow_xy_ST*x,slow_xy_ST*y,label='ST')
    pylab.plot(slow_xy_FT*x,slow_xy_FT*y,label='FT')
    pylab.plot(slow_xy_L*x,slow_xy_L*y,label='L')
    pylab.legend()
    pylab.subplot(212)
    pylab.plot(eptk_ST*x,eptk_ST*y,label='ST')
    pylab.plot(eptk_FT*x,eptk_FT*y,label='FT')
    pylab.plot(eptk_L*x,eptk_L*y,label='L')
    pylab.legend()
    pylab.show()
