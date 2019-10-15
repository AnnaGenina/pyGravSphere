import numpy as np

def multiplumsurf(r,pars):

    Mpars = pars[0:len(pars)/2]  #first half are Mi
    apars = pars[len(pars)/2:len(pars)] #second half is ai
    nplum = len(Mpars) #number of plummer components
    multplum = np.zeros(len(r))
    for i in range(len(Mpars)):
        aparsu = apars[i]
        multplum = multplum + Mpars[i]*aparsu**2.0 / (np.pi*(aparsu**2.0+r**2.0)**2.0)
    return multplum


def threeplumsurf(r,M1,M2,M3,a1,a2,a3):
    return multiplumsurf(r,np.array([M1,M2,M3,a1,a2,a3]))

def multiplummass(r,pars):
    
    Mpars = pars[0:len(pars)/2]
    apars = pars[len(pars)/2:len(pars)]
    nplum = len(Mpars)
    multplum = np.zeros(len(r))
    for i in range(len(Mpars)):
        aparsu = apars[i]
        multplum = multplum + Mpars[i]*r**3./(r**2.+aparsu**2.)**(3./2.)
    return multplum

def threeplummass(r,M1,M2,M3,a1,a2,a3):
    return multiplummass(r,np.array([M1,M2,M3,a1,a2,a3]))
