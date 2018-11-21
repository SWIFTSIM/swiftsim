#!/usr/bin/env python
import numpy as np
from hmf import MassFunction
import hmf
from astropy.cosmology import FlatLambdaCDM


def getHMF(H0=70.3, Om0=0.276, Ob0=0.0455, Tcmb0=2.725, Mmin=1e10,Mmax=1e15):
    ''' Fast function to call the HMF from hmf, this function only has 
        6 variables and will return the dn/d(log10 M) and M array.
        H0: Hubble constant
        Om0: Matter density
        Ob0: Baryon density
        Tcmb0: CMB temperature at z=0
        Mmin: minimum mass (solar masses)
        Mmax: Maximum mass (solar masses) 
    '''
    new_model = FlatLambdaCDM(H0 = H0, Om0=Om0, Ob0 = Ob0, Tcmb0=Tcmb0)
    hmff = MassFunction(cosmo_model=new_model, Mmax = np.log10(Mmax), Mmin=np.log10(Mmin),hmf_model='ST')
    return hmff.m, hmff.dndlog10m

def getHMFz(z,H0=70.3, Om0=0.276, Ob0=0.0455, Tcmb0=2.725, Mmin=1e10,Mmax=1e15):
    ''' Fast function to call the HMF from hmf, this function only has 
        7 variables and will return the dn/d(log10 M) and M array.
        z: redshift
        H0: Hubble constant
        Om0: Matter density
        Ob0: Baryon density
        Tcmb0: CMB temperature at z=0
        Mmin: minimum mass (solar masses)
        Mmax: Maximum mass (solar masses) 
    '''
    new_model = FlatLambdaCDM(H0 = H0, Om0=Om0, Ob0 = Ob0, Tcmb0=Tcmb0)
    hmff = MassFunction(cosmo_model=new_model, Mmax = np.log10(Mmax), Mmin=np.log10(Mmin),z=z,hmf_model='ST')
    return hmff.m, hmff.dndlog10m

def getHMFztinker(z,H0=70.3, Om0=0.276, Ob0=0.0455, Tcmb0=2.725, Mmin=1e10,Mmax=1e15):
    ''' Fast function to call the HMF from hmf, this function only has 
        6 variables and will return the dn/d(log10 M) and M array.
        H0: Hubble constant
        Om0: Matter density
        Ob0: Baryon density
        Tcmb0: CMB temperature at z=0
        Mmin: minimum mass (solar masses)
        Mmax: Maximum mass (solar masses) 
    '''
    new_model = FlatLambdaCDM(H0 = H0, Om0=Om0, Ob0 = Ob0, Tcmb0=Tcmb0)
    hmff = MassFunction(cosmo_model=new_model, Mmax = np.log10(Mmax), Mmin=np.log10(Mmin),z=z)
    return hmff.m, hmff.dndlog10m

'''
hmmf = MassFunction()
hmmf.parameter_info()
print('Statement 2:')
print(hmmf.hmf_model)
model_st = hmf.fitting_functions.ST([1,2])
print(model_st)

print(getHMFztinker(10))


#print(getHMFztinker(1.0))

'''
