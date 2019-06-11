# coding: utf-8

# In[5]:

#Python script to test the amount of momentum injected by stellar particles to gas particles over the star lifetime. 


from pylab import *
import h5py


# In[157]:

# Fitting formula from Agertz et al. (2013): p = p1 * (t*/t_w) * mstr * (Z/p2)**p3
# with mstr the initial mass of the star and t* the star age.

tw = 6.5       #timescale over which the star can inject momentum (Myr)
p1 = 1.8e40    #momentum per unit stellar mass that can be injected (g cm s^-1 Mo^-1)
p2 = 0.01      #Characteristic metallicity at which p1 is injected (metal mass fraciton)
p3 = 0.38      #Exponent of the metallicity-dependent momentum (dimensionless)
               

def get_momentum(t,mstr,z=0.01):
    """get momentum according to the model used in the code
        and the parameters sourced in the parameter file
    """
    
    if(z < 0.001): z = 0.001
    if(z > 0.04): z = 0.04
    
    return mstr * p1 * (z/p2)**p3 * t/tw


# In[164]:


time_factor = 3.086e21/1e5 / (60*60*24*365*1e6)

momentum_imparted_time = []
momentum_received_time = []
times = []
for i in np.arange(1,201):
    with h5py.File('output_%04d.hdf5'%i,'r') as snap:
        tcurrent = snap['Header'].attrs['Time'] * time_factor
        mgas = snap['PartType0/Masses'].value*1e10
        zstr = snap['PartType4/Metallicity'].value
        mstr = snap['PartType4/InitialMasses'].value*1e10
        tform = snap['PartType4/BirthTime'].value * time_factor
        momentum_gas = snap['PartType0/Momentum Received'].value
        momentum_str = snap['PartType4/Momentum Received'].value
                
        if(tcurrent < tw):
            total_momentum_imparted = get_momentum(tcurrent[0],mstr[0],z=zstr[0])
        else:
            total_momentum_imparted = get_momentum(tw,mstr[0],z=zstr[0])

        total_momentum_received = (np.sum(momentum_gas)) * (1.99e43) * 1e5

        momentum_imparted_time.append(total_momentum_imparted)
        momentum_received_time.append(total_momentum_received)

        times.append(tcurrent)
        
momentum_imparted_time = np.asarray(momentum_imparted_time)
momentum_received_time = np.asarray(momentum_received_time)


# In[165]:


fig = plt.figure(1, figsize=(5,5))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)


ax1.plot(times, np.log10(momentum_imparted_time), lw=2, label='Imparted')
ax1.plot(times, np.log10(momentum_received_time), lw=2, ls='--', label='Received')

ax2.plot(times, np.log10(momentum_received_time)-np.log10(momentum_imparted_time))
ax2.set_ylim(-0.1,0.1)

ax2.set_xlabel(r'$\rm Time / Myr$', size=15)
ax2.set_ylabel(r'$\rm log_{10} \left ( p_{rec} / p_{snd} \right ) $', size=15)
ax1.set_ylabel(r'$\rm log_{10} \left ( p / \ g \ cm \ s^{-1} \right ) $', size=15)

ax1.axvline(tw, ls='--')

ax1.minorticks_on()
ax2.minorticks_on()

ax1.legend(loc='lower right')

plt.savefig('check_momentum.png', bbox_inches='tight', dpi=200)

