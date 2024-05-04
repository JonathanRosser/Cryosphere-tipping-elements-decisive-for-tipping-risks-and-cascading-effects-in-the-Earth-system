#This file will take as an input the temperature and the strength and will run the network model with these parameters on each sample parameter set in the input file, outputting the results to a datafile labelled with those values

#Import packages
import numpy as np
import pycascades as pc
import sys
import os
import itertools
from pycascades.earth_system.earth import linear_coupling_earth_system

#Define classes for relevant tipping points

class cusp(pc.core.tipping_element.tipping_element):
    """Concrete class for cusp-like tipping element"""
    def __init__(self, a = -4, b = 1, c = 0, x_0 = 0.5 ):
        """Constructor with additional parameters for cusp"""
        super().__init__()
        self._type = 'cusp'
        self._par['a'] = a
        self._par['b'] = b
        self._par['c'] = c
        self._par['x_0'] = x_0

    def dxdt_diag(self):
        """returns callable of dx/dt diagonal element of cusp"""
        return lambda t, x : self._par['a'] * pow(x - self._par['x_0'],3) \
                           + self._par['b'] * (x - self._par['x_0']) \
                           + self._par['c']

    def jac_diag(self):
        """returns callable jacobian diagonal element of cusp."""
        return lambda t, x : 3 * self._par['a'] * pow(x - self._par['x_0'],2) \
                               + self._par['b']

    def tip_state(self):
        return lambda x : x > self._par['x_0']


class linear(pc.core.tipping_element.tipping_element):
    """Concrete class for cusp-like tipping element"""
    def __init__(self, a = -4, c = 0, x_0 = 0.5 ):
        """Constructor with additional parameters for cusp"""
        super().__init__()
        self._type = 'cusp'
        self._par['a'] = a
        self._par['c'] = c
        self._par['x_0'] = x_0

    def dxdt_diag(self):
        """returns callable of dx/dt diagonal element of cusp"""
        return lambda t, x : self._par['a'] * (x - self._par['x_0']) \
                           + self._par['c']

    def jac_diag(self):
        """returns callable jacobian diagonal element of cusp."""
        return lambda t, x : self._par['a']

    def tip_state(self):
        return lambda x : x > self._par['x_0']


class cubic(pc.core.tipping_element.tipping_element):
    """Concrete class for cusp-like tipping element"""
    def __init__(self, a = -4, c = 0, x_0 = 0.5 ):
        """Constructor with additional parameters for cusp"""
        super().__init__()
        self._type = 'cusp'
        self._par['a'] = a
        self._par['c'] = c
        self._par['x_0'] = x_0

    def dxdt_diag(self):
        """returns callable of dx/dt diagonal element of cusp"""
        return lambda t, x : self._par['a'] * pow(x - self._par['x_0'],3) \
                           + self._par['c']

    def jac_diag(self):
        """returns callable jacobian diagonal element of cusp."""
        return lambda t, x : 3 * self._par['a'] * pow(x - self._par['x_0'],2) \

    def tip_state(self):
        return lambda x : x > self._par['x_0']



#Define a class based on the pycascades system to set up the network and run the model

class Earth_System_6var_fixed_links():
    def __init__(self, limits_gis,limits_thc,limits_wais,limits_amaz,limits_nino,limits_assi,
            pf_wais_to_gis,pf_thc_to_gis,pf_gis_to_thc,pf_wais_to_thc,pf_assi_to_thc,pf_thc_to_wais,pf_gis_to_wais,pf_nino_to_wais,pf_thc_to_amaz,pf_nino_to_amaz,pf_thc_to_nino,pf_thc_to_assi,
            tau_gis,tau_thc,tau_wais,tau_nino,tau_amaz,tau_assi):

        #timescales
        self._gis_time = tau_gis
        self._thc_time = tau_thc
        self._wais_time = tau_wais
        self._amaz_time = tau_amaz
        self._nino_time = tau_nino
        self._assi_time = tau_assi

        #tipping limits
        self._limits_gis = limits_gis
        self._limits_thc = limits_thc
        self._limits_wais = limits_wais
        self._limits_amaz = limits_amaz
        self._limits_nino = limits_nino
        self._limits_assi = limits_assi


        #probability fractions
        self._pf_wais_to_gis = pf_wais_to_gis
        self._pf_thc_to_gis = pf_thc_to_gis
        self._pf_gis_to_thc = pf_gis_to_thc
        self._pf_wais_to_thc = pf_wais_to_thc
        self._pf_assi_to_thc = pf_assi_to_thc
        self._pf_gis_to_wais = pf_gis_to_wais
        self._pf_thc_to_wais = pf_thc_to_wais
        self._pf_nino_to_wais = pf_nino_to_wais
        self._pf_thc_to_amaz = pf_thc_to_amaz
        self._pf_nino_to_amaz = pf_nino_to_amaz
        self._pf_thc_to_nino = pf_thc_to_nino
        self._pf_thc_to_assi = pf_thc_to_assi
    def earth_network(self, effective_GMT, strength):
        gis = cusp(a=-1 / self._gis_time, b=1 / self._gis_time, c=(1 / self._gis_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_gis, effective_GMT), x_0=0)
        thc = cusp(a=-1 / self._thc_time, b=1 / self._thc_time, c=(1 / self._thc_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_thc, effective_GMT), x_0=0)
        wais = cusp(a=-1 / self._wais_time, b=1 / self._wais_time, c=(1 / self._wais_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_wais, effective_GMT), x_0=0)
        amaz = cusp(a=-1 / self._amaz_time, b=1 / self._amaz_time, c=(1 / self._amaz_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_amaz, effective_GMT), x_0=0)
        nino = linear(a=-1 / self._nino_time, c=(1 / self._nino_time) * effective_GMT/self._limits_nino, x_0=-1)
        assi = linear(a=-1 / self._assi_time, c=(1 / self._assi_time) * effective_GMT/self._limits_assi, x_0=-1)


        # set up network
        net = pc.earth_system.tipping_network_earth_system.tipping_network()
        net.add_element(gis)
        net.add_element(thc)
        net.add_element(wais)
        net.add_element(amaz)
        net.add_element(nino)
        net.add_element(assi)


        ######################################Set edges to active state#####################################
        net.add_coupling(1, 0, linear_coupling_earth_system(strength=(1 / self._gis_time) * strength * self._pf_thc_to_gis, x_0=-1.0))
        net.add_coupling(2, 0, linear_coupling_earth_system(strength=(1 / self._gis_time) * strength * self._pf_wais_to_gis, x_0=-1.0))

        net.add_coupling(0, 1, linear_coupling_earth_system(strength=(1 / self._thc_time) * strength * self._pf_gis_to_thc, x_0=-1.0))
        net.add_coupling(2, 1, linear_coupling_earth_system(strength=(1 / self._thc_time) * strength * self._pf_wais_to_thc, x_0=-1.0))
        net.add_coupling(5, 1, linear_coupling_earth_system(strength=(1 / self._thc_time) * strength * self._pf_assi_to_thc, x_0=-1.0))

        net.add_coupling(0, 2, linear_coupling_earth_system(strength=(1 / self._wais_time) * strength * self._pf_gis_to_wais, x_0=-1.0))
        net.add_coupling(1, 2, linear_coupling_earth_system(strength=(1 / self._wais_time) * strength * self._pf_thc_to_wais, x_0=-1.0))
        net.add_coupling(4, 2, linear_coupling_earth_system(strength=(1 / self._wais_time) * strength * self._pf_nino_to_wais, x_0=-1.0))


        net.add_coupling(1, 3, linear_coupling_earth_system(strength=(1 / self._amaz_time) * strength * self._pf_thc_to_amaz, x_0=-1.0))
        net.add_coupling(4, 3, linear_coupling_earth_system(strength=(1 / self._amaz_time) * strength * self._pf_nino_to_amaz, x_0=-1.0))

        net.add_coupling(1, 4, linear_coupling_earth_system(strength=(1 / self._nino_time) * strength * self._pf_thc_to_nino, x_0=-1.0))

        net.add_coupling(1, 5, linear_coupling_earth_system(strength=(1 / self._assi_time) * strength * self._pf_thc_to_assi, x_0=-1.0))




        return net





#Define a function to run the model and return the timescales:

def jon_6_var_earth_system_function(duration, timestep, strength, GMT, limits_gis,limits_thc,limits_wais,limits_amaz,limits_nino,limits_assi,pf_wais_to_gis,pf_thc_to_gis,pf_gis_to_thc,pf_wais_to_thc,pf_assi_to_thc,pf_thc_to_wais,pf_gis_to_wais,pf_nino_to_wais,pf_thc_to_amaz,pf_nino_to_amaz,pf_thc_to_nino,pf_thc_to_assi,tau_gis,tau_thc,tau_wais,tau_nino,tau_amaz,tau_assi):
    earth_system = Earth_System_6var_fixed_links(limits_gis,limits_thc,limits_wais,limits_amaz,limits_nino,limits_assi,pf_wais_to_gis,pf_thc_to_gis,pf_gis_to_thc,pf_wais_to_thc,pf_assi_to_thc,pf_thc_to_wais,pf_gis_to_wais,pf_nino_to_wais,pf_thc_to_amaz,pf_nino_to_amaz,pf_thc_to_nino,pf_thc_to_assi,tau_gis,tau_thc,tau_wais,tau_nino,tau_amaz,tau_assi)
    initial_state = [-1,-1,-1,-1,-1,-1]
    net = earth_system.earth_network(GMT, strength)
    ev = pc.evolve(net, initial_state)
    t_end = duration
    ev.integrate(timestep, t_end)
    timescales=np.zeros((6))
    timeseries=ev.get_timeseries()
    for i in range(0,6):
        if ev.get_timeseries()[1][-1, i]<0:
            timescales[i]=np.nan
        else:
            timescales[i]=np.nanmin(np.where(ev.get_timeseries()[1][:, i]>0))
    return timescales



#Load the LHC values
LHS=np.load("jon_LHS_6var.npy")

#Define the number of years to run it for, the timestep, and retrieve the input link strength and global warming values
duration = 100000.
timestep = 1
strength,GMT=np.array(sys.argv[1:],dtype=float)


limits_names=["limits_gis", "limits_thc", "limits_wais", "limits_amaz","limits_nino", "limits_assi", "pf_wais_to_gis", "pf_thc_to_gis", "pf_gis_to_thc", "pf_wais_to_thc", "pf_assi_to_thc", "pf_thc_to_wais", "pf_gis_to_wais", "pf_nino_to_wais", "pf_thc_to_amaz", "pf_nino_to_amaz", "pf_thc_to_nino", "pf_thc_to_assi", "none"]


#Create the output array

output=np.zeros((np.shape(LHS)[0],np.shape(limits_names)[0],np.shape(LHS)[1]+8))
#Run the network model on each sample leaving out/eliminating each relevant element or link in the system once
for i in range(0,np.shape(LHS)[0]):
    for k in range(0,len(limits_names)):
        limits=LHS[i,:].copy()

        #limits_name=limits_names[k]
        #print(k)
        #Remove the relevant element for the leave one out analysis
        if k<6:
            limits[k]=10**9
            if k == 0:
                limits[6:9]=0
                limits[12]=0
            elif k == 1:
                limits[7:12]=0
                limits[13]=0
                limits[15]=0
                limits[17]=0
            elif k == 2:
                limits[11:13]=0
                limits[6]=0
                limits[9]=0
            elif k == 3:
                limits[13]=0
                limits[14]=0
            elif k == 4:
                limits[14:17]=0
            elif k == 5:
                limits[10]=0
                limits[17]=0
        #elif k>17 and k<24:
        #    limits[k]=10**20
        elif k<len(limits_names)-1:
            limits[k]=0
        results_limits=limits.copy()
        output[i,k,:2]=strength,GMT
        output[i,k,2:-6]=results_limits
        output[i,k,-6:]=jon_6_var_earth_system_function(duration,timestep,strength,GMT,*list(limits))
        #Periodically save the data
    if i%100==0:
        np.save("jon_6var_loo_outputs/GMT_"+str(GMT)+"_strength_"+str(strength)+"_results_"+str(i)+".npy",output)

#Save the final datafile
np.save("jon_6var_loo_outputs/GMT_"+str(GMT)+"_strength_"+str(strength)+"_results_final.npy",output)
