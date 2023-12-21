#This file will take as an input the temperature and the strength and will run the network model for every set of sample parameter values inputed. It will then save this in a suitably labelled file.

#Import packages
import numpy as np
import pycascades as pc
import sys
from pycascades.earth_system.earth import linear_coupling_earth_system


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
        gis = pc.cusp(a=-1 / self._gis_time, b=1 / self._gis_time, c=(1 / self._gis_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_gis, effective_GMT), x_0=0)
        thc = pc.cusp(a=-1 / self._thc_time, b=1 / self._thc_time, c=(1 / self._thc_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_thc, effective_GMT), x_0=0)
        wais = pc.cusp(a=-1 / self._wais_time, b=1 / self._wais_time, c=(1 / self._wais_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_wais, effective_GMT), x_0=0)
        amaz = pc.cusp(a=-1 / self._amaz_time, b=1 / self._amaz_time, c=(1 / self._amaz_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_amaz, effective_GMT), x_0=0)
        nino = pc.cusp(a=-1 / self._nino_time, b=1 / self._nino_time, c=(1 / self._nino_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_nino, effective_GMT), x_0=0)
        assi = pc.cusp(a=-1 / self._assi_time, b=1 / self._assi_time, c=(1 / self._assi_time) * pc.earth_system.functions_earth_system.global_functions.CUSPc(0., self._limits_assi, effective_GMT), x_0=0)


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
    times=np.array((tau_gis,tau_thc,tau_wais,tau_nino,tau_amaz,tau_assi))
    time_diffs=np.zeros((6,6))
    for i in range(0,6):
        for j in range(0,6):
            if i!=j:
                time_diffs[i,j]=abs(times[i]-times[j])
            else:
                time_diffs[i,j]=np.nan
    min_time_of_interest=np.nanmin((np.nanmin(times),np.nanmin(time_diffs)))
    timestep=np.nanmax((1,int(min_time_of_interest/5)))
    ev.integrate(timestep, t_end)
    timescales=np.zeros((6))
    for i in range(0,6):
        if ev.get_timeseries()[1][-1, i]<0:
            timescales[i]=np.nan
        else:
            timescales[i]=np.nanmin(np.where(ev.get_timeseries()[1][:, i]>0))
    return timescales









#Load the samples
samples=np.load("SALib_samples_6_var.npy")

#Define the number of years to run it for, the timestep, and retrieve the input link strength and global warming values
duration = 100000.
timestep = 1
strength,GMT=np.array(sys.argv[1:],dtype=float)

output=np.zeros((np.shape(samples)[0],np.shape(samples)[1]+8))


#Iterate through and run the network model for each sample
for i in range(0,np.shape(samples)[0]):
    limits=samples[i,:].copy()
    results_limits=limits.copy()
    output[i,:2]=strength,GMT
    output[i,2:-6]=results_limits
    output[i,-6:]=jon_6_var_earth_system_function(duration,timestep,strength,GMT,*list(limits))
    if i%104==0:
        np.save("jon_SALib_outputs_6_var/GMT_"+str(GMT)+"_strength_"+str(strength)+"_SALib_test_results_6_var_"+str(i)+".npy",output)

np.save("jon_SALib_outputs_6_var/GMT_"+str(GMT)+"_strength_"+str(strength)+"_SALib_test_results_6_var_final.npy",output)




