#Use a Latin HyperCube Sampler to create sample parameter values to perform Leave One Out Analysis

import numpy as np
from pyDOE import lhs #function name >>> lhs


#Tipping limits, see Schellnhuber, et al., 2016:
limits_gis = [0.8, 3.0]
limits_thc = [1.4,8.0]
limits_wais = [1.0,3.0]
limits_amaz = [2.0,6.0]
limits_nino = [3.0,6.0]
limits_assi = [1.3,2.9]

###################################################
#TO GIS
pf_wais_to_gis = [0.1, 0.2] #Y (Destabilising, Weak)
pf_thc_to_gis = [-1.0,-0.1] #Y (Stabilising, Strong)
# TO THC
pf_gis_to_thc = [0.1, 1.] #Y (Destabilising, Strong)
pf_wais_to_thc = [-0.3, 0.3] #M (unclear, likely timescale dependent, Weak/Moderate)
pf_assi_to_thc=[0.1, 0.5] # Half of GIS Interactions Y (Destabilising, Weak/Moderate)
#pf_nino_to_thc=  Not included, (Unclear, unclear)
# TO WAIS
pf_thc_to_wais = [0.1,0.15] #N (Destabilising, unclear)
pf_gis_to_wais = [0.1, 0.5] #M (Destabilising, Moderate)
pf_nino_to_wais = [0.1,0.5] # N not included (Destabilising, unclear, weak/moderate)
#TO AMAZ
pf_nino_to_amaz = [0.1, 1.]     # from PF = [1., 10.] Y (Destabilising, strong)
# unclear link
pf_thc_to_amaz = [-0.4, 0.4] # Y (Unclear/ Moderate)
# TO NINO
pf_thc_to_nino = [0.1, 0.2]     # from PF = [1., 2.] Y (Destabilising, Weak)
# unclear link
#pf_amaz_to_nino = [-0.15, 0.15]   # from PF = [0.8, 1.5] N (not present in recent analysis?)
# TO ASSI
pf_thc_to_assi = [-0.5,-0.1]     # Half of GIS interactions, Y (Stabilising, Weak/Moderate)

#Timescales
tau_gis  = [1000, 15000]         #1000-15000(central: 10.000)      old values: [1000, 15000] 
tau_thc  = [15, 300]             #15-120 (central: 50)             old values: [15, 300]     
tau_wais = [500, 13000]          #500-13000 (central: 2000)        old values: [1000, 13000] 
tau_nino = [25, 200]             #unclear (around 100)             old values: [25, 200]     
tau_amaz = [50, 200]             #50-200 (central: 100)            old values: [50, 200]     
tau_assi = [10, 50]







limits=np.array((limits_gis,limits_thc,limits_wais,limits_amaz,limits_nino,limits_assi,pf_wais_to_gis,pf_thc_to_gis,pf_gis_to_thc,pf_wais_to_thc,pf_assi_to_thc,pf_thc_to_wais,pf_gis_to_wais,pf_nino_to_wais,pf_thc_to_amaz,pf_nino_to_amaz,pf_thc_to_nino,pf_thc_to_assi,tau_gis,tau_thc,tau_wais,tau_nino,tau_amaz,tau_assi))

#Use Latin HyperCube Sampler to sample parameter values
points = np.array(lhs(24, samples=1000))

array_limits=points*(limits[:,1]-limits[:,0])+limits[:,0]


#Save the parameter values to a file for use in the network model
np.save("jon_LHS_6var.npy",array_limits)





