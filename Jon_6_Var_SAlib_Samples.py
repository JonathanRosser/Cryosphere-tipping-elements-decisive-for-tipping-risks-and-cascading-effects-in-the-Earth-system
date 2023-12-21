#This file will take in the desired parameter ranges for the critical temperatures, link strengths, and tipping timescales and generate sample parameter values to perform a Sobol Variance Analysis on these variables

import numpy as np
from SALib import ProblemSpec


#Define the parameter ranges
sp = ProblemSpec({
    "names": ["GIS_Tc", "THC_Tc", "WAIS_Tc", "AMAZ_Tc", "NINO_Tc", "ASSI_Tc", 
        "wais_to_gis", "thc_to_gis", "gis_to_thc", "wais_to_thc", "assi_to_thc", "thc_to_wais", "gis_to_wais", "nino_to_wais", "nino_to_amaz", "thc_to_amaz", "thc_to_nino", "thc_to_assi", 
        "tau_gis", "tau_thc", "tau_wais", "tau_nino", "tau_amaz", "tau_assi"],
    "groups": None,
    "bounds": [[0.8, 3.0], [1.4,8.0], [1.0,3.0], [2.0,6.0], [3.0,6.0], [1.3,2.9], 
        [0.1, 0.2], [-1.0,-0.1], [0.1, 1.], [-0.3, 0.3], [0.1, 0.5],  [0.1,0.15], [0.1, 0.5], [0.1,0.5], [0.1, 1.], [-0.4, 0.4], [0.1, 0.2], [-0.5,-0.1], 
        [1000, 15000], [15, 300], [500, 13000], [25, 200], [50, 200], [10, 50]],
    "outputs": ["num_tipped", "GIS_tipped", "THC_tipped", "WAIS_tipped", "AMAZ_tipped", "NINO_tipped", "ASSI_tipped"],
    })

#Generate the sample parameter values
P=sp.sample_sobol(512,calc_second_order=False)

Samples=sp.samples

#Save the samples to be run through the network model
np.save("SALib_samples_6_var.npy",Samples)
