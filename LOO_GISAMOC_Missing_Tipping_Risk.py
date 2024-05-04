import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['hatch.linewidth'] = 1.0



GMT=[1.5,4.0]
strength=0.5


data_location="/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/LOO_Fixed_Signs_Fixed_Links_Updated_Structure_Linear/jon_6var_loo_outputs/"

limits_names=["limits_gis", "limits_thc", "limits_wais", "limits_amaz","limits_nino", "limits_assi", "pf_wais_to_gis", "pf_thc_to_gis", "pf_gis_to_thc", "pf_wais_to_thc", "pf_assi_to_thc", "pf_thc_to_wais", "pf_gis_to_wais", "pf_nino_to_wais",  "pf_nino_to_amaz", "pf_thc_to_amaz", "pf_thc_to_nino", "pf_thc_to_assi", "none"]
label_names=["WAIS","AMAZ","NINO","ASSI","GIS","AMOC"]

labels=["Control","GIS Removed","AMOC Removed"]

#colors=["mediumpurple","palegreen","coral","turquoise","silver","cornflowerblue"]

colors=[(202/256,178/256,214/256),(178/256,223/256,138/256),(251/256,154/256,153/256),(253/256,191/256,111/256),(255/256,255/256,153/256),(166/256,206/256,227/256)]


tip_test={}

for gmt in GMT:
        data=np.load(data_location+"GMT_"+str(gmt)+"_strength_"+str(strength)+"_results_final.npy")[:1000,:,:]
        tipped=np.isnan(data[:,:,-6:])==False
        sum_tipped=np.nansum(tipped,axis=2)
        tip_test[gmt]=np.nanmean(sum_tipped>0,axis=0)
        


for i in range(0,len(limits_names)):
    print(limits_names[i],tip_test[1.5][i],-(tip_test[1.5][-1]-tip_test[1.5][i])/tip_test[1.5][-1])

