import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl

mpl.rcParams['hatch.linewidth'] = 1.0


GMT=[1.5,4.0]
strength=0.5


data_location="/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/LOO_Fixed_Signs_Fixed_Links_Updated_Structure_Linear/jon_6var_loo_outputs/"

limits_names=["limits_gis", "limits_thc", "limits_wais", "limits_amaz","limits_nino", "limits_assi", "pf_wais_to_gis", "pf_thc_to_gis", "pf_gis_to_thc", "pf_wais_to_thc", "pf_assi_to_thc", "pf_thc_to_wais", "pf_gis_to_wais",  "pf_thc_to_amaz", "pf_nino_to_amaz", "pf_thc_to_nino", "pf_amaz_to_nino", "pf_thc_to_assi", "none"]
label_names=["GIS","AMOC","WAIS","AMAZ","NINO","ASSI"]

labels=["Node Discounted","Node Removed"]

node_values_sorted={}
node_std_sorted={}
node_names_sorted={}
link_values_sorted={}
link_std_sorted={}
link_names_sorted={}
mean_control={}
std_control={}
mean_var_missing={}
std_var_missing={}
mean_var_removed={}
std_var_removed={}
upper_var_missing={}
lower_var_missing={}
upper_var_removed={}
lower_var_removed={}
var_removed={}
percent_change_var_removed={}
var_missing={}
percent_change_var_missing={}
for gmt in GMT:
        data=np.load(data_location+"GMT_"+str(gmt)+"_strength_"+str(strength)+"_results_final.npy")[:1000,:,:]
        tipped=np.isnan(data[:,:,-6:])==False
        mean_tipped=np.nanmean(tipped,axis=0)
        num_limits=np.shape(data)[1]
        means=np.zeros((10,num_limits,6))
        interval=int(len(tipped)/10)
        for i in range(0,10):
            means[i,:,:]=np.nanmean(tipped[i*interval:(i+1)*interval],axis=0)
        std_tipped=np.nanstd(means,axis=0)/np.sqrt(10)


        mean_control[gmt]=mean_tipped[-1,:]
        std_control[gmt]=std_tipped[-1,:]


        mean_nodes=mean_tipped[:6,:]
        std_nodes=std_tipped[:6,:]
        node_names=limits_names[:6]
        var_missing[gmt]=np.zeros((6))
        for i in range(0,6):
            var_missing[gmt][i]=np.nansum(mean_tipped[-1,:])-mean_tipped[-1,i]
        mean_var_missing[gmt]=np.nansum(var_missing[gmt])/6
        std_var_missing[gmt]=np.nanstd(var_missing[gmt])
        upper_var_missing[gmt]=np.nanmax(var_missing[gmt])-mean_var_missing[gmt]
        lower_var_missing[gmt]=mean_var_missing[gmt]-np.nanmin(var_missing[gmt])
        mean_var_removed[gmt]=np.nanmean(np.nansum(mean_tipped[:6,:],axis=1))
        std_var_removed[gmt]=np.nanstd(np.nansum(mean_tipped[:6,:],axis=1))
        upper_var_removed[gmt]=np.nanmax(np.nansum(mean_tipped[:6,:],axis=1))-mean_var_removed[gmt]
        lower_var_removed[gmt]=mean_var_removed[gmt]-np.nanmin(np.nansum(mean_tipped[:6,:],axis=1))
        var_removed[gmt]=np.nansum(mean_tipped[:6,:],axis=1)
        percent_change_var_removed[gmt]=(var_removed[gmt]-np.nansum(mean_tipped[-1,:]))/np.nansum(mean_tipped[-1,:])*100
        percent_change_var_missing[gmt]=(var_missing[gmt]-np.nansum(mean_tipped[-1,:]))/np.nansum(mean_tipped[-1,:])*100




plt.close()
fig, axs = plt.subplots(1, 2,sharey="row",figsize=(8,4))

#fig.suptitle("Impact of node removal on tipping")

axs[0].set_title("1.5\N{DEGREE SIGN}C")
axs[0].boxplot((percent_change_var_missing[1.5],percent_change_var_removed[1.5]),labels=["Node discounted","Node Removed"],whis=100000)



axs[0].set_ylabel("Percentage change in \n mean number of elements tipped")
axs[0].text(0.05,0.9,"a",fontweight="bold",transform=axs[0].transAxes)
axs[0].set_xticklabels(labels,rotation=45)
axs[0].grid(axis="y",which="major",color="lightgrey")
axs[1].set_title("4.0\N{DEGREE SIGN}C")
axs[1].boxplot((percent_change_var_missing[4.0],percent_change_var_removed[4.0]),labels=["Node discounted","Node Removed"],whis=100000)
axs[1].grid(axis="y",which="major",color="lightgrey")
axs[1].text(0.05,0.9,"b",fontweight="bold",transform=axs[1].transAxes)



axs[1].set_xticklabels(labels,rotation=45)
plt.tight_layout()
plt.subplots_adjust(top=0.87)
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures_Updated_Structure_Linear/Figures/LOO_Node_Removal_Percentage_Change_Boxplot_Linear.png")
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures_Updated_Structure_Linear/Figures/LOO_Node_Removal_Percentage_Change_Boxplot_Linear.eps")



