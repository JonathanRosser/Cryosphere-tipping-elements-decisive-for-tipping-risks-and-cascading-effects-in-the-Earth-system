import numpy as np
import matplotlib.pyplot as plt

GMT=[1.5,4.0]
strength=0.5


data_location="/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/LOO_Fixed_Signs_Fixed_Links_Updated_Structure_Linear/jon_6var_loo_outputs/"

limits_names=["limits_gis", "limits_thc", "limits_wais", "limits_amaz","limits_nino", "limits_assi", "pf_wais_to_gis", "pf_thc_to_gis", "pf_gis_to_thc", "pf_wais_to_thc", "pf_assi_to_thc", "pf_thc_to_wais", "pf_gis_to_wais",  "pf_thc_to_amaz", "pf_nino_to_amaz", "pf_thc_to_nino", "pf_amaz_to_nino", "pf_thc_to_assi", "none"]
label_names=["GIS","AMOC","WAIS","AMAZ","NINO","ASSI"]

labels=["1.5","4.0"]

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
percent_change={}

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
        var_missing=np.zeros((6))
        for i in range(0,6):
            var_missing[i]=np.nansum(mean_tipped[-1,:])-mean_tipped[-1,i]
        mean_var_missing[gmt]=np.nansum(var_missing)/6
        std_var_missing[gmt]=np.nanstd(var_missing)
        upper_var_missing[gmt]=np.nanmax(var_missing)-mean_var_missing[gmt]
        lower_var_missing[gmt]=mean_var_missing[gmt]-np.nanmin(var_missing)
        mean_var_removed[gmt]=np.nanmean(np.nansum(mean_tipped[6:-1,:],axis=1))
        std_var_removed[gmt]=np.nanstd(np.nansum(mean_tipped[6:-1,:],axis=1))
        upper_var_removed[gmt]=np.nanmax(np.nansum(mean_tipped[6:-1,:],axis=1))-mean_var_removed[gmt]
        lower_var_removed[gmt]=mean_var_removed[gmt]-np.nanmin(np.nansum(mean_tipped[6:-1,:],axis=1))
        var_removed[gmt]=np.nansum(mean_tipped[6:-1,:],axis=1)
        percent_change[gmt]=(var_removed[gmt]-np.nansum(mean_tipped[-1,:]))/np.nansum(mean_tipped[-1,:])*100



max1p5=1-(upper_var_removed[1.5]+mean_var_removed[1.5])/np.nansum(mean_control[1.5])
min1p5=1-(mean_var_removed[1.5]-lower_var_removed[1.5])/np.nansum(mean_control[1.5])
max4p0=1-(upper_var_removed[4.0]+mean_var_removed[4.0])/np.nansum(mean_control[4.0])
min4p0=1-(mean_var_removed[4.0]-lower_var_removed[4.0])/np.nansum(mean_control[4.0])


plt.close()

fig, ax = plt.subplots(figsize=(5,3))
ax.boxplot((percent_change[1.5],percent_change[4.0]),labels=["1.5\N{DEGREE SIGN}C","4.0\N{DEGREE SIGN}C"],whis=100000)
ax.set_ylabel("Percentage change in \n mean number of elements tipped")
ax.set_xlabel("Global Mean Temperature")
#ax.set_title("Impact of link removal at different temperatures")
ax.grid(axis="y",which="major",color="lightgrey")
plt.tight_layout()
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures_Updated_Structure_Linear/Figures/LOO_Link_Removal_Boxplot.png")
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures_Updated_Structure_Linear/Figures/LOO_Link_Removal_Boxplot.eps")


"""
fig, axs = plt.subplots(1, 2,sharey="row")

fig.suptitle("Impact of link removal on tipping")

axs[0].set_title("1.5")
axs[0].bar(labels,(np.nansum(mean_control[1.5]),mean_var_removed[1.5]))
axs[0].errorbar(labels,(np.nansum(mean_control[1.5]),mean_var_removed[1.5]),yerr=((0,lower_var_removed[1.5]),(0,upper_var_removed[1.5])),fmt="none",color="k",capsize=3.0)
axs[0].set_ylabel("Mean Number of Elements Tipped")
axs[0].set_xticklabels(labels,rotation=45)
axs[1].set_title("4.0")
axs[1].bar(labels,(np.nansum(mean_control[4.0]),mean_var_removed[4.0]))
axs[1].errorbar(labels,(np.nansum(mean_control[4.0]),mean_var_removed[4.0]),yerr=((0,lower_var_removed[4.0]),(0,upper_var_removed[4.0])),fmt="none",color="k",capsize=3.0)

axs[1].set_xticklabels(labels,rotation=45)
plt.tight_layout()
plt.subplots_adjust(wspace=0,hspace=0)
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures/Figures/LOO_Link_Removal.png")
"""


