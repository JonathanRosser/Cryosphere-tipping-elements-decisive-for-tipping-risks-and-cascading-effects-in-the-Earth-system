import numpy as np
import matplotlib.pyplot as plt

GMT=[1.5,4.0]
strength=0.5


data_location="/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/LOO_Fixed_Signs_Fixed_Links_Updated_Structure_Linear/jon_6var_loo_outputs/"

limits_names=["limits_gis", "limits_thc", "limits_wais", "limits_amaz","limits_nino", "limits_assi", "pf_wais_to_gis", "pf_thc_to_gis", "pf_gis_to_thc", "pf_wais_to_thc", "pf_assi_to_thc", "pf_thc_to_wais", "pf_gis_to_wais",  "pf_thc_to_amaz", "pf_nino_to_amaz", "pf_thc_to_nino", "pf_amaz_to_nino", "pf_thc_to_assi", "none"]
old_label_names=["WAIS","AMAZ","NINO","ASSI","GIS","AMOC"]
label_names=["WAIS","GIS","AMOC","AMAZ","ASSI","ENSO"]

labels=["Control","AMOC"+r'$\rightarrow$'+"GIS Removed","GIS"+r'$\rightarrow$'+"WAIS Removed"]

old_colors=[(202/256,178/256,214/256),(178/256,223/256,138/256),(251/256,154/256,153/256),(253/256,191/256,111/256),(255/256,255/256,153/256),(166/256,206/256,227/256)]

colors=[(202/256,178/256,214/256),(255/256,255/256,153/256),(166/256,206/256,227/256),(178/256,223/256,138/256),(253/256,191/256,111/256),(251/256,154/256,153/256)]



node_values_sorted={}
node_std_sorted={}
node_names_sorted={}
link_values_sorted={}
link_std_sorted={}
link_names_sorted={}
mean_control={}
std_control={}
mean_missing={}
std_missing={}
mean_removed={}
std_removed={}
amoc_removed={}
amocstd_removed={}
gisfrac_removed={}

rearrange_indices=[0,4,5,1,3,2]


for gmt in GMT:
        data=np.load(data_location+"GMT_"+str(gmt)+"_strength_"+str(strength)+"_results_final.npy")[:1000,:,:]
        tipped=np.isnan(data[:,:,-6:])==False
        mean_tipped=np.nanmean(tipped,axis=0)
        num_limits=np.shape(data)[1]
        means=np.zeros((10,num_limits,6))
        mean_missing_calc=np.zeros((10,6))
        interval=int(len(tipped)/10)
        for i in range(0,10):
            means[i,:,:]=np.nanmean(tipped[i*interval:(i+1)*interval],axis=0)
            missing=tipped[i*interval:(i+1)*interval,-1]
            missing[:,0]=0
            mean_missing_calc[i]=np.nanmean(missing,axis=0)
        std_tipped=np.nanstd(means,axis=0)/np.sqrt(10)
        std_missing[gmt]=np.nanstd(mean_missing_calc,axis=0)/np.sqrt(10)

        mean_control[gmt]=mean_tipped[-1,:]
        std_control[gmt]=std_tipped[-1,:]
        mean_removed[gmt]=mean_tipped[12,:]
        amoc_removed[gmt]=mean_tipped[7,:]
        amocstd_removed[gmt]=std_tipped[7,:]
        std_removed=std_tipped[12,:]
        mean_missing[gmt]=mean_tipped[-1,:].copy()
        mean_missing[gmt][0]=0




mean_control[1.5]=np.roll(mean_control[1.5],-2)[rearrange_indices]
mean_removed[1.5]=np.roll(mean_removed[1.5],-2)[rearrange_indices]
amoc_removed[1.5]=np.roll(amoc_removed[1.5],-2)[rearrange_indices]



mean_control[4.0]=np.roll(mean_control[4.0],-2)[rearrange_indices]
mean_removed[4.0]=np.roll(mean_removed[4.0],-2)[rearrange_indices]
amoc_removed[4.0]=np.roll(amoc_removed[4.0],-2)[rearrange_indices]



plt.close()
fig, axs = plt.subplots(3, 1,sharey="row",sharex="none",figsize=(8,10),layout="constrained")

gridspec = axs[0].get_subplotspec().get_gridspec()

axs[2].remove()

#fig.suptitle("Impact of link removal on tipping")
x=np.arange(len(labels))
width=0.12
axs[0].set_title("1.5\N{DEGREE SIGN}C")
for i in range(0,6):
    axs[0].bar(x+(i-3)*width+0.06,(mean_control[1.5][i],amoc_removed[1.5][i],mean_removed[1.5][i]),width,label=label_names[i],color=colors[i])

axs[0].vlines(x[:-1]+0.5,ymin=0,ymax=1,colors="k",linewidth=1)


#axs[0].errorbar(labels,(np.nansum(mean_control[1.5]),mean_var_missing[1.5],mean_var_removed[1.5]),yerr=(std_control[1.5],std_missing[1.5],std_removed[1.5]),fmt="none",color="k",capsize=3.0)
axs[1].set_title("4.0\N{DEGREE SIGN}C")
for i in range(0,6):
    axs[1].bar(x+(i-3)*width+0.06,(mean_control[4.0][i],amoc_removed[4.0][i],mean_removed[4.0][i]),width,label=label_names[i],color=colors[i])

axs[1].vlines(x[:-1]+0.5,ymin=0,ymax=1,colors="k",linewidth=1)


#axs[1].errorbar(labels,(np.nansum(mean_control[4.0]),mean_var_missing[4.0],mean_var_removed[4.0]),yerr=(std_control[4.0],std_missing[4.0],std_removed[4.0]),fmt="none",color="k",capsize=3.0)
axs[0].legend(reverse=False)
axs[0].set_ylim([0,1])
axs[1].set_ylim([0,1])
axs[0].set_ylabel("Fraction Tipped")
axs[1].set_ylabel("Fraction Tipped")

axs[0].set_xticks([])
axs[1].set_xticks(x,labels,rotation=15)


#axs[0].set_xticklabels(labels,rotation=45)
#axs[1].set_xticklabels(labels,rotation=45)

axs[0].text(0.02,0.95,"a",fontweight="bold",transform=axs[0].transAxes)
axs[1].text(0.02,0.95,"b",fontweight="bold",transform=axs[1].transAxes)


subfig=fig.add_subfigure(gridspec[2,0])
subfig_axs=subfig.subplots(1,2,sharey="row")


subfig_axs[0].bar(labels,(0,-100*(1-np.nansum(amoc_removed[1.5])/np.nansum(mean_control[1.5])),-100*(1-np.nansum(mean_removed[1.5])/np.nansum(mean_control[1.5]))),color="c",fill=False,hatch="//",edgecolor="c")
subfig_axs[0].axhline(y=0,color="k",linewidth=0.8)
subfig_axs[1].bar(labels,(0,-100*(1-np.nansum(amoc_removed[4.0])/np.nansum(mean_control[4.0])),-100*(1-np.nansum(mean_removed[4.0])/np.nansum(mean_control[4.0]))),color="c",fill=False,hatch="//",edgecolor="c")
subfig_axs[0].set_ylabel("Percentage Change")

subfig_axs[0].set_xticklabels(labels,rotation=45)
subfig_axs[1].set_xticklabels(labels,rotation=45)
subfig_axs[1].axhline(y=0,color="k",linewidth=0.8)
subfig_axs[0].set_title("1.5\N{DEGREE SIGN}C")
subfig_axs[1].set_title("4.0\N{DEGREE SIGN}C")
subfig_axs[0].text(0.05,0.9,"c",fontweight="bold",transform=subfig_axs[0].transAxes)
subfig_axs[1].text(0.05,0.9,"d",fontweight="bold",transform=subfig_axs[1].transAxes)


#plt.tight_layout()
#plt.subplots_adjust(wspace=0,hspace=0.1,top=0.9)
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures_Updated_Structure_Linear/Figures/LOO_Links_Removal_Expanded_Bottom_Panel.png")
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures_Updated_Structure_Linear/Figures/LOO_Links_Removal_Expanded_Bottom_Panel.eps")



