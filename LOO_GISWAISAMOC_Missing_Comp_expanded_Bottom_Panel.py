import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['hatch.linewidth'] = 1.0



GMT=[1.5,4.0]
strength=0.5


data_location="/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/LOO_Fixed_Signs_Fixed_Links_Updated_Structure_Linear/jon_6var_loo_outputs/"

limits_names=["limits_gis", "limits_thc", "limits_wais", "limits_amaz","limits_nino", "limits_assi", "pf_wais_to_gis", "pf_thc_to_gis", "pf_gis_to_thc", "pf_wais_to_thc", "pf_assi_to_thc", "pf_thc_to_wais", "pf_gis_to_wais",  "pf_thc_to_amaz", "pf_nino_to_amaz", "pf_thc_to_nino", "pf_amaz_to_nino", "pf_thc_to_assi", "none"]
label_names=["WAIS","AMAZ","NINO","ASSI","GIS","AMOC"]

labels=["Control","GIS Removed","WAIS Removed", "AMOC Removed"]

#colors=["mediumpurple","palegreen","coral","turquoise","silver","cornflowerblue"]

colors=[(202/256,178/256,214/256),(178/256,223/256,138/256),(251/256,154/256,153/256),(253/256,191/256,111/256),(255/256,255/256,153/256),(166/256,206/256,227/256)]

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
gis_removed={}
gisstd_removed={}
wais_removed={}
waisstd_removed={}

amoc_removed={}
amocstd_removed={}
gisfrac_removed={}


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
        gis_removed[gmt]=mean_tipped[0,:]
        amoc_removed[gmt]=mean_tipped[1,:]
        amocstd_removed[gmt]=std_tipped[1,:]
        gisstd_removed=std_tipped[0,:]
        wais_removed[gmt]=mean_tipped[2,:]
        waisstd_removed[gmt]=std_tipped[2,:]

        mean_missing[gmt]=mean_tipped[-1,:].copy()
        mean_missing[gmt][0]=0




mean_control[1.5]=np.roll(mean_control[1.5],-2)
gis_removed[1.5]=np.roll(gis_removed[1.5],-2)
amoc_removed[1.5]=np.roll(amoc_removed[1.5],-2)
wais_removed[1.5]=np.roll(wais_removed[1.5],-2)



mean_control[4.0]=np.roll(mean_control[4.0],-2)
gis_removed[4.0]=np.roll(gis_removed[4.0],-2)
amoc_removed[4.0]=np.roll(amoc_removed[4.0],-2)
wais_removed[4.0]=np.roll(wais_removed[4.0],-2)



plt.close()

fig, axs = plt.subplots(3, 1,sharey="row",sharex="none",figsize=(8,10),layout="constrained")
gridspec = axs[0].get_subplotspec().get_gridspec()


#fig.suptitle("Impact of node removal on tipping")
x=np.arange(len(labels))
axs[0].set_title("1.5\N{DEGREE SIGN}C")
width=0.12
for i in range(0,6):
    axs[0].bar(x+(i-3)*width+0.06,(mean_control[1.5][i],gis_removed[1.5][i],wais_removed[1.5][i],amoc_removed[1.5][i]),width,label=label_names[i],color=colors[i])

axs[0].vlines(x[:-1]+0.5,ymin=0,ymax=1,colors="k",linewidth=1)

#axs[0].errorbar(labels,(np.nansum(mean_control[1.5]),mean_var_missing[1.5],mean_var_removed[1.5]),yerr=(std_control[1.5],std_missing[1.5],std_removed[1.5]),fmt="none",color="k",capsize=3.0)
axs[0].set_ylabel("Fraction Tipped")
axs[0].set_ylim([0,1])
axs[1].set_title("4.0\N{DEGREE SIGN}C")

for i in range(0,6):
    axs[1].bar(x+(i-3)*width+0.06,(mean_control[4.0][i],gis_removed[4.0][i],wais_removed[4.0][i],amoc_removed[4.0][i]),width,label=label_names[i],color=colors[i])

axs[1].vlines(x[:-1]+0.5,ymin=0,ymax=1,colors="k",linewidth=1)


#axs[1].errorbar(labels,(np.nansum(mean_control[4.0]),mean_var_missing[4.0],mean_var_removed[4.0]),yerr=(std_control[4.0],std_missing[4.0],std_removed[4.0]),fmt="none",color="k",capsize=3.0)
#axs[0,0].legend(reverse=True)
axs[0].legend(loc=1)
axs[1].set_ylabel("Fraction Tipped")
axs[1].set_ylim([0,1])
axs[0].set_xticks([])
axs[1].set_xticks(x,labels,rotation=15)

axs[0].text(0.02,0.95,"A",fontweight="bold",transform=axs[0].transAxes)
axs[1].text(0.02,0.95,"B",fontweight="bold",transform=axs[1].transAxes)
axs[2].remove()

location=gridspec[2,0]

subfig=fig.add_subfigure(gridspec[2,0])
subfig_axs=subfig.subplots(1,2,sharey="row")
subfig_axs[0].bar(labels,(0,-100*(1-np.nansum(gis_removed[1.5])/np.nansum(mean_control[1.5])),-100*(1-np.nansum(wais_removed[1.5])/np.nansum(mean_control[1.5])),-100*(1-np.nansum(amoc_removed[1.5])/np.nansum(mean_control[1.5]))),color="c",fill=False,hatch="//",edgecolor="c")
subfig_axs[1].bar(labels,(0,-100*(1-np.nansum(gis_removed[4.0])/np.nansum(mean_control[4.0])),-100*(1-np.nansum(wais_removed[4.0])/np.nansum(mean_control[4.0])),-100*(1-np.nansum(amoc_removed[4.0])/np.nansum(mean_control[4.0]))),color="c",fill=False,hatch="//",edgecolor="c")
subfig_axs[0].set_ylabel("Percentage Change")
subfig_axs[0].set_xticklabels(labels,rotation=45)
subfig_axs[1].set_xticklabels(labels,rotation=45)
subfig_axs[0].set_title("1.5\N{DEGREE SIGN}C")
subfig_axs[1].set_title("4.0\N{DEGREE SIGN}C")
subfig_axs[0].text(0.05,0.9,"C",fontweight="bold",transform=subfig_axs[0].transAxes)
subfig_axs[1].text(0.05,0.9,"D",fontweight="bold",transform=subfig_axs[1].transAxes)




#subfig.subplots_adjust(wspace=0,hspace=0,left=0)



#plt.tight_layout()
#plt.subplots_adjust(wspace=0,hspace=0.3,top=0.95,left=0.4,right=0.05)
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures_Updated_Structure_Linear/Figures/LOO_GISWAISAMOC_Removal_Expanded.png")
plt.savefig("/home/jonathan/Documents/Coding/ERA/pycascades/earth_system/Paper_Figures_Updated_Structure_Linear/Figures/LOO_GISWAISAMOC_Removal_Expanded.eps")



