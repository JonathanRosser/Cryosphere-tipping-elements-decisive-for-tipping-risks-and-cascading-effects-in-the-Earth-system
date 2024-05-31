Code used to generate data for tipping point uncertainty analysis for the "Cryosphere tipping elements decisive for tipping risks and cascading effects in the Earth system" paper

The packages required to run this analysis are:
numpy
pycascades
pyDOE
SALib


Files and outputs:
Sobol Analysis, Figures 1 & 2: 
Jon_6_Var_SAlib_Sample.py - generates a set of sample matrices of input parameters such as tipping thresholds and interaction strengths based on distributions estimated from the literature. 
Jon_6_Var_SAlib_Run.py - Runs the PyCascades model over the set of sample matrices of input parameters and then performs Sobol analysis on these outputs. This generates the total effect indices for the inputs which are used to generate Figures 1 & 2.

Leave One Out Analysis, Figures 3,4,5 & 6:
jon_LHS_6Var.py - Generates samples of input parameters evenly distributed across distributions estimated from the literature using the latin hypercube sampling technique.
Jon_6_Var_LOO_run.py - Runs the Pycascades model over the input samples
LOO_GISAMOC_Missing_Tipping_Risk.py - Calculates the risk of elements tipping based on the outputs of the Leave One Out model runs
LOO_Node_Missing_percentage_change_boxplot.py - Uses leave one out model runs to calculate the percentage change in the mean number of elements tipped when different elements are removed and create Figure 3.
LOO_GISWAISAMOC_Missing_Comp_expanded_Bottom_Panel.py - Uses leave one out model runs to calculate the fraction of runs in which tipping elements tipped in different scenarios when the GIS, WAIS, and AMOC are removed and a control scenario where none were removed. This is plotted in Figure 4.
LOO_Link_Missing_Comp_percentage_change_boxplot.py - Uses leave one out model runs to calculate the percentage change in the mean number of elements tipped when different links between elements are removed and create Figure 5.
LOO_Missing_Links_Compressed_Bottom_Panel.py - Uses leave one out model runs to calculate the fraction of runs in which tipping elements tipped in different scenarios when the AMOC -> GIS link and the GIS -> WAIS link were removed compared to a control run. This is plotted in Figure 6.
