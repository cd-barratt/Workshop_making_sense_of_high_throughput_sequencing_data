import matplotlib
matplotlib.use('Agg')

import sys
import os
import numpy as np
import dadi
from datetime import datetime
import Models_3D
import pylab
import matplotlib.pyplot as plt

'''
usage: python dadi_3D_04_plotting_functions.py
Requires the Models_3D.py script to be in same working directory.
This is where all the population model functions are stored. 
This is written for the models specifically found in that script. 
Script will simulate the model using a fixed set of input parameters. That is,
no optimizations will occur here, and the assumption is you've searched
thoroughly for the parameter set with the highest likelihood value. The goal
of this script is to produce plots of the 3D-JSFS for the data and model,
plus the residuals. This will take the form of three two-population comparisons (2D-JSFS plots)
per plot. A pdf file will be written for each model to the current working directory.
Requires user to edit sections of code marked with #**************
You'll absolutely need to provide the path to your SNPs input file
along with your specific projections and population labels. 
###########################################################
**Note, if you see this error when plotting:
"ValueError: Data has no positive values, and therefore can not be log-scaled."
You will need to change the vmin in the plotting routine from None to something 0<1.
I have changed vmin to 0.005-0.01 with good results. So, at the bottom of this script:
dadi.Plotting.plot_3d_comp_multinom(model, data, resid_range = 3)
becomes:
dadi.Plotting.plot_3d_comp_multinom(model, data, resid_range = 3, vmin=0.005)
This should eliminate that particular error.
###########################################################
###########################################################
*Note: This may be version specific to both dadi and numpy, but
I kept getting a plotting error and so changed some lines in the
dadi Plotting.py module starting with line 395:
    ax4 = pylab.subplot(2,2,4)
    flatresid = numpy.compress(numpy.logical_not(resid.mask.ravel()), 
                               resid.ravel())
    hister, bin_edges = numpy.histogram(flatresid, bins = 30, range=None)
    #ax4.hist(flatresid, bins=20, normed=True)
    ax4.bar(bin_edges[:-1], hister, width=0.2)
    ax4.set_title('residuals')
    ax4.set_yticks([])
    
This solved an issue with numpy incompatibility (instead of using the 
ax4.hist syntax). An updated numpy should however work fine as is.
###########################################################
############################################
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-Matplotlib
-dadi
############################################
Dan Portik
daniel.portik@uta.edu
April 2017
'''
t_begin = datetime.now()

#===========================================================================
#get snps file 
#**************
snps = "/public/barratt/work/Lflavomaculatus/Inputs/dadi/dadi_3pops_north_south_taratibu_snps.txt"

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#projection sizes, in ALLELES not individuals
proj = [36,29,15]
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=['north','south','taratibu']

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs_1 = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

print '\n', '\n', "Data for spectrum:"
print "projection", proj
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'


#======================================================================================
#create function to run models with fixed user input parameters (presumably from best runs)

def Three_Pop_Plot(pts, fs, model_name, params):
    
    ######################################################################################################################
    #Basic models
    ######################################################################################################################
    
    if model_name == "split_nomig":
        
        #####################################
        #Split with no migration
        #####################################
        print "---------------------------------------------------"
        print "Split with No Migration",'\n'
        
        #first call a predefined model
        model_call = Models_3D.split_nomig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        print "starting parameters = ", params
        #simulate the model with the optimized parameters
        sim_model = func_exec(params, fs.sample_sizes, pts)

        #calculate likelihood
        ll = dadi.Inference.ll_multinom(sim_model, fs)
        ll = np.around(ll, 2)
        print "likelihood = ", ll

        #calculate theta
        theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
        theta = np.around(theta, 2)
        print "Theta = ", theta

        #calculate AIC 
        aic = ( -2*( float(ll))) + (2*6)
        print "AIC = ", aic, '\n', '\n'
        
        return sim_model

      
#======================================================================================
# Finally, execute model with appropriate arguments, including name for returned model
# returned_model = Three_Pop_Plot(pts, fs, model_name, params):

# returned_model:  name the simulated model returned by function (to store for plotting later)

# pts:  grid choice (list of three numbers, ex. [20,30,40]

# fs:  spectrum object name

# model_name:  "split_nomig", "split_symmig_all", "split_symmig_adjacent", "refugia_1",
#        "refugia_2", "refugia_3", "ancmig_3", "ancmig_2", "ancmig_1", 

# params:  list of best parameters to start optimizations from
#		ex. some_params = [4.787,0.465,7.071,1.879,0.181,0.635]


#===========================================================================
# enter best param values for each model here, presumably you will get these
# from the outputs of the previous script

#************** "split_nomig"
# 6 Values
split_nomig_params = [3.4995,1.9552,1.6082,0.3107,0.3913,0.2069]


#**************
#Input some of the basic reusable arguments here
pts = [50,60,70]
fs = fs_1
outfile = "north_south_taratibu"


#===========================================================================
# returned_model = Three_Pop_Plot(pts, fs, model_name, params):

# Each model is executed with one replicate using fixed parameter values.
# The simulated model is stored as an object to be called on in the
# subsequent plotting function.

#**************
#Models from the Models_3D.py script
split_nomig = Three_Pop_Plot(pts, fs, "split_nomig", split_nomig_params)

#======================================================================================
#write a plotting function for data and model comparison

def plot_all(sim_model, data, outfile, model_name):
    print '{0}_{1}.pdf'.format(outfile,model_name), '\n'
    outname = '{0}_{1}.pdf'.format(outfile,model_name)
    fig = pylab.figure(1)
    fig.clear()
    dadi.Plotting.plot_3d_comp_multinom(sim_model, data, resid_range = 3, vmin=0.005)
    fig.savefig(outname)
    
#======================================================================================
# plot_all(sim_model, fs, outfile, model_name)

# sim_model:  the simulated model returned from function Three_Pop_Plot

# fs:  spectrum object name (ie the "data")

# outfile:  same prefix for output naming of pdf files -> "{0}_{1}.pdf".format(outfile,model_name)

# model_name:  "split_nomig"

# call plot function for each model

plot_all(split_nomig, fs, outfile, "split_nomig")

#===========================================================================
#clock it!

t_finish = datetime.now()
elapsed = t_finish - t_begin

print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'

#===========================================================================