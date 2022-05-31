import sys
import os
import numpy
import dadi
#import pylab
from datetime import datetime
import Optimize_Functions

'''
Usage: python dadi_Run_Optimizations.py - #### I've used only V5 here, not V1-4 examples! #####

This is meant to be a general use script to run dadi to fit any model on an
afs/jsfs with one to three populations. The user will have to edit information about
their allele frequency dadi.Spectrum and provide a custom model. The instructions are
annotated below, with a #************** marking sections that will have to be edited.
Several examples of how to use various arguments to control optimizations are shown.

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary.

If you'd like to use this script for larger sets of models already available, please
look on the other repositories to see how to import models from external model scripts.


General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.

 
Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688


Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject. 

 The chi-squared test statistic is calculated assuming the sfs is folded. If this is not
 true for your data set, this number will not be accurate. This could be edited in the
 'collect_results' function in the Optimize_Functions.py script for an unfolded dadi.Spectrum.

Citations:
 If you use these scripts for your work, please cite the following publication:
    Portik, D.M., Leach, A.D., Rivera, D., Blackburn, D.C., Rdel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K.Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 52455263.
    doi: 10.1111/mec.14266

-------------------------
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated May 2018
'''

#===========================================================================
# Import data to create joint-site frequency dadi.Spectrum
#===========================================================================

#**************
snps = "/public/barratt/work/Lflavomaculatus/Inputs/dadi/dadi_3pops_north_south_taratibu_snps.txt"

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=['north','south','taratibu']

#**************
#projection sizes, in ALLELES not individuals
proj = [36,29,15]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded dadi.Spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#print some useful information about the afs or jsfs
print "\n\n============================================================================\nData for site frequency dadi.Spectrum\n============================================================================\n"
print "projection", proj
print "sample sizes", fs.sample_sizes
sfs_sum = numpy.around(fs.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

#================================================================================
# Here is an example of using a custom model within this script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script:
 
 Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" ")
 
   Mandatory Arguments =
    fs:  dadi.Spectrum object name
    pts: grid size for extrapolation, list of three values
    outfile:  prefix for output naming
    model_name: a label to slap on the output files; ex. "no_mig"
    func: access the model function from within 'dadi_Run_Optimizations.py' or from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)

   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order
'''


#================================================================================
# Let's start by defining our models!
#================================================================================

##########################################################################################
#Basic models of (no gene flow / gene flow) between (all / some) population pairs
##########################################################################################
def split_no_mig(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    """
    #6 parameters	
    nu1, nuA, nu2, nu3, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def split_symmig_all(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def split_symmig_adjacent(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3. Assume 2 occurs
    in between populations 1 and 3, which do not come in to contact with one another.
    Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
    """
    #9 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def refugia_1(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, gene flow does not occur. Period of symmetric secondary contact occurs between 
    adjacent populations (ie 1<->2, 2<->3, but not 1<->3) after all splits are complete. 
    longest isolation
    """
    #9 parameters
    nu1, nuA, nu2, nu3, m1, m2, T1, T2, T3 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = dadi.Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def refugia_2(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, with gene flow. After appearance of 2 and 3, gene flow also occurs between 1 
    and 2.
    shorter isolation
    """
    #8 parameters
    nu1, nuA, nu2, nu3, m1, m2, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def refugia_3(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur, but then 
    secondary contact occurs. Split between pops 2 and 3 occurs with gene flow, and gene flow
    happens between 1 and 2 as well.
    shortest isolation
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1a, T1b, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=0, m21=0)
    phi = dadi.Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def ancmig_3(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), with gene flow, which then stops. Split 
    between pops 2 and 3, gene flow does not occur at all.
    longest isolation
    """
    #8 parameters
    nu1, nuA, nu2, nu3, mA, T1a, T1b, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    phi = dadi.Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def ancmig_2(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3, and all gene flow ceases.
    shorter isolation
    """
    #7 parameters
    nu1, nuA, nu2, nu3, mA, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def ancmig_1(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3 with gene flow, then all gene flow ceases.
    shortest isolation
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2, T3 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    phi = dadi.Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def sim_split_no_mig(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3, gene flow does not occur. 
    """
    #4 parameters
    nu1, nu2, nu3, T1 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def sim_split_sym_mig_all(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    """
    #7 parameters
    nu1, nu2, nu3, m1, m2, m3, T1 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def sim_split_sym_mig_adjacent(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3
    Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
    """
    #6 parameters
    nu1, nu2, nu3, m1, m2, T1 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs
    
def sim_split_refugia_sym_mig_all(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3 followed by isolation. Period of symmetric secondary contact occurs between 
    all populations after all splits are complete. 
 
    """
    #8 parameters
    nu1, nu2, nu3, m1, m2, m3, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def sim_split_refugia_sym_mig_adjacent(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3 followed by isolation. Period of symmetric secondary contact occurs between 
    adjacent populations (ie 1<->2, 2<->3, but not 1<->3) after all splits are complete. 
 
    """
    #7 parameters
    nu1, nu2, nu3, m1, m2, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def split_nomig_human(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair, size change
    """
    #10 parameters	
    nu1a, nuA, nu2a, nu3a, nu1b, nu2b, nu3b, T1, T2, T3 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1a, nu2a=nuA, m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1a, nu2a, nu3a, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = dadi.Integration.three_pops(phi, xx, T3, nu1b, nu2b, nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs
   
'''
Example 5. It's also good run the optimization routine multiple times. Let's write a short
loop to do the above optimization routine five times. We will name the prefix based
on which point we are at, and include it within the loops. Note that when you use
the range argument in python it will go up to, but not include, the final number.
That's why I have written a range of 1-6 to perform this 5 times.
'''
pts = [50,60,70]
p_labels = "nu1, nuA, nu2, nu3, T1, T2"
upper = [20,20,20,20,15,15]
lower = [0.01,0.01,0.01,0.01,0.1,0.1]
reps = [5,7,10]
maxiters = [3,5,5]
folds = [3,2,1]

for i in range(1,2):
    prefix = "full_3D_dadi_job_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_no_mig", split_no_mig, 3, 6,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2"
#upper = [20,20,20,20,20,20,20,20,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_symmig_all", split_symmig_all, 3, 10,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2"
#upper = [20,20,20,20,20,20,20,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_symmig_adjacent", split_symmig_adjacent, 3, 9,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nuA, nu2, nu3, m1, m2, T1, T2, T3"
#upper = [20,20,20,20,20,20,15,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "refugia_1", refugia_1, 3, 9,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nuA, nu2, nu3, m1, m2, T1, T2"
#upper = [20,20,20,20,20,20,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "refugia_2", refugia_2, 3, 8,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, T1a, T1b, T2"
#upper = [20,20,20,20,20,20,20,15,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "refugia_3", refugia_3, 3, 10,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nuA, nu2, nu3, mA, T1a, T1b, T2"
#upper = [20,20,20,20,20,15,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "ancmig_3", ancmig_3, 3, 8,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nuA, nu2, nu3, mA, T1, T2"
#upper = [20,20,20,20,20,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "ancmig_2", ancmig_2, 3, 7,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2, T3"
#upper = [20,20,20,20,20,20,20,15,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "ancmig_1", ancmig_1, 3, 10,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nu2, nu3, T1"
#upper = [20,20,20,15]
#lower = [0.01,0.01,0.01,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_no_mig", sim_split_no_mig, 3, 4,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, T1, T2"
#upper = [20,20,20,20,20,20,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_no_mig_size", sim_split_no_mig_size, 3, 8,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nu2, nu3, m1, m2, m3, T1"
#upper = [20,20,20,20,20,20,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_sym_mig_all", sim_split_sym_mig_all, 3, 7,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nu2, nu3, m1, m2, T1"
#upper = [20,20,20,20,20,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_sym_mig_adjacent", sim_split_sym_mig_adjacent, 3, 6,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nu2, nu3, m1, m2, T1, T2"
#upper = [20,20,20,20,20,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_refugia_sym_mig_adjacent", sim_split_refugia_sym_mig_adjacent, 3, 7,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)

#pts = [50,60,70]
#p_labels = "nu1, nu2, nu3, m1, m2, m3, T1, T2"
#upper = [20,20,20,20,20,20,15,15]
#lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1]
#reps = [10,20,50]
#maxiters = [5,10,20]
#folds = [3,2,1]

#for i in range(1,2):
#    prefix = "full_3D_dadi_job_{}".format(i)
#    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_refugia_sym_mig_all", asym_mig, 3, 8,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)





