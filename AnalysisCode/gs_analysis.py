import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import time
import h5py
import sys
from scipy.integrate import simps
from scipy.interpolate import interp1d
import os
import analysis_func




workdir = sys.argv[1]
codedir = sys.argv[2]
project_name = sys.argv[3]
galaxy = sys.argv[4]
samples = int(sys.argv[5])
cut_off = int(sys.argv[6])
chi_cut = float(sys.argv[7])
min_r = float(sys.argv[8])
max_r = float(sys.argv[9])
points = int(sys.argv[10])


print 'Starting'

foptions = open(workdir + project_name + '/options.txt', 'r')
input_opt = (foptions.readline()).split()
dm_option = input_opt[0] 
beta_option = input_opt[1]
plummer_option = input_opt[2]
nsteps = int(input_opt[3])
nwalkers = int(input_opt[4])
options = [dm_option, beta_option, plummer_option]

priors = np.loadtxt(workdir + project_name + '/Submissions/priors.txt', dtype = 'str')

print 'Reading output chains. This might take a while.'

try:
	chains = np.loadtxt(workdir + '/' + project_name + '/%s' % galaxy +'/%s_Chains' %galaxy + project_name + ".txt")
except IOError:
	print "Chains were not found. Quitting."
	quit()

print 'Finished reading chains'


print 'Leaving the last %d iterations' % cut_off
chains = chains[cut_off*nwalkers:]


print 'Selecting better chi squared chains'
min_chisq = np.max(chains[:,-1])
index, = np.where(chains[:,-1] > min_chisq*chi_cut)
chains = chains[index]
print '%d chains remain' %len(chains)

ch_org = len(chains)
like = chains[:,-1]
fin, = np.where(np.isfinite(like) == True)	
chains = chains[fin]
if len(chains) < samples:
	print 'Less finite likelihood chains than required samples.'
	print 'There are %d' %len(chains) + ' finite chains. Proceeding.'
	samples = len(chains)

print 'Selecting random sample of chains'
np.random.shuffle(chains)
num_ch = np.random.choice(np.arange(len(chains)),samples, replace = False)
chains = chains[num_ch]
ch_org = len(chains)

print 'Running beta'
analysis_func.return_beta(chains, options, priors, min_r, max_r, points,codedir,workdir,project_name, galaxy)
print 'Running Mass/Density/Slope'
analysis_func.return_mass(chains, options, priors, min_r, max_r, points,codedir,workdir,project_name,galaxy)
print 'Running Plummer'
analysis_func.return_plummer(chains, options, priors, min_r, max_r, points,codedir,workdir,project_name,galaxy)
print 'Running sigR and VSPs'
analysis_func.return_sigma_vsp(chains, options, priors, min_r, max_r, points,codedir,workdir,project_name,galaxy)
print 'Finished, hurray!'




