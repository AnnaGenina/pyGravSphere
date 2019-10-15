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






workdir = sys.argv[1]
project_name = sys.argv[2]
first_gal = sys.argv[3]
col = int(sys.argv[4])
samples = int(sys.argv[5])



galaxies = [int(first_gal)]


for galaxy in galaxies:


	print galaxy

	r = np.logspace(-2,1.5, 25)


	try:
		chains = np.loadtxt(workdir + project_name + '/Galaxy_%d/' % galaxy + 'Galaxy_%d' %galaxy + '_Chains' + project_name + ".txt", )
	except IOError:
		continue

	orglen = (len(chains)/1000) / 2
	chains = chains[-1000*orglen:]

	min_chisq = np.max(chains[:,-1])
    	index, = np.where(chains[:,-1] > min_chisq*10.0)
   
    
        chains = chains[index]

	
	ch_org = len(chains)


	like = chains[:,-1]
	fin, = np.where(np.isfinite(like) == True)	
	chains = chains[fin]

	

	np.random.shuffle(chains)
        num_ch = np.random.choice(np.arange(len(chains)),samples, replace = False)
        chains = chains[num_ch]




	betas = [[] for i in range(len(r))]
	betas = np.zeros((len(chains), len(r)))


	for i in range(0, len(r)):
		betas[:,i] = chains[:,col]




	f_handle = file(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_Beta'%galaxy + '.txt' , 'w')
	np.savetxt(f_handle,betas, delimiter = '\t')
       	f_handle.close()
