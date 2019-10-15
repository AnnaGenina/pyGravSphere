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
starter = int(sys.argv[4])
samples = int(sys.argv[5])



galaxies = [int(first_gal)]


def plummer_proj_comp(x,mass,shape):

	return 10**mass/(np.pi * (10**shape)**2.) * (1. + x**2./(10**shape)**2)**(-2.)
def plummer_proj_sum(args, x_data, n_comp):

	sums_proj = 0

	start = 0
	for i in range(0, n_comp):
		sums_proj = sums_proj + plummer_proj_comp(x_data,args[start], args[start + 1])

		start = start + 2

	return np.log10(sums_proj)


for galaxy in galaxies:

	rhdat = np.loadtxt(workdir + '/GalaxyData/Galaxy_%d_Rhalf.txt' % int(galaxy))

    	rh = float(rhdat)


	r = np.logspace(-2,1.5, 25)

	try:

		chains = np.loadtxt(workdir + project_name + '/Galaxy_%d/' % galaxy + 'Galaxy_%d' %galaxy + '_Chains' + project_name + ".txt")
	except IOError:
		print 'Not here'
		continue

	chains = chains[-2000*1000:]
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
	chains = chains[:,[starter,starter+1,starter+2,starter +3,starter+4,starter+5]]
	ch_org = len(chains)



	


	plum_l = [[] for i in range(len(r))]

	t0 = time.time()


	c = 0
	bins = np.array([0.25,0.5,1.0,2.0,4.0])



	for chain in chains:

		m1,a1,m2,a2,m3,a3 = chain[0:6]
		a1 = np.log10(a1)
		a2 = np.log10(a2)
		a3 = np.log10(a3)



		for b in range(0, len(r)):
			plum_l[b].append(plummer_proj_sum([m1,a1,m2,a2,m3,a3], r[b], 3))




	plum_l = np.concatenate([plum_l])
	plum_l = plum_l.T

	f = open(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d' % galaxy +  '_Plummer.txt' , 'w')
	np.savetxt(f, plum_l, delimiter = '\t')
	f.close()
