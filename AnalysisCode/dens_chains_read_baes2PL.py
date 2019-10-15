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
samples = int(sys.argv[4])



galaxies = [int(first_gal)]

def zhao(x, rhos, rs,alpha,beta,gamma):
	return 10**rhos * (x/(10**rs))**(-gamma) * (1 + (x/(10**rs))**alpha)**(-(beta-gamma)/alpha)
def mass_zhao(x,rhos, rs,alpha,beta,gamma):
	return 4. * np.pi * x**2. * zhao(x,rhos, rs,alpha,beta,gamma)
def baes(x, b0, binf, ra, eta):
	return (b0 + binf * (x/ra)**eta) / ( 1 + (x/ra)**eta)


huge_scale = np.logspace(-3.01,1.6, 10000)
inter_zhao = (huge_scale[1:] + huge_scale[:-1])/2.
diff_zhao_dist = huge_scale[1:] - huge_scale[:-1]





for galaxy in galaxies:

	r = np.logspace(-2,1.5, 25)

	try:
		chains = np.loadtxt(workdir + project_name + '/Galaxy_%d/' % galaxy + 'Galaxy_%d' %galaxy + '_Chains' + project_name + ".txt")
	except IOError:
		continue

	print 'Got chains'
	
	chains = chains[-2000*1000:]

	orglen = (len(chains)/1000) / 2
	print 'Reducing to ', orglen

	chains = chains[-1000*orglen:]
	
	min_chisq = np.max(chains[:,-1])
    	index, = np.where(chains[:,-1] > min_chisq*10.0)
   
    
        chains = chains[index]

	
	ch_org = len(chains)

	print 'Best chi', ch_org

	like = chains[:,-1]
	fin, = np.where(np.isfinite(like) == True)	
	chains = chains[fin]

	

	np.random.shuffle(chains)
        num_ch = np.random.choice(np.arange(len(chains)),samples, replace = False)
        chains = chains[num_ch]


	
	chains = chains[:,[6,7,8,9]]
	ch_org = len(chains)

	print 'Beginning with samples: ', ch_org


	masses =   [[] for i in range(len(r))]

	betas = [[] for i in range(len(r))]

	t0 = time.time()

	chains[:,0] = 2.*chains[:,0] / (1. + chains[:,0]) #converte to actual values
	chains[:,1] = 2.*chains[:,1] / (1. + chains[:,1])
	chains[:,2] = 10**(chains[:,2])



	c = 0


	for chain in chains:

		b0,binf,ra,eta = chain

		bet_prof = baes(r, b0,binf,ra,eta)

		for b in range(0, len(r)):

			betas[b].append(float(bet_prof[b]))

	print 'Done'

	betas = np.concatenate([betas])
	betas = betas.T
	f_handle = file(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_Beta'%galaxy + '.txt' , 'w')
	np.savetxt(f_handle,betas, delimiter = '\t')
       	f_handle.close()
