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

huge_scale = np.logspace(-3.01,1.6, 10000)
inter_zhao = (huge_scale[1:] + huge_scale[:-1])/2.
diff_zhao_dist = huge_scale[1:] - huge_scale[:-1]



for galaxy in galaxies:


	print galaxy

	r = np.logspace(-2,1.5, 25)











	try:
		chains = np.loadtxt(workdir + project_name + '/Galaxy_%d/' % galaxy + 'Galaxy_%d' %galaxy + '_Chains' + project_name + ".txt")
	except IOError:
		print 'Not here'
		continue

	orglen = (len(chains)/1000) / 2

	chains = chains[-1000*orglen:]


	rhdat = np.loadtxt(workdir + '/GalaxyData/Galaxy_%d_Rhalf.txt' % int(galaxy))

    	rh = float(rhdat)


	#above_rs, = np.where(chains[:,1] >= np.log10(rh))
	#below_gamma, = np.where(chains[:,4] <= 1.)
	#together = np.intersect1d(above_rs, below_gamma)

	#chains = chains[together]


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
	chains = chains[:,[0,1,2,3,4]]
	ch_org = len(chains)

	
	


	js = [[] for i in range(len(r))]
	masses =   [[] for i in range(len(r))]
	grad_rho =   [[] for i in range(len(r))]
	betas = []

	t0 = time.time()


	#betas.append(chains[:,5])

	c = 0




	for chain in chains:

		rhos,rs,alpha,beta,gamma = chain[0:5]
		big_j = zhao(r, rhos,rs,alpha,beta,gamma)

		zhao_huge =  zhao(huge_scale, rhos,rs,alpha,beta,gamma)
		zhao_inter = zhao(inter_zhao, rhos,rs,alpha,beta,gamma)
		grad_zhao = (zhao_huge[1:] - zhao_huge[:-1])/diff_zhao_dist

		loglog = inter_zhao/zhao_inter * grad_zhao

		grad_prof = interp1d(inter_zhao, loglog, kind = 'linear')

		spec_loglog = grad_prof(r)



		big_j = np.log10(big_j)

		c = c+1
		

		for b in range(0, len(big_j)):

			lin = np.logspace(-3,np.log10(r[b]),2000)
			mass_try = simps(mass_zhao(lin,rhos,rs,alpha,beta,gamma), lin, even = 'avg')
			js[b].append(float(big_j[b]))
			masses[b].append(float(np.log10(mass_try)))
			grad_rho[b].append(float(spec_loglog[b]))

	js = np.concatenate([js])
	js = js.T
	masses = np.concatenate([masses])
	masses = masses.T
	grad_rho = np.concatenate([grad_rho])
	grad_rho = grad_rho.T



	f = open(workdir + project_name + '/Analysis/Output/'  + 'Galaxy_%d' % galaxy +  '_Density2.txt' , 'w')
	np.savetxt(f, js, delimiter = '\t')
	f.close()
	f2 = open(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d' % galaxy +  '_Mass2.txt' , 'w')
	np.savetxt(f2,masses, delimiter = '\t')
	f2.close()
	f3 = open(workdir + project_name + '/Analysis/Output/'  + 'Galaxy_%d' % galaxy +  '_LogLog2.txt' , 'w')
	np.savetxt(f3,grad_rho, delimiter = '\t')
	f3.close()
