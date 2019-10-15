import matplotlib
matplotlib.use('Agg')
#from mpi4py import MPI
#import mpi4py
import numpy as np
import matplotlib.pyplot as plt
import time
import h5py
import sys
from scipy.integrate import simps
from scipy.interpolate import interp1d
import os


#comm = MPI.COMM_WORLD
#size = comm.Get_size()
#rank = comm.Get_rank()



def dens(r,rho0, bins, gammas,rh):
    output = np.zeros_like(r)
    bins = bins*rh
    rho0 = float(10**rho0)

    x1 = [(r < bins[0])]
    x2 = [((r >= bins[0])& (r < bins[1]))]
    x3 = [((r >= bins[1])& (r < bins[2]))]
    x4 = [((r >= bins[2])& (r < bins[3]))]
    x5 = [((r >= bins[3])& (r < bins[4]))]
    x6 = [(r >= bins[4])]


    output[x1] =  rho0 * (r[x1]/bins[0])**(-gammas[0])

    output[x2] =  rho0 * (r[x2]/bins[1])**(-gammas[1]) * (bins[1]/bins[0]) **  (-gammas[1] )

    output[x3] =  rho0 * (r[x3]/bins[2])**(-gammas[2]) * (bins[2]/bins[1]) **  (-gammas[2]) * (bins[1]/bins[0]) **  (-gammas[1] )

    output[x4] = rho0 * (r[x4]/bins[3])**(-gammas[3])* (bins[3]/bins[2]) **  (-gammas[3]) * (bins[2]/bins[1]) **  (-gammas[2]) * (bins[1]/bins[0]) **  (-gammas[1])


    output[x5] =  rho0 * (r[x5]/bins[4])**(-gammas[4]) * (bins[4]/bins[3])**(-gammas[4])* (bins[3]/bins[2]) **  (-gammas[3]) * (bins[2]/bins[1]) **  (-gammas[2]) * (bins[1]/bins[0]) **  (-gammas[1])


    output[x6] = rho0 * (r[x6]/bins[4])**(-gammas[4]) * (bins[4]/bins[3])**(-gammas[4])* (bins[3]/bins[2]) **  (-gammas[3]) * (bins[2]/bins[1]) **  (-gammas[2]) * (bins[1]/bins[0]) **  (-gammas[1])

    return output


def mass_dens(x,rho0, bins, gammas,rh):


    output = np.empty_like(x)

    bins = bins*rh


    rho0 = pow(10, rho0)* 4 * np.pi

    x1 = [(x < bins[0])]
    x2 = [((x >= bins[0])& (x < bins[1]))]
    x3 = [((x >= bins[1])& (x < bins[2]))]
    x4 = [((x >= bins[2])& (x < bins[3]))]
    x5 = [((x >= bins[3])& (x < bins[4]))]
    x6 = [(x >= bins[4])]


    first_integral = rho0 * pow(bins[0], gammas[0]) * pow(bins[0], 3 - gammas[0])/(3 - gammas[0])
    second_integral = rho0 *pow(bins[1], gammas[1]) * pow(bins[1]/bins[0], - gammas[1]) * ( pow(bins[1], 3-gammas[1]) -  pow(bins[0], 3-gammas[1])  )/(3. -gammas[1])
    third_integral = rho0 * pow(bins[2], gammas[2]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2]) * ( pow(bins[2], 3-gammas[2]) -  pow(bins[1], 3-gammas[2])  )/(3. -gammas[2])
    fourth_integral = rho0 *pow(bins[3], gammas[3]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * ( pow(bins[3], 3-gammas[3]) -  pow(bins[2], 3-gammas[3])  )/(3. -gammas[3])



    output[x1]= rho0 * pow(bins[0], gammas[0]) * pow(x[x1], 3 - gammas[0])/(3 - gammas[0])

    output[x2] =  first_integral + rho0 *pow(bins[1], gammas[1]) * pow(bins[1]/bins[0], - gammas[1]) * ( pow(x[x2], 3-gammas[1]) -  pow(bins[0], 3-gammas[1])  )/(3. -gammas[1])

    output[x3] = first_integral + second_integral + rho0 * pow(bins[2], gammas[2]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2]) * ( pow(x[x3], 3-gammas[2]) -  pow(bins[1], 3-gammas[2])  )/(3. -gammas[2])

    output[x4] =  first_integral + second_integral + third_integral + rho0 *pow(bins[3], gammas[3]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * ( pow(x[x4], 3-gammas[3]) -  pow(bins[2], 3-gammas[3])  )/(3. -gammas[3])

    output[x5] =  first_integral + second_integral + third_integral + fourth_integral + rho0 *pow(bins[4], gammas[4]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * pow(bins[4]/bins[3], - gammas[4]) * (pow( x[x5], 3-gammas[4]) -  pow(bins[3], 3-gammas[4])  )/(3. -gammas[4])

    output[x6] = first_integral + second_integral + third_integral + fourth_integral + rho0 *pow(bins[4], gammas[4]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * pow(bins[4]/bins[3], - gammas[4]) * (pow( x[x6], 3-gammas[4]) -  pow(bins[3], 3-gammas[4])  )/(3. -gammas[4])




    return output



workdir = sys.argv[1]
project_name = sys.argv[2]
first_gal = sys.argv[3]
samples = int(sys.argv[4])



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

    print galaxy

    r = np.logspace(-2,1.5, 25)

    huge_scale = np.logspace(-3,1.6, 10000)
    inter_zhao = (huge_scale[1:] + huge_scale[:-1])/2.
    diff_zhao_dist = huge_scale[1:] - huge_scale[:-1]

    
    
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
    chains = chains[:,[0,1,2,3,4,5]]
    #chains = np.array_split(chains, size)
    #chains = chains[rank]
   
    print len(chains), 'Split chains'


    bin_edges = np.array([0.125, 0.25, 0.50, 1, 2 , 4 ])*rh
    mid_bin = (np.log10(bin_edges[1:]) + np.log10(bin_edges[:-1]))/2.
    mid_bin = 10**mid_bin

    tot_bins = np.array([bin_edges[0], mid_bin[0], bin_edges[1], bin_edges[1], mid_bin[1], bin_edges[2],bin_edges[2], mid_bin[2], bin_edges[3], bin_edges[3], mid_bin[3], bin_edges[4], bin_edges[4], mid_bin[4], bin_edges[5]])


    js = [[] for i in range(len(r))]
    masses =   [[] for i in range(len(r))]
    grad_rho =   [[] for i in range(len(tot_bins))]
    betas = []

    t0 = time.time()


    #betas.append(chains[:,5])

    c = 0
    bins = np.array([0.25,0.5,1.0,2.0,4.0])



    for chain in chains:

    	rho0,gamma0,gamma1,gamma2,gamma3,gamma4= chain[0:6]
    	gammas = [gamma0,gamma1,gamma2,gamma3,gamma4]
    	big_j = dens(r, rho0, bins, gammas,rh)

    	##zhao_huge =  dens(huge_scale, rho0, bins, gammas,rh)
    	##zhao_inter = dens(inter_zhao, rho0, bins, gammas,rh)
    	#grad_zhao = (zhao_huge[1:] - zhao_huge[:-1])/diff_zhao_dist

    	#loglog = inter_zhao/zhao_inter * grad_zhao

    	#grad_prof = interp1d(inter_zhao, loglog, kind = 'linear')

    	#spec_loglog = grad_prof(mid_bin)



    	big_j = np.log10(big_j)

    	#c = c+1
    	#if c%1000 == 0:
    	#	print rank, '%d' %c + '/%d' % len(chains), big_j, spec_loglog, mid_bin

    	for b in range(0, len(big_j)):


    		#mass_try = simps(mass_zhao(lin,rhos,rs,alpha,beta,gamma), lin, even = 'avg')
    		mass_try = mass_dens(np.array([r[b]]),rho0, bins, gammas,rh)
    		js[b].append(float(big_j[b]))
    		masses[b].append(float(np.log10(mass_try)))
    	count = 0
    	for b in range(0, len(mid_bin)):
    		for b2 in range(0,3):
    			grad_rho[3*b + b2].append(-gammas[b])





    js = np.concatenate([js])
    js = js.T
    masses = np.concatenate([masses])
    masses = masses.T
    grad_rho = np.concatenate([grad_rho])
    grad_rho = grad_rho.T



    f = open(workdir + project_name + '/Analysis/Output/'  + 'Galaxy_%d' % galaxy +  '_Density.txt', 'w')
    np.savetxt(f, js, delimiter = '\t')
    f.close()
    f2 = open(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d' % galaxy +  '_Mass.txt' , 'w')
    np.savetxt(f2,masses, delimiter = '\t')
    f2.close()
    f3 = open(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d' % galaxy +  '_LogLog.txt' , 'w')
    np.savetxt(f3,grad_rho, delimiter = '\t')
    f3.close()
