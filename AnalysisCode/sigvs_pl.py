import numpy as np
from mpi4py import MPI
import mpi4py
import sys
sys.path.append("/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/Read/")
from GravSphere import gal_input
import simps
from scipy.optimize import minimize
import time








workdir = sys.argv[1]
project_name = sys.argv[2]
first_gal = sys.argv[3]
samples = int(sys.argv[4])

galaxies = [int(first_gal)]



def plummer(x,norm,rc):
        return np.log10(   (1./(np.pi * rc * rc))    / (  1 + x*x/(rc*rc)  )**2. ) + norm

def min_plummer(args,x_data, y_data,err_data):
        norm,rc = args
        return np.sum(  (  10**plummer(x_data,norm,rc) - y_data  )**2. / err_data**2. )


for galaxy in galaxies:

    
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

    

    kindat,lightpower,vir_shape,surfden,r_c, stellar_mass = gal_input.galaxy_data_read(galaxy, workdir + '/GalaxyData/')


    min_rad = np.log10(r_c/100)
    max_rad = np.log10(r_c*50)

    vs1arr = []
    vs2arr = []
    res_arr = []


    n_params = len(chains[0])
    
    chlen = len(chains)

    print 'N params', n_params
    print 'n chains', chlen
	
    plum_fit =minimize(min_plummer, [np.log10(np.amax(surfden[:,1])),0.5], args=(surfden[:,0], surfden[:,1],surfden[:,2],), bounds = ((0, 10),(0.1,5),))
    m1,a1 = plum_fit.x

    count = 0  	
    t0 = time.time()	
    for chain in chains:


	count = count + 1
	if count%100 == 0:
		print count, chlen, time.time() - t0
		t0 = time.time()

	if n_params == 18:
		rho0, gamma0, gamma1, gamma2, gamma3, gamma4, beta0, betainf, ra, eta, m1,a1,m2,a2,m3,a3,mstar,like = chain
	if n_params == 15:
		rho0, gamma0, gamma1, gamma2, gamma3, gamma4, beta0, m1,a1,m2,a2,m3,a3,mstar,like = chain
		betainf = beta0
		ra = 0
		eta = 1
	if n_params == 9:
		rho0, gamma0, gamma1, gamma2, gamma3, gamma4, beta0,mstar,like = chain
		betainf = beta0
		ra = 0
		eta = 1
		m2 = -5
		a2 = 1
		m3 = -5
		a3 = 1
	if n_params == 12:
		rho0, gamma0, gamma1, gamma2, gamma3, gamma4, beta0,betainf, ra, eta,mstar,like = chain
		m2 = -5
		a2 = 1
		m3 = -5
		a3 = 1
	

        rho_params = np.array([r_c, rho0, gamma0,gamma1,gamma2,gamma3,gamma4 ])
        beta_params = np.array([ beta0, betainf  ,10**ra, eta])
        plum_params = np.array([10**m1,a1,10**m2,a2,10**m3,a3])
        

	results = np.zeros_like(kindat[:,0])
	vsp1 = np.zeros((1), dtype = np.float64)
        vsp2 = np.zeros((1), dtype = np.float64)

	simps.PowerLawFitfunc(kindat[:,0], kindat[:,1], kindat[:,2], rho_params, beta_params, plum_params, vir_shape,10**mstar, r_c, 100, results, vsp1, vsp2, min_rad, max_rad)

	vs1arr.append(float(vsp1))
	vs2arr.append(float(vsp2))
	res_arr.append(results)


    vs1arr = np.array(vs1arr)
    vs2arr = np.array(vs2arr)
    res_arr = np.concatenate([res_arr])


    f1 = open(workdir + project_name + '/Analysis/Output/'  + 'Galaxy_%d' % galaxy +  '_SigLos.txt', 'w')
    f2 = open(workdir + project_name + '/Analysis/Output/'  + 'Galaxy_%d' % galaxy +  '_VSP1.txt', 'w')
    f3 = open(workdir + project_name + '/Analysis/Output/'  + 'Galaxy_%d' % galaxy +  '_VSP2.txt', 'w')

    np.savetxt(f2, vs1arr, delimiter = '\t')
    np.savetxt(f3, vs2arr, delimiter = '\t')
    np.savetxt(f1, res_arr, delimiter = '\t')




