import os
import time
import sys
import numpy as np

#Baes + PL + Plummer_Sum + VSP

# Just removed the physical constraint on montonically decreasing density slopes
workdir = sys.argv[1]
codedir = sys.argv[2]
project_name = sys.argv[3]
galaxy_number = str(sys.argv[4])
cores = int(sys.argv[5])
num_walkers = int(sys.argv[6])
burn_in = int(sys.argv[7])
steps = int(sys.argv[8])
bins = int(sys.argv[9])
darkmatter = sys.argv[10]
anisotropy = sys.argv[11]
vsps = sys.argv[12]
plummer = sys.argv[13]




with open(workdir + "/" + project_name + "/Submissions/" + "script_bin_%s" %project_name + "_%s" %(galaxy_number) +  ".py", "w") as f:
	f.write(r"""import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy import integrate
from scipy.interpolate import interp1d
import emcee
import corner
from scipy import special
import time
from schwimmbad import MPIPool
""")
	f.write(r'sys.path.append("%s")' %codedir + '\n')
	f.write(r"""from GSpro import gal_input
import gravsphere
from scipy.optimize import minimize


def plummer(x,norm,rc):
	return np.log10(   (1./(np.pi * rc * rc))    / (  1 + x*x/(rc*rc)  )**2. ) + norm

def min_plummer(args,x_data, y_data,err_data):
	norm,rc = args
	return np.sum(  (  10**plummer(x_data,norm,rc) - y_data  )**2. / err_data**2. )


def check_beta(beta):
	if beta < -0.97:
		beta = -0.97
	if beta > 0.98:
		beta = 0.98
	return beta

""")
	f.write(r"galaxy_number = '%s'" % (galaxy_number) + "\n")
	
	f.write(r"burn_in = %d" %burn_in + "\n")
	f.write(r"steps = %d" %steps + "\n")
	f.write(r"bins = %d" %bins + "\n")
	f.write(r"workdir = '%s'" % workdir + "\n")
	f.write(r"project_name = '%s'" % project_name + "\n")
	f.write(r"completed = 0" +  "\n")
        f.write(r"""kindat,lightpower,vir_shape,surfden,r_c, stellar_mass = gal_input.galaxy_data_read(galaxy_number, workdir + '/GalaxyData/')


r_bins = np.array([0.25,0.5,1,2,4]) * r_c

min_rad = np.log10(r_c/100)
max_rad = np.log10(r_c*50)

gamsmooth = 1



""")
	priors = np.loadtxt(workdir + "/" + project_name + "/Submissions/priors.txt", dtype = "str")
	for p in range(0, len(priors)):
		f.write(priors[p,1] + ',' + priors[p,2] + ',' + priors[p,3] +' = ' + priors[p,4] + ',' + priors[p,5] + ',' + priors[p,6] + '\n')

	const_var, = np.where(priors[:,6] == 'True')
	noconst_var, = np.where(priors[:,6] == 'False')
	f.write('priors = np.loadtxt(workdir + "/" + project_name + "/Submissions/priors.txt", dtype = "str")' + '\n')
	

	
	f.write(r"""data = (kindat, lightpower,vir_shape,r_c, surfden,stellar_mass)





G = 4.302 * 10.**(-6.)

def plummer_proj_comp(x,mass,shape):
	return 10**mass/(np.pi * (10**shape)**2.) * (1. + x**2./(10**shape)**2)**(-2.)

def plummer_proj_sum(args, x_data, n_comp):
	sums_proj = 0
	start = 0
	for i in range(0, n_comp):
		sums_proj = sums_proj + plummer_proj_comp(x_data,args[start], args[start + 1])
		#sums_3d = sums_3d + plummer_dens_comp(x_data,args[start], args[start + 1])
		start = start + 2
	return sums_proj""" + '\n')

	num_params = 0
	free_param_list = ''


#gRho0,gGamma0,gGamma1,gGamma2,gGamma3,gGamma4,gBeta0, gBeta1,gRa, gEta,gM1,gA1,gM2,gA2,gM3,gA3

	if darkmatter == 'PL':
		priors_dark = priors[num_params:num_params+6]
		notconst, = np.where(priors_dark[:,6] == 'False') 
		params_dark = ''
		for i in range(0, len(notconst)):
			params_dark = params_dark + priors_dark[notconst[i],0] + ','
		num_params = num_params + 6
		
	elif darkmatter == 'Zhao':
		priors_dark = priors[num_params:num_params+5]
		notconst, = np.where(priors_dark[:,6] == 'False') 
		params_dark = ''
		for i in range(0, len(notconst)):
			params_dark = params_dark + priors_dark[notconst[i],0] + ','
		num_params = num_params + 5
		
	if anisotropy == 'Baes':
		priors_anis = priors[num_params:num_params+4]
		notconst, = np.where(priors_anis[:,6] == 'False') 
		params_anis = ''
		for i in range(0, len(notconst)):
			params_anis = params_anis + priors_anis[notconst[i],0] + ','
		num_params = num_params + 4
		
	elif anisotropy	== 'Const':
		priors_anis = priors[num_params:num_params+1].reshape(6)
		 
		params_anis = ''
		if priors_anis[6] == 'False':
		
			params_anis = params_anis + priors_anis[0] + ','
		num_params = num_params + 1
		
	if plummer == 'Plummer':
		params_plummer = ''

	elif plummer == 'Plummer3':
		priors_plummer = priors[num_params:num_params+6]
		notconst, = np.where(priors_plummer[:,6] == 'False') 
		params_plummer = ''
		for i in range(0, len(notconst)):
			params_plummer = params_plummer + priors_plummer[notconst[i],0] + ','
		num_params = num_params + 6
		
	if vsps == 'NOVSP':
		f.write("vsps = np.zeros((4,))" + '\n')

	if priors[-1,6] != 'True':
		params_mass = 'mstar'
	else:
		params_mass = ''
	
	num_params = num_params + 1

	
	if plummer == 'Plummer':

		f.write("res = minimize(min_plummer, [np.log10(np.amax(surfden[:,1])),0.5], args=(surfden[:,0], surfden[:,1],surfden[:,2],), bounds = ((0, 10),(0.1,5),))" + '\n')
		f.write("m1,a1 = res.x" + '\n')
		f.write("print 'Plummer fit', m1, a1" + '\n')
		
	f.write(r'def lnlike(params):' + '\n')
	for i in range(0, len(const_var)):
		f.write("\t" + priors[const_var,0] + '=' + priors[const_var,4] + '\n')
	f.write("\t" + params_dark + params_anis + params_plummer + params_mass +  " = params" + '\n')

	if anisotropy == 'Const':
		f.write("\t" + """beta_t0 = check_beta(beta_t0)

""")
	else:
		f.write("\t" + "beta_t0 = check_beta(beta_t0)" + '\n')
		f.write("\t" + "beta_t1 = check_beta(beta_t1)" + '\n')


	if anisotropy == 'Const':
		f.write("\t" + "beta0 = (2. * beta_t0)/(1. + beta_t0)"  + '\n' + '\t' + "kindat,lightpower, vir_shape,r_c,surfden,stellar_mass = data" + '\n')
	else:
		f.write("\t" + "beta0 = (2. * beta_t0)/(1. + beta_t0)" + '\n' + '\t' + "beta1 = (2. * beta_t1)/(1. + beta_t1)" + '\n' + '\t' + "kindat,lightpower, vir_shape,r_c,surfden,stellar_mass = data" + '\n')

	if vsps == 'NOVSP':
		f.write("\t" + "vir_shape = [0,1,0,1]" + "\n")

	if darkmatter == "PL":
		f.write("\t" + "rho_param = [r_c, rho0, gamma0, gamma1, gamma2, gamma3,gamma4]" + "\n")
	else:
		f.write("\t" + "rho_param = [r_c, rhos, rs,alpha,beta,gamma]" + "\n")
	if plummer == 'Plummer':
		f.write("\t" + "plum_param = [10**m1,a1,0,1,0,1]" + "\n")
	else:
		f.write("\t" + "plum_param = [10**m1,a1,10**m2,a2,10**m3,a3]" + "\n")
	if anisotropy == 'Const':
		f.write("\t" + "beta_param = [beta0,beta0,1,0]" + "\n")
	else:	
		f.write("\t" + "beta_param = [beta0,beta1,10**ra,eta]" + "\n")

	if darkmatter == "PL":
		if plummer == 'Plummer':
			f.write("\t" + "like = gravsphere.GetLikeBinAnisSphVSPPowerLawfunc(kindat[:,0],kindat[:,1],kindat[:,2],rho_param,beta_param,plum_param,vir_shape, 10**mstar, r_c,bins, min_rad, max_rad)" + '\n')
		else:
			f.write("\t" + "like = gravsphere.GetLikeBinAnisSphVSPPowerLawfunc(kindat[:,0],kindat[:,1],kindat[:,2],rho_param,beta_param,plum_param,vir_shape, 10**mstar, r_c,bins,min_rad, max_rad)" + '\n')
	else:
		if plummer == 'Plummer':
			f.write("\t" + "like = gravsphere.GetLikeBinAnisSphVSPfunc(kindat[:,0],kindat[:,1],kindat[:,2],rho_param,beta_param,plum_param,vir_shape, 10**mstar, r_c,bins, min_rad, max_rad)" + "\n")
		else:
			f.write("\t" + "like = gravsphere.GetLikeBinAnisSphVSPfunc(kindat[:,0],kindat[:,1],kindat[:,2],rho_param,beta_param,plum_param,vir_shape, 10**mstar, r_c,bins, min_rad, max_rad)" + "\n")

	if plummer == 'Plummer3':	    
		f.write("\t" + "like = like -0.5 * np.sum((plummer_proj_sum([m1,np.log10(a1),m2,np.log10(a2),m3,np.log10(a3)],surfden[:,0],3) - surfden[:,1] )**2./surfden[:,2]**2. )" + '\n' )
		
	f.write("\t" + "if np.isnan(like) == True:" + "\n")		
	f.write("\t" + "\t" + "return -np.inf" + "\n")    			    	
	f.write("\t" + "return like" + "\n" + "\n" + "\n")    			    
	


	f.write("def lnprior(params):" + '\n')
	f.write("\t" + params_dark + params_anis + params_plummer + params_mass + " = params" + '\n')
	f.write("\t" + "kindat,lightpower,vir_shape,r_c,surfden,stellar_mass = data" + '\n')

	if darkmatter == 'PL':
		f.write("\t" + "gammas = np.array([gamma0,gamma1,gamma2,gamma3,gamma4])" + '\n')
		
		f.write("\t" + "min_gammas = np.array([gamma0min,gamma1min,gamma2min,gamma3min,gamma4min])" + '\n')
		f.write("\t" + "max_gammas = np.array([gamma0max,gamma1max,gamma2max,gamma3max,gamma4max])" + '\n')
		f.write("\t" + "min_gammas[1:] =  params[1:5] - 0.01" + "\n")
		f.write("\t" + "max_gammas[1:] =  params[1:5] + gamsmooth" + "\n")
		f.write("\t" + "min_gammas[min_gammas < min_gammas0] = min_gammas0[min_gammas < min_gammas0]" + "\n")
		f.write("\t" + "max_gammas[max_gammas > max_gammas0] = max_gammas0[max_gammas > max_gammas0]" + "\n")



	if darkmatter == 'PL':
		prior_dm = "(rho0min < rho0 < rho0max) and all(minarr < thetau < maxarr for minarr,thetau,maxarr in zip(min_gammas,gammas,max_gammas)) "


	elif darkmatter == 'Zhao':
		prior_dm = "(rhomin < rhos < rhomax) and (rsmin < rs < rsmax) and (alphamin < alpha < alphamax ) and  (betamin < beta < betamax) and  (gammamin < gamma < gammamax) "


	if anisotropy == 'Baes':
		prior_anis = "(beta0min < beta_t0 < beta0max) and (betainfmin < beta_t1 < betainfmax) and (ramin < ra < ramax) and (etamin < eta < etamax) "
	elif anisotropy	== 'Const':
		prior_anis = "(beta0min < beta_t0 < beta0max) "


	if plummer == 'Plummer3':
		prior_plummer = "and (m1min < m1 < m1max)  and (m2min < m2 < m2max)  and (m3min < m3 < m3max) and (a1min < a1 < a1max) and (a2min < a2 < a2max) and (a3min < a3 < a3max) "
		
	elif plummer == 'Plummer':
		prior_plummer = ""

	prior_mass = " and (mstarmin < mstar < mstarmax)"


	f.write("\t" + "if " + prior_dm + "and" + prior_anis + prior_plummer + prior_mass + ":" +  "\n")	
	f.write("\t" + "\t" + "return 0.0" + "\n")
	f.write("\t" + "return -np.inf" + '\n')


	f.write("""def logprob(params):
	lp = lnprior(params)

	    
	if not np.isfinite(lp):
		return -np.inf
	    
	like = lnlike(params)

	return like + lp""" + '\n')

	f.write("""
with MPIPool() as pool:
	if not pool.is_master():
		pool.wait()
		sys.exit(0)


	size = pool.size
	print size + 1


	t_org = time.time()
""")
	all_params = params_dark + params_anis + params_plummer + params_mass
	num_params = len([x.strip() for x in all_params.split(',')])
	f.write("\t" + "ndim, nwalkers = {:d}, {:d}".format(num_params, num_walkers) + '\n')

	f.write("\t" + "pos = np.zeros((nwalkers,ndim))")
	if darkmatter == 'PL':
		f.write("""
	partrack = 6
	pos[:,0] = np.random.uniform(rho0min, rho0max, nwalkers)
	pos[:,1] = np.random.uniform(gamma0min,3, nwalkers)
	pos[:,2] = np.random.uniform(gamma1min, 3, nwalkers)
	pos[:,3] = np.random.uniform(gamma2min, 3, nwalkers)
	pos[:,4] = np.random.uniform(gamma3min, 3, nwalkers)
	pos[:,5] = np.random.uniform(gamma4min, 3, nwalkers)""")
	else:
		f.write("""
	partrack = 5
	pos[:,0] = np.random.uniform(rhomin, rhomax, nwalkers)
	pos[:,1] = np.random.uniform(rsmin, rsmax, nwalkers)
	pos[:,2] = np.random.uniform(alphamin ,alphamax, nwalkers)
	pos[:,3] = np.random.uniform(betamin, betamax, nwalkers)
	pos[:,4] = np.random.uniform(gammamin , gammamax, nwalkers)""")
	if anisotropy == 'Const':

		f.write("""
	pos[:,partrack] = np.random.uniform(beta0min, beta0max, nwalkers)
	partrack = partrack + 1""")
	else:
		f.write("""
	pos[:,partrack] = np.random.uniform(beta0min, beta0max, nwalkers)
	pos[:,partrack + 1] = np.random.uniform(betainfmin, betainfmax, nwalkers)
	pos[:,partrack + 2] = np.random.uniform(ramin, ramax, nwalkers)
	pos[:,partrack + 3] = np.random.uniform(etamin, etamax, nwalkers)
	partrack = partrack + 4""")
	if plummer == 'Plummer3':
		f.write("""
	pos[:,partrack] = np.random.uniform(m1min, m1max, nwalkers)
	pos[:,partrack + 1]= np.random.uniform(a1min, a1max, nwalkers)
	pos[:,partrack + 2]= np.random.uniform(m2min, m2max, nwalkers)
	pos[:,partrack + 3]= np.random.uniform(a2min, a2max, nwalkers)
	pos[:,partrack + 4]= np.random.uniform(m3min, m3max, nwalkers)
	pos[:,partrack + 5]= np.random.uniform(a3min, a3max, nwalkers)

	partrack = partrack + 6""" + '\n' + '\n')


	



	f.write("\t" + "pos[:,-1] = np.random.uniform(mstarmin, mstarmax, nwalkers)" + "\n")
	

	if darkmatter == 'PL':
		f.write("""
	vec_lnprior = np.apply_along_axis(lnprior, 1, pos)
	junk, = np.where(np.isfinite(vec_lnprior) == False)
	not_junk = nwalkers - np.size(junk)
	

	while not_junk != nwalkers:
		

		pos[junk,0] = np.random.uniform(rho0min, rho0max, nwalkers - not_junk)
		pos[junk,1] = np.random.uniform(gamma0min, 3, nwalkers - not_junk)
		pos[junk,2] = np.random.uniform(gamma1min, 3, nwalkers - not_junk)
		pos[junk,3] = np.random.uniform(gamma2min, 3, nwalkers - not_junk)
		pos[junk,4] = np.random.uniform(gamma3min, 3, nwalkers - not_junk)
		pos[junk,5] = np.random.uniform(gamma4min, 3, nwalkers - not_junk)
		
		vec_lnprior = np.apply_along_axis(lnprior, 1, pos)
		junk, = np.where(np.isfinite(vec_lnprior) == False)
		not_junk = nwalkers - np.size(junk)
		
		
	f.write("\t" + "pos = pos[noconst_var]")

""")



	

	
	# log(ra/rc) = 0.1 ;  log(ra) - log(rc) = 0.1 ; 0.1 + log(rc)
		
	f.write("\t"+ """sampler = emcee.EnsembleSampler(nwalkers, ndim, logprob, pool = pool)

	try:
		pos = np.loadtxt(workdir  + project_name + '/%s' % galaxy_number +'/%s_LastWalkerPos' % galaxy_number + project_name + ".txt" )
		print 'Got pos'
		for result in sampler.sample(pos, 1):
			pos = result[0]
			prob = result[1]
			state = result[2]

		chains = np.genfromtxt(workdir + '/' + project_name + '/%s' % galaxy_number +'/%s_' %galaxy_number + "Chains" + project_name + ".txt")

		completed = int(float(len(chains))/float(nwalkers))
		tot_iter = steps -  completed


	except IOError:
		try:
			chains = np.genfromtxt(workdir + '/' + project_name + '/%s' % galaxy_number +'/%s_' %galaxy_number + "Chains" + project_name + ".txt")
			last_chains = chains[-nwalkers*100:]
			split_ch = np.array_split(last_chains, nwalkers)
			pos = np.zeros((nwalkers,ndim))
			for c in range(0, nwalkers):
				pos[c] = split_ch[c][-1,0:ndim]
			for result in sampler.sample(pos, 1):
				pos = result[0]
				prob = result[1]
				state = result[2]

			completed = int(float(len(chains))/float(nwalkers))
			tot_iter = steps -  completed



		except ValueError:
			f_chains = open(workdir + '/' + project_name + '/%s' % galaxy_number +'/%s_' %galaxy_number + "Chains" + project_name + ".txt", "w")
			f_chains.close()

			#pos = [free_param_list + 1e-2 * np.random.randn(ndim) for i in range(nwalkers)]

			pos, prob,state = sampler.run_mcmc(pos, 1)
		        tot_iter = steps
			
		 		

				
		except IOError:
			#pos = [free_param_list+ 1e-2 * np.random.randn(ndim) for i in range(nwalkers)]

			pos, prob,state = sampler.run_mcmc(pos, 1)
			tot_iter = steps
			
				

	sampler.reset()	

		
	
		
		

		


	for i in range(0, (tot_iter)/100):
		print 'Starting', i
		pos, prob, state = sampler.run_mcmc(pos, 100, lnprob0 = prob, rstate0 = state)
		print 'Finished', i 
		samples = sampler.flatchain
		pp = sampler.lnprobability

		if completed + i*100 > burn_in:
			out = np.zeros((len(samples), ndim + 1))
			out[:,0:ndim] = samples
			out[:,ndim] = pp.reshape((len(samples)))
			f_handle = file(workdir + '/' + project_name + '/%s' % galaxy_number +'/%s_Chains' %galaxy_number + project_name + ".txt", 'a+')
			np.savetxt(f_handle,out)
			f_handle.close()
		sampler.reset()
		
	
	


	np.savetxt(workdir + '/' + project_name + '/%s' % galaxy_number + '/%s_LastWalkerPos' % galaxy_number + project_name + ".txt" , pos)
		
		
	dt = time.time() - t_org
	print 'Finished running chains in %f seconds ' %(dt)

	

""")

t = time.time()
mpi_time = os.system("mpiexec -n %d" %(cores) +  " python " + workdir + "/" + project_name + "/Submissions/" + "script_bin_%s" %project_name + "_%s" %(galaxy_number) +  ".py")
#os.system("python " + "script_bin_%s" %volume + "_%d" %int(subhalo) + "_%d" %int(observer) +  ".py")
mpi_time = time.time() - t
print("MPI took {0:.1f} seconds".format(mpi_time))
#print("{0:.1f} times faster than serial".format(serial_time / mpi_time))





