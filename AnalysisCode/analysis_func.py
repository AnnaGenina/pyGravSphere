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
from scipy.optimize import minimize

def plummer(x,norm,rc):
        return np.log10(   (1./(np.pi * rc * rc))    / (  1 + x*x/(rc*rc)  )**2. ) + norm

def min_plummer(args,x_data, y_data,err_data):
        norm,rc = args
        return np.sum(  (  10**plummer(x_data,norm,rc) - y_data  )**2. / err_data**2. )


def plummer_proj_comp(x,mass,shape):

	return 10**mass/(np.pi * (10**shape)**2.) * (1. + x**2./(10**shape)**2)**(-2.)
def plummer_proj_sum(args, x_data, n_comp):

	sums_proj = 0

	start = 0
	for i in range(0, n_comp):
		sums_proj = sums_proj + plummer_proj_comp(x_data,args[start], args[start + 1])

		start = start + 2

	return np.log10(sums_proj)


def baes(x, b0, binf, ra, eta):
	return (b0 + binf * (x/ra)**eta) / ( 1 + (x/ra)**eta)

def zhao(x, rhos, rs,alpha,beta,gamma):
	return 10**rhos * (x/(10**rs))**(-gamma) * (1 + (x/(10**rs))**alpha)**(-(beta-gamma)/alpha)

def mass_zhao(x,rhos, rs,alpha,beta,gamma):
	return 4. * np.pi * x**2. * zhao(x,rhos, rs,alpha,beta,gamma)

def dens(r,rho0, bins, gammas,r_c):
    output = np.zeros_like(r)
    bins = bins*r_c
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


def mass_dens(x,rho0, bins, gammas,r_c):

    output = np.empty_like(x)

    bins = bins*r_c

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




def return_beta(chains, options, priors, min_r, max_r, points,codedir,workdir, project_name, galaxy):
	sys.path.append(codedir)
	from GStools import gsTools
	from GSpro import gal_input
	kindat,lightpower,vir_shape,surfden,r_c, stellar_mass = gal_input.galaxy_data_read(galaxy, workdir + '/GalaxyData/')

	
	r = np.logspace(np.log10(min_r), np.log10(max_r), points)

	beta_opt = options[1]
	dark_opt = options[0]
	
	fixed_params, = np.where(priors[:,6] == 'True')

	if dark_opt == 'PL':
		fixed, = np.where(priors[0:6, 6] == 'True')
		starter = 6 
		effective = starter - np.size(fixed)

	if dark_opt == 'Zhao':
		fixed, = np.where(priors[0:5, 6] == 'True')
		starter = 5 
		effective = starter - np.size(fixed)

	if beta_opt == 'Baes':
		
		beta_priors = priors[starter: starter + 4, :]
		fixed_beta_params, = np.where(beta_priors[:, 6] == 'True')
		free_beta_params, = np.where(beta_priors[:, 6] == 'False')
		samples = np.zeros((len(chains), 4))
		for p in range(0, np.size(fixed_beta_params)):
			if beta_priors[fixed_beta_params[p], 0] == 'ra':
				samples[:,fixed_beta_params[p]] = np.ones((len(samples))) * np.log10(float(beta_priors[fixed_beta_params[p], 5]) * r_c)
			else:
				samples[:,fixed_beta_params[p]] = np.ones((len(samples))) * float(beta_priors[fixed_beta_params[p], 5])
		for p in range(0, np.size(free_beta_params)):
			samples[:,free_beta_params[p]] = chains[:,effective]
			effective = effective + 1	


		betas = [[] for i in range(len(r))]

		samples[:,0] = 2.*samples[:,0] / (1. + samples[:,0]) #convert to actual values
		samples[:,1] = 2.*samples[:,1] / (1. + samples[:,1])
		samples[:,2] = 10**(samples[:,2])


		for sample in samples:

			b0,binf,ra,eta = sample

			bet_prof = baes(r, b0,binf,ra,eta)

			for b in range(0, len(r)):

				betas[b].append(float(bet_prof[b]))



		betas = np.concatenate([betas])
		betas = betas.T
		f_handle = file(workdir + project_name + '/Analysis/Output/' + '%s_Beta'%galaxy + '.txt' , 'w')
		np.savetxt(f_handle,betas, delimiter = '\t')
	       	f_handle.close()

		res = gsTools.get_lims(betas,r)
		np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Beta'% galaxy +  'Lims.txt', res)


	elif beta_opt == 'Const':
		beta_priors = priors[starter: starter + 1, :]
		fixed_beta_params, = np.where(beta_priors[:, 6] == 'True')
		free_beta_params, = np.where(beta_priors[:, 6] == 'False')
		samples = np.zeros((len(chains), 1))
		for p in range(0, np.size(fixed_beta_params)):
			samples[:,fixed_beta_params[p]] = np.ones((len(samples))) * float(beta_priors[fixed_beta_params[p], 5])
		for p in range(0, np.size(free_beta_params)):
			samples[:,free_beta_params[p]] = chains[:,effective]
			effective = effective + 1			

		betas = [[] for i in range(len(r))]

		samples[:,0] = 2.*samples[:,0] / (1. + samples[:,0]) #convert to actual values
		

		for sample in samples:


			bet_prof = np.ones((len(r))) * sample

			for b in range(0, len(r)):

				betas[b].append(float(bet_prof[b]))

		betas = np.concatenate([betas])
		betas = betas.T
		f_handle = file(workdir + project_name + '/Analysis/Output/' + '%s_Beta'%galaxy + '.txt' , 'w')
		np.savetxt(f_handle,betas, delimiter = '\t')
	       	f_handle.close()

		res = gsTools.get_lims(betas,r)
		np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Beta'% galaxy +  'Lims.txt', res)


def return_mass(chains, options, priors, min_r, max_r, points,codedir,workdir,project_name, galaxy):
	sys.path.append(codedir)
	from GStools import gsTools
	from GSpro import gal_input
	kindat,lightpower,vir_shape,surfden,r_c, stellar_mass = gal_input.galaxy_data_read(galaxy, workdir + '/GalaxyData/')

	r = np.logspace(np.log10(min_r), np.log10(max_r), points)

	dark_opt = options[0]
	
	effective = 0

	if dark_opt == 'PL':
		
		dark_priors = priors[0:6, :]
		fixed_dm_params, = np.where(dark_priors[:, 6] == 'True')
		free_dm_params, = np.where(dark_priors[:, 6] == 'False')
		samples = np.zeros((len(chains), 6))
		for p in range(0, np.size(fixed_dm_params)):
			samples[:,fixed_dm_params[p]] = np.ones((len(samples))) * float(dark_priors[fixed_dm_params[p], 5])
		for p in range(0, np.size(free_dm_params)):
			samples[:,free_dm_params[p]] = chains[:,effective]
			effective = effective + 1	

		bins = np.array([0.25,0.5,1.0,2.0,4.0])

		bin_edges = np.array([0.125, 0.25, 0.50, 1, 2 , 4 ])*r_c
		mid_bin = (np.log10(bin_edges[1:]) + np.log10(bin_edges[:-1]))/2.
		mid_bin = 10**mid_bin

		tot_bins = np.array([bin_edges[0], mid_bin[0], bin_edges[1], bin_edges[1], mid_bin[1], bin_edges[2],bin_edges[2], mid_bin[2], bin_edges[3], bin_edges[3], mid_bin[3], bin_edges[4], bin_edges[4], mid_bin[4], bin_edges[5]])


		js = [[] for i in range(len(r))]
		masses =   [[] for i in range(len(r))]
		grad_rho =   [[] for i in range(len(tot_bins))]


		for sample in samples:

			rho0,gamma0,gamma1,gamma2,gamma3,gamma4= sample
		    	gammas = [gamma0,gamma1,gamma2,gamma3,gamma4]
		    	big_j = dens(r, rho0, bins, gammas,r_c)
			big_j = np.log10(big_j)

			for b in range(0, len(big_j)):

				mass_try = mass_dens(np.array([r[b]]),rho0, bins, gammas,r_c)
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

		f = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_Density.txt', 'w')
		np.savetxt(f, js, delimiter = '\t')
		f.close()
		f2 = open(workdir + project_name + '/Analysis/Output/' + '%s' % galaxy +  '_Mass.txt' , 'w')
		np.savetxt(f2,masses, delimiter = '\t')
		f2.close()
		f3 = open(workdir + project_name + '/Analysis/Output/' + '%s' % galaxy +  '_LogLog.txt' , 'w')
		np.savetxt(f3,grad_rho, delimiter = '\t')
		f3.close()

		res = gsTools.get_lims_loglog(grad_rho, tot_bins)
		np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_LogLog'% galaxy + 'Lims.txt', res)
		
		res = gsTools.get_lims(js,r)
		np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_Density'% galaxy +  'Lims.txt', res)

		res = gsTools.get_lims(masses,r)
		np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_Mass'% galaxy +  'Lims.txt', res)
		


	if dark_opt == 'Zhao':
		huge_scale = np.logspace(np.log10(min_r)-0.01,np.log10(max_r)+0.01, 10000)
		inter_zhao = (huge_scale[1:] + huge_scale[:-1])/2.
		diff_zhao_dist = huge_scale[1:] - huge_scale[:-1]


		dark_priors = priors[0:5,:]
		fixed_dm_params, = np.where(dark_priors[:, 6] == 'True')
		free_dm_params, = np.where(dark_priors[:, 6] == 'False')
		samples = np.zeros((len(chains), 5))
		for p in range(0, np.size(fixed_dm_params)):
			samples[:,fixed_dm_params[p]] = np.ones((len(samples))) * float(dark_priors[fixed_dm_params[p], 5])
		for p in range(0, np.size(free_dm_params)):
			samples[:,free_dm_params[p]] = chains[:,effective]
			effective = effective + 1

		js = [[] for i in range(len(r))]
		masses =   [[] for i in range(len(r))]
		grad_rho =   [[] for i in range(len(r))]
	

		for sample in samples:

			rhos,rs,alpha,beta,gamma = sample
			big_j = zhao(r, rhos,rs,alpha,beta,gamma)

			zhao_huge =  zhao(huge_scale, rhos,rs,alpha,beta,gamma)
			zhao_inter = zhao(inter_zhao, rhos,rs,alpha,beta,gamma)
			grad_zhao = (zhao_huge[1:] - zhao_huge[:-1])/diff_zhao_dist

			loglog = inter_zhao/zhao_inter * grad_zhao

			grad_prof = interp1d(inter_zhao, loglog, kind = 'linear')

			spec_loglog = grad_prof(r)



			big_j = np.log10(big_j)

			
			

			for b in range(0, len(big_j)):

				lin = np.logspace(np.minimum(np.log10(min_r)-1, -3),np.log10(r[b]),1000)
				mass_try = simps(np.log(10) * lin * mass_zhao(lin,rhos,rs,alpha,beta,gamma), np.log10(lin), even = 'avg')
				js[b].append(float(big_j[b]))
				masses[b].append(float(np.log10(mass_try)))
				grad_rho[b].append(float(spec_loglog[b]))

		js = np.concatenate([js])
		js = js.T
		masses = np.concatenate([masses])
		masses = masses.T
		grad_rho = np.concatenate([grad_rho])
		grad_rho = grad_rho.T



		f = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_Density.txt' , 'w')
		np.savetxt(f, js, delimiter = '\t')
		f.close()
		f2 = open(workdir + project_name + '/Analysis/Output/' + '%s' % galaxy +  '_Mass.txt' , 'w')
		np.savetxt(f2,masses, delimiter = '\t')
		f2.close()
		f3 = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_LogLog.txt' , 'w')
		np.savetxt(f3,grad_rho, delimiter = '\t')
		f3.close()

		res = gsTools.get_lims(grad_rho,r)
		np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_LogLog'% galaxy + 'Lims.txt', res)

		res = gsTools.get_lims(js,r)
		np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_Density'% galaxy +  'Lims.txt', res)

		res = gsTools.get_lims(masses,r)
		np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_Mass'% galaxy +  'Lims.txt', res)


def return_plummer(chains, options, priors, min_r, max_r, points,codedir,workdir, project_name, galaxy):
	sys.path.append(codedir)
	from GStools import gsTools
	from GSpro import gal_input
	kindat,lightpower,vir_shape,surfden,r_c, stellar_mass = gal_input.galaxy_data_read(galaxy, workdir + '/GalaxyData/')


	r = np.logspace(np.log10(min_r), np.log10(max_r), points)

	beta_opt = options[1]
	dark_opt = options[0]
	plummer_opt = options[2]
	
	fixed_params, = np.where(priors[:,6] == 'True')

	effective = 0

	if dark_opt == 'PL':
		fixed, = np.where(priors[0:6, 6] == 'True')
		effective = 6 - np.size(fixed)
		starter = 6
	elif dark_opt == 'Zhao':
		fixed, = np.where(priors[0:5, 6] == 'True')
		effective = 5 - np.size(fixed)	
		starter = 5
	if beta_opt == 'Baes':
		fixed, = np.where(priors[starter:starter+4, 6] == 'True')
		starter = starter + 4
		effective = effective + 4 - np.size(fixed)
	elif beta_opt == 'Const':	
		fixed, = np.where(priors[starter:starter+1, 6] == 'True')
		starter = starter + 1
		effective = effective + 1 - np.size(fixed)
		
	if plummer_opt == 'Plummer3':
		plummer_priors = priors[starter:starter + 6,:]
		fixed_plummer_params, = np.where(plummer_priors[:, 6] == 'True')
		free_plummer_params, = np.where(plummer_priors[:, 6] == 'False')
		samples = np.zeros((len(chains), 6))
		for p in range(0, np.size(fixed_plummer_params)):
			if plummer_priors[fixed_plummer_params[p], 0] == 'm1':
				samples[:,fixed_plummer_params[p]] = np.ones((len(samples))) * np.log10(lightpower[0] * float(plummer_priors[fixed_plummer_params[p], 5]))
			elif plummer_priors[fixed_plummer_params[p], 0] == 'm2':
				samples[:,fixed_plummer_params[p]] = np.ones((len(samples))) * np.log10(lightpower[1] * float(plummer_priors[fixed_plummer_params[p], 5]))
			elif plummer_priors[fixed_plummer_params[p], 0] == 'm3':
				samples[:,fixed_plummer_params[p]] = np.ones((len(samples))) * np.log10(lightpower[2] * float(plummer_priors[fixed_plummer_params[p], 5]))
			elif plummer_priors[fixed_plummer_params[p], 0] == 'a1':
				samples[:,fixed_plummer_params[p]] = np.ones((len(samples))) * (lightpower[3] * float(plummer_priors[fixed_plummer_params[p], 5]))
			elif plummer_priors[fixed_plummer_params[p], 0] == 'a2':
				samples[:,fixed_plummer_params[p]] = np.ones((len(samples))) * (lightpower[4] * float(plummer_priors[fixed_plummer_params[p], 5]))
			elif plummer_priors[fixed_plummer_params[p], 0] == 'a3':
				samples[:,fixed_plummer_params[p]] = np.ones((len(samples))) * (lightpower[5] * float(plummer_priors[fixed_plummer_params[p], 5]))


			
		for p in range(0, np.size(free_plummer_params)):
			samples[:,free_plummer_params[p]] = chains[:,effective]
			effective = effective + 1
	
		plum_l = [[] for i in range(len(r))]		

		for sample in samples:

			m1,a1,m2,a2,m3,a3 = sample
			a1 = np.log10(a1)
			a2 = np.log10(a2)
			a3 = np.log10(a3)



			for b in range(0, len(r)):
				plum_l[b].append(plummer_proj_sum([m1,a1,m2,a2,m3,a3], r[b], 3))




		plum_l = np.concatenate([plum_l])
		plum_l = plum_l.T

		f = open(workdir + project_name + '/Analysis/Output/' + '%s' % galaxy +  '_Plummer.txt' , 'w')
		np.savetxt(f, plum_l, delimiter = '\t')
		f.close()

		res = gsTools.get_lims(plum_l,r)
		np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Plummer'% galaxy + 'Lims.txt', res)


		
	else:
		
		print 'No limits for single Plummer profile'

def return_sigma_vsp(chains, options, priors, min_r, max_r, points,codedir,workdir,project_name, galaxy):

	sys.path.append(codedir)
	from GSpro import gal_input
	import gravsphere
	from GStools import gsTools

	kindat,lightpower,vir_shape,surfden,r_c, stellar_mass = gal_input.galaxy_data_read(galaxy, workdir + '/GalaxyData/')


	r = np.logspace(np.log10(min_r), np.log10(max_r), points)

	beta_opt = options[1]
	dark_opt = options[0]
	plummer_opt = options[2]

	nparams = 0
	if dark_opt == 'PL':
		nparams = nparams + 6
	elif dark_opt == 'Zhao':
		nparams = nparams + 5
	if beta_opt == 'Baes':
		nparams = nparams + 4
	elif beta_opt == 'Const':
		nparams = nparams + 1
	if plummer_opt == 'Plummer3':
		nparams = nparams + 6
	else:
		nparams = nparams
	nparams = nparams + 1

	if dark_opt == 'PL':
		effective = 0
		fixed_params, = np.where(priors[:,6] == 'True')
		free_params, = np.where(priors[:,6] == 'False')
		samples = np.zeros((len(chains), nparams))
		for p in range(0, np.size(fixed_params)):

			if priors[fixed_params[p], 0] == 'ra':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*r_c)
			elif priors[fixed_params[p], 0] == 'm1':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*lightpower[0])
			elif priors[fixed_params[p], 0] == 'm2':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*lightpower[1])
			elif priors[fixed_params[p], 0] == 'm3':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*lightpower[2])
			elif priors[fixed_params[p], 0] == 'a1':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * (float(priors[fixed_params[p], 5])*lightpower[3])
			elif priors[fixed_params[p], 0] == 'a2':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * (float(priors[fixed_params[p], 5])*lightpower[4])
			elif priors[fixed_params[p], 0] == 'a3':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * (float(priors[fixed_params[p], 5])*lightpower[5])
			elif priors[fixed_params[p], 0] == 'mstar':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*stellar_mass)
			else:
				samples[:,fixed_params[p]] = np.ones((len(samples))) * float(priors[fixed_params[p], 5])
			


		for p in range(0, np.size(free_params)):
			samples[:,free_params[p]] = chains[:,effective]
			effective = effective + 1

		


		min_rad = np.log10(r_c/100)
		max_rad = np.log10(r_c*50)

		vs1arr = []
		vs2arr = []
		res_arr = []


		n_params = len(samples[0])

		chlen = len(samples)

		plum_fit =minimize(min_plummer, [np.log10(np.amax(surfden[:,1])),0.5], args=(surfden[:,0], surfden[:,1],surfden[:,2],), bounds = ((0, 10),(0.1,5),))
		m1,a1 = plum_fit.x
	
		for sample in samples:

			if n_params == 17:
				rho0, gamma0, gamma1, gamma2, gamma3, gamma4, beta0, betainf, ra, eta, m1,a1,m2,a2,m3,a3,mstar = sample
			if n_params == 14:
				rho0, gamma0, gamma1, gamma2, gamma3, gamma4, beta0, m1,a1,m2,a2,m3,a3,mstar = sample
				betainf = beta0
				ra = 0
				eta = 1
			if n_params == 8:
				rho0, gamma0, gamma1, gamma2, gamma3, gamma4, beta0,mstar = sample
				betainf = beta0
				ra = 0
				eta = 1
				m2 = -5
				a2 = 1
				m3 = -5
				a3 = 1
			if n_params == 11:
				rho0, gamma0, gamma1, gamma2, gamma3, gamma4, beta0,betainf, ra, eta,mstar = sample
				m2 = -5
				a2 = 1
				m3 = -5
				a3 = 1

			

			rho_params = np.array([r_c, rho0, gamma0,gamma1,gamma2,gamma3,gamma4])
			beta_params = np.array([ 2*beta0/(1+beta0), 2*betainf/(1+betainf)  ,10**ra, eta])
			plum_params = np.array([10**m1,a1,10**m2,a2,10**m3,a3])


			results = np.zeros_like(kindat[:,0])
			vsp1 = np.zeros((1), dtype = np.float64)
			vsp2 = np.zeros((1), dtype = np.float64)

			gravsphere.PowerLawFitfunc(kindat[:,0], kindat[:,1], kindat[:,2], rho_params, beta_params, plum_params, vir_shape,10**mstar, r_c, 100, results, vsp1, vsp2, min_rad, max_rad)

			vs1arr.append(float(vsp1))
			vs2arr.append(float(vsp2))
			res_arr.append(results)


		vs1arr = np.array(vs1arr)
		vs2arr = np.array(vs2arr)
		res_arr = np.concatenate([res_arr])


		f1 = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_SigLos.txt', 'w')
		f2 = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_VSP1.txt', 'w')
		f3 = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_VSP2.txt', 'w')

		np.savetxt(f2, vs1arr, delimiter = '\t')
		np.savetxt(f3, vs2arr, delimiter = '\t')
		np.savetxt(f1, res_arr, delimiter = '\t')

		res = gsTools.get_lims_sig(res_arr, workdir, galaxy)
		np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_SigLos'% galaxy + 'Lims.txt', res)


	elif dark_opt == 'Zhao':
		effective = 0
		fixed_params, = np.where(priors[:,6] == 'True')
		free_params, = np.where(priors[:,6] == 'False')
		samples = np.zeros((len(chains), nparams))
		for p in range(0, np.size(fixed_params)):
			if priors[fixed_params[p], 0] == 'ra':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*r_c)
			elif priors[fixed_params[p], 0] == 'm1':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*lightpower[0])
			elif priors[fixed_params[p], 0] == 'm2':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*lightpower[1])
			elif priors[fixed_params[p], 0] == 'm3':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*lightpower[2])
			elif priors[fixed_params[p], 0] == 'a1':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * (float(priors[fixed_params[p], 5])*lightpower[3])
			elif priors[fixed_params[p], 0] == 'a2':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * (float(priors[fixed_params[p], 5])*lightpower[4])
			elif priors[fixed_params[p], 0] == 'a3':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * (float(priors[fixed_params[p], 5])*lightpower[5])
			elif priors[fixed_params[p], 0] == 'mstar':
				samples[:,fixed_params[p]] = np.ones((len(samples))) * np.log10(float(priors[fixed_params[p], 5])*stellar_mass)

			else:
				samples[:,fixed_params[p]] = np.ones((len(samples))) * float(priors[fixed_params[p], 5])
		for p in range(0, np.size(free_params)):
			samples[:,free_params[p]] = chains[:,effective]
			effective = effective + 1


		kindat,lightpower,vir_shape,surfden,r_c, stellar_mass = gal_input.galaxy_data_read(galaxy, workdir + '/GalaxyData/')

		plum_fit =minimize(min_plummer, [np.log10(np.amax(surfden[:,1])),0.5], args=(surfden[:,0], surfden[:,1],surfden[:,2],), bounds = ((0, 10),(0.1,5),))
		m1,a1 = plum_fit.x		

		min_rad = np.log10(r_c/100)
		max_rad = np.log10(r_c*50)

		vs1arr = []
		vs2arr = []
		res_arr = []

		n_params = len(samples[0])

		chlen = len(samples)

    	

		for sample in samples:
			if n_params == 16:
				rho0, rs, alpha, beta, gamma , beta0, betainf, ra, eta, m1,a1,m2,a2,m3,a3,mstar = sample
			if n_params == 13:
				rho0, rs, alpha, beta, gamma , beta0, m1,a1,m2,a2,m3,a3,mstar = sample
				betainf = beta0
				ra = 0
				eta = 1
			if n_params == 7:
				rho0, rs, alpha, beta, gamma , beta0,mstar = sample
				betainf = beta0
				ra = 0
				eta = 1
				m2 = -5
				a2 = 1
				m3 = -5
				a3 = 1
			if n_params == 10:
				rho0, rs, alpha, beta, gamma , beta0,betainf, ra, eta,mstar = sample
				m2 = -5
				a2 = 1
				m3 = -5
				a3 = 1


			


			rho_params = np.array([r_c, rho0, rs, alpha, beta, gamma ])
			beta_params = np.array([ 2*beta0/(1+beta0), 2*betainf/(1+betainf)  ,10**ra, eta])
			plum_params = np.array([10**m1,a1,10**m2,a2,10**m3,a3])


			results = np.zeros_like(kindat[:,0])
			vsp1 = np.zeros((1), dtype = np.float64)
			vsp2 = np.zeros((1), dtype = np.float64)

			gravsphere.ZhaoFitfunc(kindat[:,0], kindat[:,1], kindat[:,2], rho_params, beta_params, plum_params, vir_shape,10**mstar, r_c, 100, results, vsp1, vsp2, min_rad, max_rad)

			vs1arr.append(float(vsp1))
			vs2arr.append(float(vsp2))
			res_arr.append(results)


		vs1arr = np.array(vs1arr)
		vs2arr = np.array(vs2arr)
		res_arr = np.concatenate([res_arr])


		f1 = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_SigLos.txt', 'w')
		f2 = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_VSP1.txt', 'w')
		f3 = open(workdir + project_name + '/Analysis/Output/'  + '%s' % galaxy +  '_VSP2.txt', 'w')

		np.savetxt(f2, vs1arr, delimiter = '\t')
		np.savetxt(f3, vs2arr, delimiter = '\t')
		np.savetxt(f1, res_arr, delimiter = '\t')

		res = gsTools.get_lims_sig(res_arr, workdir, galaxy)
		np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_SigLos'% galaxy + 'Lims.txt', res)


