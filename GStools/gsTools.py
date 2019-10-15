import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def checkdirs(workdir, codedir):
	kinphot = os.path.isdir(workdir + "/KinPhotDat")
	if kinphot == False:
		os.system("mkdir " + workdir + "/KinPhotDat")	
	galdat = os.path.isdir(workdir + "/GalaxyData")
	if galdat == False:
		os.system("mkdir " + workdir + "/GalaxyData")
	gal_list = os.path.exists(workdir + "/galaxy_list.txt")
	if gal_list == False:
		os.system("touch " + workdir + "/galaxy_list.txt")
	projects =  os.path.exists(workdir + "/projects.txt")
	if projects == False:
		os.system("touch " + workdir + "/projects.txt")	
	code =  os.path.exists(workdir + "/pygravsphere.py")
	if code == False:
		os.system("cp " + codedir + "/pygravsphere.py " + workdir + " /pygravsphere.py")

def get_cmap(n, name='hsv'):
   
	return plt.cm.get_cmap(name, n)

def plot_chains(samples,  param_names, nparams, workdir, project_name, galaxy):
	nwalkers = 1000
	ll = len(samples)/(100*nwalkers)
	sample_split = np.array_split(samples, ll)			
	sample_new = [[] for i in range(nwalkers)]

	for s in range(0, ll):
		second_split = np.array_split(sample_split[s], nwalkers)
		for s2 in range(nwalkers):				
			sample_new[s2].append(second_split[s2]) 

	sample_new = np.concatenate(np.concatenate(sample_new))
	
	for p in range(0, len(param_names)):

		fig, ax1 = plt.subplots(1,1, figsize = (16,4))
	
		cmap = get_cmap(nwalkers)

		plot_sample = np.array_split(sample_new, nwalkers)
		for i in range(nwalkers):
			if i%25 == 0:
				ax1.plot(plot_sample[i][:,p],  color =cmap(i), lw = 1.0, alpha = 0.5)
		ax1.set_ylabel(param_names[p], fontsize = 12)	
		ax1.set_xlabel('Step', fontsize = 12)	

		
		fig.savefig(workdir + project_name + '/%s/' % galaxy + 'Param_%d.png' % p, bbox_inches = "tight")
		plt.close()


def plot_triangle(samples,  param_names, nparams, workdir, project_name, galaxy):


	like = samples[:,-1]
	fin, = np.where(np.isfinite(like) == True)	
	samples = samples[fin]

	if len(param_names) == 10:
		
		fig = corner.corner(samples[:,0:10], labels=param_names, quantiles=[0.16, 0.5, 0.84],title_kwargs={"fontsize": 14})
		fig.savefig(workdir + project_name + '/%s/' % galaxy + 'Triangle.png')
		plt.close()
	
		param_names.extend([r"$m_1$", r"$a_1$", r"$m_2$", r"$a_2$",r"$m_3$", r"$a_3$"])
		fig = corner.corner(samples[:,0:-2], labels=param_names, quantiles=[0.16, 0.5, 0.84],title_kwargs={"fontsize": 14})
		fig.savefig(workdir + project_name + '/%s/' % galaxy + 'TriangleFull.png')
		plt.close()

	elif len(param_names) == 9:

		fig = corner.corner(samples[:,0:9], labels=param_names, quantiles=[0.16, 0.5, 0.84],title_kwargs={"fontsize": 14})
		fig.savefig(workdir + project_name + '/%s/' % galaxy + 'Triangle.png')
		plt.close()

def get_lims(data,r):


	
	
	out = np.zeros((len(r), 6))
	

	out[:,0] = r
	out[:,1] = np.median(data, axis = 0)
	out[:,2] = np.percentile(data, 16, axis = 0)
	out[:,3]= np.percentile(data, 84, axis = 0)
	out[:,4]= np.percentile(data, 2.5, axis = 0)
	out[:,5]= np.percentile(data, 97.5, axis = 0)

	return out


def get_lims_all(data,r):


	
	
	out = np.zeros((len(r), 102))
	out[:,0] = r
	for i in range(0, 101):
		
		out[:,i+1] = np.percentile(data, i, axis = 0)

	

		
	return out


def get_lims_loglog(data, tot_bins):


	
	out = np.zeros((len(tot_bins), 6))

	out[:,0] = tot_bins
	out[:,1] = np.median(data, axis = 0)
	out[:,2] = np.percentile(data, 16, axis = 0)
	out[:,3]= np.percentile(data, 84, axis = 0)
	out[:,4]= np.percentile(data, 2.5, axis = 0)
	out[:,5]= np.percentile(data, 97.5, axis = 0)

	return out


def get_lims_sig(data, workir, galaxy):
	pos = np.loadtxt(workdir + '/GalaxyData/%s_KinDat'% int(galaxy) + '.txt')
	pos = pos[:,0]

	out = np.zeros((len(pos), 6))

	out[:,0] = pos
	out[:,1] = np.median(data, axis = 0)
	out[:,2] = np.percentile(data, 16, axis = 0)
	out[:,3]= np.percentile(data, 84, axis = 0)
	out[:,4]= np.percentile(data, 2.5, axis = 0)
	out[:,5]= np.percentile(data, 97.5, axis = 0)

	return out

def create_sub(project_name, num_cores, timevar, workdir,codedir, anis, darkmatter, vsps, plummer, num_walkers, burn_in, steps, int_points, mpi_opt):
	if mpi_opt == 'y':
		core = int(num_cores)
		galaxies = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		
		sub_script = np.loadtxt("sub_script.txt", dtype = 'str')

		prestr1 = workdir + project_name + '/Submissions/'
		prestr2 = workdir + project_name + '/OutErr/'

		
		with open("sub_script.txt") as forg:
			newText=f.read()
			newText.replace('CORENUM', "%d" %core)
			newText.replace('-J GALID', '-J %s_' % galaxy  + '%s' % project_name)
			newText.replace('GALID.err', prestr2 + "%s_" %galaxy + "%s.err" % project_name)
			newText.replace('GALID.out', prestr2 + "%s_" %galaxy + "%s.out" % project_name )
			newText.replace('TIME', "%02d" %timevar + ":00:00")
			
		for galaxy in galaxies:
			with open(prestr1 + '%s.sh' %galaxy, "w") as f:
				f.write(newText + '\n')
				f.write("python " + codedir + '/write_script_mpi.py' + " %s"%(workdir) + " %s"%(codedir) + " %s" %(project_name) +  " %s"%(galaxy)  + " %s"%(core) + " %s"%(num_walkers) + " %s"%(burn_in) + " %s"%(steps)  + " %s"%(int_points) + " %s"%(darkmatter) + " %s"%(anis) +  " %s"%(vsps) + " %s"%(plummer)  )
				f.close()


		forg.close()
	else:
		galaxies = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')

		prestr1 = workdir + project_name + '/Submissions/'
		prestr2 = workdir + project_name + '/OutErr/'
		for galaxy in galaxies:
			with open(prestr1 + '%s.sh' %galaxy, "w") as f:
				f.write("#!/bin/sh" + '\n')
				f.write("python " + codedir + '/write_script.py' + " %s"%(workdir) + " %s"%(codedir) + " %s" %(project_name) +  " %s"%(galaxy)  + " %s"%(num_walkers) + " %s"%(burn_in) + " %s"%(steps)  + " %s"%(int_points) + " %s"%(darkmatter) + " %s"%(anis) +  " %s"%(vsps) + " %s"%(plummer)  )
				f.close()
				os.system("chmod u+x " + prestr1 + '%s.sh' % galaxy)
	
	


def create_ana_sub(project_name, workdir, codedir, dm_option, beta_option, plummer_option, samples, mpi_opt):

	if mpi_opt == 'y':
		core = 1
		galaxies = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')

		sub_script = np.loadtxt("sub_script.txt", dtype = 'str')

		prestr1 = workdir + project_name + '/Analysis/Submissions/'
		prestr2 = workdir + project_name + '/Analysis/OutErr/'

		
		with open("sub_script.txt") as forg:
			newText=f.read()
			newText.replace('CORENUM', "%d" %core)
			newText.replace('-J GALID', '-J %s_' % galaxy  + '%s' % project_name)
			newText.replace('GALID.err', prestr2 + "%s_" %galaxy + "%s.err" % project_name)
			newText.replace('GALID.out', prestr2 + "%s_" %galaxy + "%s.out" % project_name )
			newText.replace('TIME', "%02d" %timevar + ":00:00")
			
		for galaxy in galaxies:
			with open(prestr1 + '%s.sh' %galaxy, "w") as f:
				f.write(newText + '\n')


			
				if (dm_option == 'PL') and (beta_option == 'Baes') and (plummer_option == 'Plummer3'):
					starter = 10
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_plummer.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) + ' ' + str(int(starter)) + ' ' + str(int(samples)) + '\n')
				if (dm_option == 'Zhao') and (beta_option == 'Baes') and (plummer_option == 'Plummer3'):
					starter = 9
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_plummer.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) + ' ' + str(int(starter)) + ' ' + str(int(samples)) + '\n')
					

				if (dm_option == 'Zhao') and (beta_option == 'Baes'):
					f.write("python " + codedir + '/AnalysisCode/' + 'dens_chains_read_baes2.py ' + workdir +' ' +  project_name  + ' '+ str(galaxy) + ' ' + str(int(samples)) + '\n')
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_chains_read_baes.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy)  + ' ' + str(int(samples)) + '\n')
				elif (dm_option == 'PL') and (beta_option == 'Baes'):
					f.write("python " + codedir +   '/AnalysisCode/' + 'dens_chains_read_baes2PL.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) + ' ' + str(int(samples)) + '\n')
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_chains_read_baesPL.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy)  + ' ' + str(int(samples)) + '\n')
				elif (dm_option == 'Zhao') and (beta_option == 'Const'):
					f.write("python " + codedir + '/AnalysisCode/' + 'dens_chains_read.py ' + workdir +' ' +  project_name  + ' '+ str(galaxy) + ' 5 ' + ' ' + str(int(samples)) +  '\n')
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_chains_read_baes.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy)  + ' ' + str(int(samples)) + '\n')
				elif (dm_option == 'PL') and (beta_option == 'Const'):
					
					f.write("python " + codedir +   '/AnalysisCode/' + 'dens_chains_read.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) +' 6 ' + ' ' + str(int(samples)) +  '\n')
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_chains_read_baesPL.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy)  + ' ' + str(int(samples)) + '\n')

				if dm_option == 'PL':
					f.write("python " + codedir +  '/AnalysisCode/' + 'sigvs_pl.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) +' ' + str(int(samples)) +  '\n')
				else:
					f.write("python " + codedir +  '/AnalysisCode/' + 'sigvs_zhao.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) + ' ' + str(int(samples)) +  '\n')			
				f.close()

	else:

		core = 1
		galaxies = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')

		sub_script = np.loadtxt("sub_script.txt", dtype = 'str')

		prestr1 = workdir + project_name + '/Analysis/Submissions/'
		prestr2 = workdir + project_name + '/Analysis/OutErr/'

			
		for galaxy in galaxies:
			with open(prestr1 + '%s.sh' %galaxy, "w") as f:
				
			
				if (dm_option == 'PL') and (beta_option == 'Baes') and (plummer_option == 'Plummer3'):
					starter = 10
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_plummer.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) + ' ' + str(int(starter)) + ' ' + str(int(samples)) + '\n')
				if (dm_option == 'Zhao') and (beta_option == 'Baes') and (plummer_option == 'Plummer3'):
					starter = 9
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_plummer.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) + ' ' + str(int(starter)) + ' ' + str(int(samples)) + '\n')
					

				if (dm_option == 'Zhao') and (beta_option == 'Baes'):
					f.write("python " + codedir + '/AnalysisCode/' + 'dens_chains_read_baes2.py ' + workdir +' ' +  project_name  + ' '+ str(galaxy) + ' ' + str(int(samples)) + '\n')
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_chains_read_baes.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy)  + ' ' + str(int(samples)) + '\n')
				elif (dm_option == 'PL') and (beta_option == 'Baes'):
					f.write("python " + codedir +   '/AnalysisCode/' + 'dens_chains_read_baes2PL.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) + ' ' + str(int(samples)) + '\n')
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_chains_read_baesPL.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy)  + ' ' + str(int(samples)) + '\n')
				elif (dm_option == 'Zhao') and (beta_option == 'Const'):
					f.write("python " + codedir + '/AnalysisCode/' + 'dens_chains_read.py ' + workdir +' ' +  project_name  + ' '+ str(galaxy) + ' 5 ' + ' ' + str(int(samples)) +  '\n')
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_chains_read_baes.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy)  + ' ' + str(int(samples)) + '\n')
				elif (dm_option == 'PL') and (beta_option == 'Const'):
					
					f.write("python " + codedir +   '/AnalysisCode/' + 'dens_chains_read.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) +' 6 ' + ' ' + str(int(samples)) +  '\n')
					f.write("python " + codedir +  '/AnalysisCode/' + 'dens_chains_read_baesPL.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy)  + ' ' + str(int(samples)) + '\n')

				if dm_option == 'PL':
					f.write("python " + codedir +  '/AnalysisCode/' + 'sigvs_pl.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) +' ' + str(int(samples)) +  '\n')
				else:
					f.write("python " + codedir +  '/AnalysisCode/' + 'sigvs_zhao.py ' + workdir +' ' +  project_name  + ' ' + str(galaxy) + ' ' + str(int(samples)) +  '\n')			
				f.close()

				os.system("chmod u+x " + prestr1 + '%s.sh' %galaxy)


def submit_jobs(workdir, project_name, all_gals, sub_com):
	for galaxy in all_gals:
		if sub_com == False:		
			os.system('./' + workdir + project_name + '/Submissions/' + '%s.sh' % galaxy)
		else:
			os.system(sub_com + ' ' + workdir + project_name + '/Submissions/' + '%s.sh' % galaxy)


def preprocess(workdir, codedir, all_gals):
	sys.path.append(codedir)
	from GSpro import gal_input
	for galaxy in all_gals:
		gal_input.galaxy_data(galaxy, workdir + '/GalaxyData/', workdir + '/KinPhotDat/')

	return 0
