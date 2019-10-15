import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import corner
import glob




def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
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

		
		fig.savefig(workdir + project_name + '/Galaxy_%d/' % galaxy + 'Param_%d.png' % p, bbox_inches = "tight")
		plt.close()


def plot_triangle(samples,  param_names, nparams, workdir, project_name, galaxy):

	#index, = np.where(samples[:,-1] < np.amin(samples[:,-1])*10.0)
	#samples = samples[index]

	like = samples[:,-1]
	fin, = np.where(np.isfinite(like) == True)	
	samples = samples[fin]

	if len(param_names) == 10:
		
		fig = corner.corner(samples[:,0:10], labels=param_names, quantiles=[0.16, 0.5, 0.84],title_kwargs={"fontsize": 14})
		fig.savefig(workdir + project_name + '/Galaxy_%d/' % galaxy + 'Triangle.png')
		plt.close()
	
		param_names.extend([r"$m_1$", r"$a_1$", r"$m_2$", r"$a_2$",r"$m_3$", r"$a_3$"])
		fig = corner.corner(samples[:,0:-2], labels=param_names, quantiles=[0.16, 0.5, 0.84],title_kwargs={"fontsize": 14})
		fig.savefig(workdir + project_name + '/Galaxy_%d/' % galaxy + 'TriangleFull.png')
		plt.close()

	elif len(param_names) == 9:

		fig = corner.corner(samples[:,0:9], labels=param_names, quantiles=[0.16, 0.5, 0.84],title_kwargs={"fontsize": 14})
		fig.savefig(workdir + project_name + '/Galaxy_%d/' % galaxy + 'Triangle.png')
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
	pos = np.loadtxt(workdir + '/GalaxyData/Galaxy_%d_KinDat'% int(galaxy) + '.txt')
	pos = pos[:,0]

	out = np.zeros((len(pos), 6))

	out[:,0] = pos
	out[:,1] = np.median(data, axis = 0)
	out[:,2] = np.percentile(data, 16, axis = 0)
	out[:,3]= np.percentile(data, 84, axis = 0)
	out[:,4]= np.percentile(data, 2.5, axis = 0)
	out[:,5]= np.percentile(data, 97.5, axis = 0)

	return out

def create_sub(project_name, num_cores, timevar, workdir, anis, darkmatter, vsps, plummer, num_walkers, burn_in, steps, int_points):
	core = int(num_cores)
	num_allgal = glob.glob(workdir + '/KinPhotData/Galaxy*hdf5')
	galaxies = np.arange(len(num_allgal))

	for galaxy in galaxies:
		time = timevar

		prestr1 = workdir + project_name + '/Submissions/'
		prestr2 = workdir + project_name + '/OutErr/'

		f = open(prestr1 + 'Read_%d.sh' % galaxy, 'w')
		f.write("#!/bin/bash -l"+ '\n')
		f.write("#SBATCH --ntasks=%d" %core + '\n')
		f.write("#SBATCH -J Gal%d_" % galaxy  + "%s" % project_name + '\n')
		f.write("#SBATCH -e " + prestr2 + "Gal_%d_" %galaxy + "%s.err" % project_name + '\n')
		f.write("#SBATCH -o " + prestr2 + "Gal_%d_" %galaxy + "%s.out" % project_name + '\n')
		f.write("#SBATCH -p cosma-prince" + '\n')
		f.write("#SBATCH -A durham"+ '\n')
		f.write("#SBATCH --time=%02d" %time + ":00:00"+ '\n')
		f.write("module purge"+ '\n')
		f.write("module load python"+ '\n')
		f.write("module load gnu_comp/7.3.0 openmpi/3.0.1"+ '\n')
		f.write("module load gsl"+ '\n')

		f.write("python " + workdir +  project_name + '/Submissions/' + 'write_script.py' + " %s"%(workdir) + " %s" %(project_name) +  " %d"%(galaxy)  + " %s"%(core) + " %s"%(num_walkers) + " %s"%(burn_in) + " %s"%(steps)  + " %s"%(int_points) + " %s"%(darkmatter) + " %s"%(anis) +  " %s"%(vsps) + " %s"%(plummer)  )
		f.close()


def create_ana_sub(project_name, workdir, dm_option, beta_option, plummer_option, samples):
	core = 1
	num_allgal = glob.glob(workdir + '/KinPhotData/Galaxy*hdf5')
	galaxies = np.arange(len(num_allgal))
	

	for galaxy in galaxies:
		

		prestr1 = workdir + project_name + '/Analysis/Submissions/'
		prestr2 = workdir + project_name + '/Analysis/Submissions/'

		f = open(prestr1 + 'Read_%d.sh' % galaxy, 'w')
		f.write("#!/bin/bash -l"+ '\n')
		f.write("#SBATCH --ntasks=%d" %core + '\n')
		f.write("#SBATCH -J Ana%d_" % galaxy  + "%s" % project_name + '\n')
		f.write("#SBATCH -e " + prestr2 + "Ana_%d_" %galaxy + "%s.err" % project_name + '\n')
		f.write("#SBATCH -o " + prestr2 + "Ana_%d_" %galaxy + "%s.out" % project_name + '\n')
		f.write("#SBATCH -p cosma-prince" + '\n')
		f.write("#SBATCH -A durham"+ '\n')
		f.write("#SBATCH --time=03:00:00"+ '\n')
		f.write("module purge"+ '\n')
		f.write("module load python"+ '\n')
		f.write("module load gnu_comp/7.3.0 openmpi/3.0.1"+ '\n')
		f.write("module load gsl"+ '\n')


		
		if (dm_option == 'PL') and (beta_option == 'Baes') and (plummer_option == 'Plummer3'):
			starter = 10
			f.write("python " + code_dir +  '/AnalysisCode/' + 'dens_plummer.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy)) + ' ' + str(int(starter)) + ' ' + str(int(samples)) + '\n')
		if (dm_option == 'Zhao') and (beta_option == 'Baes') and (plummer_option == 'Plummer3'):
			starter = 9
			f.write("python " + code_dir +  '/AnalysisCode/' + 'dens_plummer.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy)) + ' ' + str(int(starter)) + ' ' + str(int(samples)) + '\n')
			

		if (dm_option == 'Zhao') and (beta_option == 'Baes'):
			f.write("python " + code_dir + '/AnalysisCode/' + 'dens_chains_read_baes2.py ' + workdir +' ' +  project_name  + ' '+ str(int(galaxy)) + ' ' + str(int(samples)) + '\n')
			f.write("python " + code_dir +  '/AnalysisCode/' + 'dens_chains_read_baes.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy))  + ' ' + str(int(samples)) + '\n')
		elif (dm_option == 'PL') and (beta_option == 'Baes'):
			f.write("python " + code_dir +   '/AnalysisCode/' + 'dens_chains_read_baes2PL.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy)) + ' ' + str(int(samples)) + '\n')
			f.write("python " + code_dir +  '/AnalysisCode/' + 'dens_chains_read_baesPL.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy))  + ' ' + str(int(samples)) + '\n')
		elif (dm_option == 'Zhao') and (beta_option == 'Const'):
			f.write("python " + code_dir + '/AnalysisCode/' + 'dens_chains_read.py ' + workdir +' ' +  project_name  + ' '+ str(int(galaxy)) + ' 5 ' + ' ' + str(int(samples)) +  '\n')
			f.write("python " + code_dir +  '/AnalysisCode/' + 'dens_chains_read_baes.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy))  + ' ' + str(int(samples)) + '\n')
		elif (dm_option == 'PL') and (beta_option == 'Const'):
			
			f.write("python " + code_dir +   '/AnalysisCode/' + 'dens_chains_read.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy)) +' 6 ' + ' ' + str(int(samples)) +  '\n')
			f.write("python " + code_dir +  '/AnalysisCode/' + 'dens_chains_read_baesPL.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy))  + ' ' + str(int(samples)) + '\n')

		if dm_option == 'PL':
			f.write("python " + code_dir +  '/AnalysisCode/' + 'sigvs_pl.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy)) +' ' + str(int(samples)) +  '\n')
		else:
			f.write("python " + code_dir +  '/AnalysisCode/' + 'sigvs_zhao.py ' + workdir +' ' +  project_name  + ' ' + str(int(galaxy)) + ' ' + str(int(samples)) +  '\n')			
		f.close()



def submit_jobs(workdir, project_name, all_gals):
	for galaxy in all_gals:
		os.system('sbatch ' + workdir + project_name + '/Submissions/' + 'Read_%d.sh' % galaxy)




def preprocess(workdir, all_gals):
	sys.path.append(workdir)
	from GravSphere import gal_input
	for galaxy in all_gals:
		gal_input.galaxy_data(int(galaxy), workdir + '/GalaxyData/', workdir + '/KinPhotData/')

	return 0





#workdir = open("/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/Read/workdir.txt", 'r')
#workdir = workdir.read()
#workdir = workdir.strip()

cwd = str(os.getcwd())

try:
	workdir = open(cwd + '/workdir.txt', 'r')
	workdir = workdir.read()
	workdir = workdir.strip()
except IOError:
	workdir = cwd 



code_dir = workdir

program = True

while program == True:
	print("""
		                                                   
,---.     |                   ,---.                
`---.,---.|---.,---.,---.,---.|  _.,---.,---..    ,
    ||   ||   ||---'|    |---'|   ||    ,---| \  / 
`---'|---'`   '`---'`    `---'`---'`    `---^  `'  
     |                                            


	""")



	print "Hi Anna, what would you like to do?"
	print "0) Preprocess"
	print "1) Create a new project" 
	print "2) Submit jobs"
	print "3) Analysis"
	print "4) Get limits"
	print "5) Plot chains"
	print "6) Plot triangle"
	print "quit -- a self-explanatory option"

	option = raw_input("Option : ")
	print option.strip()

	if option.strip() == "0":
		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ")
			exists = os.path.isdir(workdir)
			while exists == False:
				print "Directory does not exist, please type again"
				workdir = raw_input("Enter new name : ")
				exists = os.path.isdir(workdir)

				
			
			workdir = workdir				
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()


		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			pass

		else:
			print 'Directory name is not changed. Type correctly next time. :P'
			program = False
			continue
			
		print 'Current working directory is', workdir
		print "Type galaxies, separated by commas"
		res = raw_input('Type here: ')
		all_gals2 = res.split(',')
		all_gals = []
		for g in all_gals2:
			all_gals.append(int(g))
		all_gals = np.array(all_gals)

		print 'Generating input files'

		preprocess(workdir, all_gals)		

		print 'Done'
		program = True



	if option.strip() == "1":
		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ")
			exists = os.path.isdir(workdir)
			while exists == False:
				print "Directory does not exist, please type again"
				workdir = raw_input("Enter new name : ")
				exists = os.path.isdir(workdir)

				
			
				
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()

			project_name = raw_input("Please name your new project: " )
			exists = os.path.isdir(workdir + '/' + project_name)
			if exists == True:
				print 'This project already exists! Quitting!'
				break
			project_description = raw_input("Please write a short summary (or long if you want): ")		


			print 'Creating project directory'
			
			print 'Creating Submission files'
			num_cores = raw_input('How many cores? ')
			timevar = raw_input('How much time do you need (hours only)? ')
			plzh = raw_input('1) PowerLaw or 2) Zhao ? ')
			anis = raw_input('1) Baes or 2) Const ? ')
			vsps = raw_input('VSP ? y or n?')
			plummer = raw_input('1) Plummer Single  2) Plummer Components ')
			num_walkers = raw_input('How many walkers? ')
			burn_in = raw_input('Burn-in? ')
			steps = raw_input('Steps? ')
			int_points = raw_input('Integration points? ')

			if plzh == '1':
				plzh = 'PL'
			else:
				plzh = 'Zhao'
			if anis == '1':
				anis = 'Baes'
			else:
				anis = 'Const'
			if vsps == 'y':
				vsps = 'VSP'
			else:
				vsps = 'NOVSP'
			if plummer == '1':
				plummer = 'Plummer'
			else:
				plummer = 'Plummer3'

			os.system("mkdir " + workdir + '/' + project_name) 
		


			print 'Creating submissions directory'
			os.system("mkdir " + workdir + project_name + '/Submissions')
			print 'Copying codes'
			os.system("cp " + workdir + '/simps.py ' + workdir + project_name + '/Submissions/')
			os.system("cp " + workdir + '/_simps.so ' + workdir + project_name + '/Submissions/')
			os.system("cp " + workdir + '/priors.txt ' + workdir + project_name + '/Submissions/')
			os.system("cp " + workdir + '/write_script.py ' + workdir + project_name + '/Submissions/')

			print 'Creating OutErr directory'
			os.system("mkdir " + workdir  + project_name + '/OutErr/')
			print 'Creating Triangle directory'
			os.system("mkdir " + workdir  + project_name + '/Triangle/')


			
			create_sub(project_name, int(num_cores), int(timevar), workdir,  anis, plzh, vsps, plummer, num_walkers, burn_in, steps, int_points)		

			p_list = open(workdir + '/projects.txt', 'a+')
			p_list.write(project_name + '\n')
			p_list.close()
			d_file = open(workdir + project_name + '/ReadMe.txt', 'w')
			d_file.write(project_description)
			d_file.close()	



		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			project_name = raw_input("Please name your new project: " )
			exists = os.path.isdir(workdir + '/' + project_name)
			if exists == True:
				print 'This project already exists! Quitting!'
				break
			
			num_cores = raw_input('How many cores? ')
			timevar = raw_input('How much time do you need (hours only)? ')

			plzh = raw_input('1) PowerLaw or 2) Zhao ? ')
			anis = raw_input('1) Baes or 2) Const ? ')
			vsps = raw_input('VSP ? y or n? ')
			plummer = raw_input('1) Plummer Single  2) Plummer Components ')
			num_walkers = raw_input('How many walkers? ')
			burn_in = raw_input('Burn-in? ')
			steps = raw_input('Steps? ')
			int_points = raw_input('Integration points? ')
			if plzh == '1':
				plzh = 'PL'
			else:
				plzh = 'Zhao'
			if anis == '1':
				anis = 'Baes'
			else:
				anis = 'Const'
			if vsps == 'y':
				vsps = 'VSP'
			else:
				vsps = 'NOVSP'
			if plummer == '1':
				plummer = 'Plummer'
			else:
				plummer = 'Plummer3'


			print 'Creating project directory'
			print "mkdir " + workdir + '/' + project_name
			os.system("mkdir " + workdir + '/' + project_name) 

			project_description = raw_input("Please write a short summary (or long if you want): ")	
			
			os.system("mkdir " + workdir + project_name + '/Submissions')

			os.system("cp " + workdir + '/simps.py ' + workdir + project_name + '/Submissions/')
			os.system("cp " + workdir + '/_simps.so ' + workdir + project_name + '/Submissions/')
			os.system("cp " + workdir + '/priors.txt ' + workdir + project_name + '/Submissions/')
			os.system("cp " + workdir + '/write_script.py ' + workdir + project_name + '/Submissions/')

			print 'Creating submissions directory'
			
			print 'Creating OutErr directory'
			os.system("mkdir " + workdir  + project_name + '/OutErr/')
			print 'Creating Triangle directory'
			os.system("mkdir " + workdir  + project_name + '/Triangle/')
			print 'Creating Submission files'


			
			create_sub(project_name, int(num_cores), int(timevar), workdir,  anis, plzh, vsps, plummer, num_walkers, burn_in, steps, int_points)		

			p_list = open(workdir + '/projects.txt', 'a+')
			p_list.write(project_name + '\n')
			p_list.close()
			d_file = open(workdir + project_name + '/ReadMe.txt', 'w')
			d_file.write(project_description)
			d_file.close()	

				

		else:
			print "Fucked it, try again."
			program  = True 


	elif option.strip() == "2":
		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ")
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()
			valid = False
			while valid == False:
				print "Current projects:"
				try:
					project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
					
					if np.size(project_list) == 1:
						project_list = np.array([project_list], dtype = 'str')
					for pl in project_list:
						d_file = open(workdir + pl + '/ReadMe.txt', 'r')
						description = d_file.read()
						print pl
						print description
						print "********************************************"
				#print (contents)
				except IOError:
					print 'There are no projects. Lazy!'
				
				project_name = raw_input("What's the name of the project? " )
			
				if not os.listdir(workdir + project_name + '/Submissions'):
					print 'Submission directory is empty.'
					dec = raw_input("Try again? y or n? ")
					if dec == 'y':
						valid = False
					else:
						print 'Quitting'
						program = False
						break 	
				else:
					valid = True

					
					

					

					print 'Which galaxies would you like to submit?'
					print '1) All 	2) Specify'
					opt = raw_input('Option: ')
					if opt == '1':
						num_allgal = glob.glob(workdir + '/KinPhotData/Galaxy*hdf5')
						all_gals = np.arange(len(num_allgal))
						
						for gal in all_gals:
							exists = os.path.isdir(workdir  + project_name + '/Galaxy_%d/' %gal)
							if exists == False:
								os.system("mkdir " + workdir  + project_name + '/Galaxy_%d/' %gal)
						submit_jobs(workdir, project_name, all_gals)
						
						program = False
					elif opt == '2':
						print 'Please type galaxies, separated by commas'
						res = raw_input('Type here: ')
						all_gals2 = res.split(',')
						all_gals = []
						for g in all_gals2:
							all_gals.append(int(g))
						all_gals = np.array(all_gals)
						for gal in all_gals:
							exists = os.path.isdir(workdir  + project_name + '/Galaxy_%d/' %gal)
							if exists == False:
								os.system("mkdir " + workdir  + project_name + '/Galaxy_%d/' %gal)

						submit_jobs(workdir, project_name, all_gals)
						
						program = False
					else:
						print 'Fucked it , bye!'
						
						program = False

		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			valid = False
			while valid == False:
				print "Current projects:"
				try:
					project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
					
					if np.size(project_list) == 1:
						project_list = np.array([project_list], dtype = 'str')
					for pl in project_list:
						d_file = open(workdir + pl + '/ReadMe.txt', 'r')
						description = d_file.read()
						print pl
						print description
						print "********************************************"
				#print (contents)
				except IOError:
					print 'There are no projects. Lazy!'
				
				project_name = raw_input("What's the name of the project? " )
			
				if not os.listdir(workdir + project_name + '/Submissions'):
					print 'Submission directory is empty.'
					dec = raw_input("Try again? y or n? ")
					if dec == 'y':
						valid = False
					else:
						print 'Quitting'
						program = False
						break 	
				else:
					valid = True
					print 'Which galaxies would you like to submit?'
					print '1) All 	2) Specify'
					opt = raw_input('Option: ')
					if opt == '1':
						num_allgal = glob.glob(workdir + '/KinPhotData/Galaxy*hdf5')
						all_gals = np.arange(len(num_allgal))
						for gal in all_gals:
							exists = os.path.isdir(workdir  + project_name + '/Galaxy_%d/' %gal)
							if exists == False:
								os.system("mkdir " + workdir  + project_name + '/Galaxy_%d/' %gal)
						submit_jobs(workdir, project_name, all_gals)
						
						program = False
					elif opt == '2':
						print 'Please type galaxies, separated by commas'
						res = raw_input('Type here: ')
						all_gals2 = res.split(',')
						all_gals = []
						for g in all_gals2:
							all_gals.append(int(g))
						all_gals = np.array(all_gals)
						for gal in all_gals:
							exists = os.path.isdir(workdir  + project_name + '/Galaxy_%d/' %gal)
							if exists == False:
								os.system("mkdir " + workdir  + project_name + '/Galaxy_%d/' %gal)

						submit_jobs(workdir, project_name, all_gals)
						
						program = False
					else:
						print 'Fucked it , bye!'
						
						program = False
						
			

	

		else:
			print "Fucked it, try again."
			program  = True 
	elif option.strip() == "3":
		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ")
							
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()
			valid = False
			while valid == False:
				print "Current projects:"
				try:
					project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
					
					if np.size(project_list) == 1:
						project_list = np.array([project_list], dtype = 'str')
					for pl in project_list:
						d_file = open(workdir + pl + '/ReadMe.txt', 'r')
						description = d_file.read()
						print pl
						print description
						print "********************************************"
				#print (contents)
				except IOError:
					print 'There are no projects. Lazy!'
				
				project_name = raw_input("What's the name of the project? " )
			
				if not os.listdir(workdir + project_name + '/Submissions'):
					print 'Submission directory is empty.'
					dec = raw_input("Try again? y or n? ")
					if dec == 'y':
						valid = False
					else:
						print 'Quitting'
						program = False
						break 	
				else:
					valid = True

		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			valid = False
			while valid == False:
				print "Current projects:"
				try:
					project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
					
					if np.size(project_list) == 1:
						project_list = np.array([project_list], dtype = 'str')
					for pl in project_list:
						d_file = open(workdir + pl + '/ReadMe.txt', 'r')
						description = d_file.read()
						print pl
						print description
						print "********************************************"
				#print (contents)
				except IOError:
					print 'There are no projects. Lazy!'
				
				project_name = raw_input("What's the name of the project? " )
			
				if not os.listdir(workdir + project_name + '/Submissions'):
					print 'Submission directory is empty.'
					dec = raw_input("Try again? y or n? ")
					if dec == 'y':
						valid = False
					else:
						print 'Quitting'
						program = False
						break 	
				else:
					valid = True

		else:
			program = False
			continue

		if os.path.isdir(workdir+project_name + '/Analysis') == False:
			os.system("mkdir " + workdir  + project_name + '/Analysis')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Submissions')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Output')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Limits')


		samples = raw_input("How many samples?" )

		print "Generating submission scripts"
		
		
		prestr1 = workdir + project_name + '/Submissions/'	
		f_sub = open(prestr1 + 'Read_0.sh' , 'r')
		line = f_sub.readline()
		while line:
			input_opt = line.split()
			if input_opt[0] == 'python':
				dm_option = input_opt[-4]
				beta_option = input_opt[-3]
				plummer_option = input_opt[-1]
			line = f_sub.readline()		
		
		create_ana_sub(project_name, workdir, dm_option, beta_option, plummer_option, samples)		

		opt_sub = raw_input('Would you like to submit a job? y or n?')
		if opt_sub.strip() == 'y':
			opt_gal = raw_input("Which galaxies?  1) Specify    2) All")
			if opt_gal.strip() == '1':
				print "Type galaxies, separated by commas"
				res = raw_input('Type here: ')
				all_gals2 = res.split(',')
				all_gals = []
				for g in all_gals2:
					all_gals.append(int(g))
				all_gals = np.array(all_gals)

				for galaxy in all_gals:
					os.system('sbatch ' + workdir + project_name + '/Analysis/Submissions/' + 'Read_%d.sh' % galaxy)
			elif opt_gal.strip() == '2':
				num_allgal = glob.glob(workdir + '/KinPhotData/Galaxy*hdf5')
				all_gals = np.arange(len(num_allgal))
				for galaxy in all_gals:
					os.system('sbatch ' + workdir + project_name + '/Analysis/Submissions/' + 'Read_%d.sh' % galaxy)
			else:
				print "Option doesn't exist. Quitting."	
				program = False
				continue
		elif opt_sub.strip() == 'y':
			print 'Oh, I see how it is. Fine then.'	
			program = True
			continue
		else:
			print "Option doesn't exist. Quitting."	
			program = False
			continue


	elif option.strip() == '4':
		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ")
							
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()
		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			pass
		else:
			print 'Not valid'
			continue
			

		
		valid = False
		while valid == False:
			print "Current projects:"
			try:
				project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
				
				if np.size(project_list) == 1:
					project_list = np.array([project_list], dtype = 'str')
				for pl in project_list:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
					print "********************************************"
			#print (contents)
			except IOError:
				print 'There are no projects. Lazy!'
			
			project_name = raw_input("What's the name of the project? " )
		
			if not os.listdir(workdir + project_name + '/Submissions'):
				print 'Submission directory is empty.'
				dec = raw_input("Try again? y or n? ")
				if dec == 'y':
					valid = False
					#continue
				else:
					print 'Quitting'
					program = False
					valid = True
					#continue 	
			else:
				valid = True

		if program == False:
			break 
		
		
		print 'Which galaxies would you like to submit?'
		print '1) All 	2) Specify'
		opt = raw_input('Option: ')
		if opt == '1':
			num_allgal = glob.glob(workdir + '/KinPhotData/Galaxy*hdf5')
			all_gals = np.arange(len(num_allgal))
						
						
						
						
		elif opt == '2':
			print 'Please type galaxies, separated by commas'
			res = raw_input('Type here: ')


		
			all_gals2 = res.split(',')
			all_gals = []
			for g in all_gals2:
				all_gals.append(int(g))
			all_gals = np.array(all_gals)
		
		prestr1 = workdir + project_name + '/Submissions/'	
			
		f_sub = open(prestr1 + 'Read_%d.sh' %int(all_gals[0]) , 'r')
                line = f_sub.readline()
                while line:
                	input_opt = line.split()
                	if input_opt[0] == 'python':
                               dm_option = input_opt[-4]
                
                        line = f_sub.readline()


		for galaxy in all_gals:
			rhs = np.loadtxt(workdir + '/GalaxyData/Galaxy_%d_Rhalf.txt' %int(galaxy))
			#rsel, = np.where(rhs[:,0] == int(galaxy))
			rh = rhs

			bin_edges = np.array([0.125, 0.25, 0.50, 1, 2 , 4 , 8])*rh
			mid_bin = (np.log10(bin_edges[1:]) + np.log10(bin_edges[:-1]))/2.
			mid_bin = 10**mid_bin

			tot_bins = np.array([bin_edges[0], mid_bin[0], bin_edges[1], bin_edges[1], mid_bin[1], bin_edges[2],bin_edges[2], mid_bin[2], bin_edges[3], bin_edges[3], mid_bin[3], bin_edges[4], bin_edges[4], mid_bin[4], bin_edges[5]])

			r = np.logspace(-2,1.5, 25)

			
			if dm_option == 'PL':
				try:
					f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_LogLog'%galaxy + '.txt', invalid_raise = False)
					res = get_lims_loglog(f, tot_bins)
					np.savetxt(workdir + project_name + '/Analysis/Limits/' + 'Galaxy_%d_LogLog'% galaxy + 'Lims.txt', res)
				except IOError:
					print 'Slope output not found'
			else:
				try:
					f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_LogLog'%galaxy + '.txt' , invalid_raise = False)
					res = get_lims(f,r)
					np.savetxt(workdir + project_name + '/Analysis/Limits/' + 'Galaxy_%d_LogLog'% galaxy + 'Lims.txt', res)
				except IOError:
					print 'Slope output not found'
			print 'Galaxy', galaxy, 'Slopes done'
			
			try:	
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_Beta'%galaxy + '.txt', invalid_raise = False)
				res = get_lims(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ 'Galaxy_%d_Beta'% galaxy +  'Lims.txt', res)
			except IOError:
				print 'Beta output not found'
			print 'Galaxy', galaxy, 'Beta done'

			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_Density'%galaxy  +'.txt' , invalid_raise = False)
					
				res = get_lims(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/' + 'Galaxy_%d_Density'% galaxy +  'Lims.txt', res)
			except IOError:
				print 'Density output not found'
			print 'Galaxy', galaxy, 'Density done'

			
			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_Mass'%galaxy + '.txt', invalid_raise = False)
				res = get_lims(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ 'Galaxy_%d_Mass'% galaxy + 'Lims.txt', res)
			except IOError:
				print 'Mass output not found'

			print 'Galaxy', galaxy, 'Mass done'

			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_Plummer'%galaxy + '.txt', invalid_raise = False)
				res = get_lims(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ 'Galaxy_%d_Plummer'% galaxy + 'Lims.txt', res)
			except IOError:
				print 'Plummer output not found'

			print 'Galaxy', galaxy, 'Plummer done'	

			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_SigLos'%galaxy + '.txt', invalid_raise = False)
				res = get_lims_sig(f, workdir, galaxy)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ 'Galaxy_%d_SigLos'% galaxy + 'Lims.txt', res)
			except IOError:
				print 'SigP output not found'

			print 'Galaxy', galaxy, 'SigmaP + VSP done'


	elif (option.strip() == '5'):
		valid = False
		while valid == False:
			print "Current projects:"
			try:
				project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
				
				if np.size(project_list) == 1:
					project_list = np.array([project_list], dtype = 'str')
				for pl in project_list:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
					print "********************************************"
			#print (contents)
			except IOError:
				print 'There are no projects. Lazy!'
			
			project_name = raw_input("What's the name of the project? " )
		
			if not os.listdir(workdir + project_name + '/Submissions'):
				print 'Submission directory is empty.'
				dec = raw_input("Try again? y or n? ")
				if dec == 'y':
					valid = False
				else:
					print 'Quitting'
					program = False
					break 	
			else:
				valid = True

		print "Type galaxies, separated by commas"
		res = raw_input('Type here: ')
		all_gals2 = res.split(',')
		all_gals = []
		for g in all_gals2:
			all_gals.append(int(g))
		all_gals = np.array(all_gals)

		prestr1 = workdir + project_name + '/Submissions/'	
		f_sub = open(prestr1 + 'Read_0.sh' , 'r')
		line = f_sub.readline()
		while line:
			input_opt = line.split()
			if input_opt[0] == 'python':
				dm_option = input_opt[-4]
				beta_option = input_opt[-3]
				plummer_option = input_opt[-1]
			line = f_sub.readline()		

		param_names = []
		nparams = 0
		if dm_option == 'PL':
			param_names.extend([r"$\rho_0$", r"$\gamma_0$", r"$\gamma_1$", r"$\gamma_2$",r"$\gamma_3$",r"$\gamma_4$"])
			nparams = nparams + 6
		else:
			param_names.extend([r"$\rho_s$", r"$r_s$", r"$\alpha$", r"$\beta$", r"$\gamma$"])
			nparams = nparams + 5
		if beta_option == 'Baes':
			nparams = nparams + 4
			param_names.extend([r"$\beta_0$", r"$\beta_{\infty}$", r"$r_a$", r"$\eta$"])
		else:
			param_names.extend([r"$\beta_0$"])
			nparams = nparams + 1
		if plummer_option == 'Plummer3':
			nparams = nparams + 6
			param_names.extend([r"$m_1$", r"$a_1$", r"$m_2$", r"$a_2$", r"$m_3$", r"$a_3$"])
		else:
			
			nparams = nparams + 0

		for galaxy in all_gals:
			data = np.loadtxt(workdir + project_name + '/Galaxy_%d/' %galaxy + 'Galaxy_%d_Chains' % galaxy + project_name + '.txt')
			plot_chains(data, param_names, nparams, workdir, project_name, galaxy)
			




	elif (option.strip() == '6'):
		valid = False
		while valid == False:
			print "Current projects:"
			try:
				project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
				
				if np.size(project_list) == 1:
					project_list = np.array([project_list], dtype = 'str')
				for pl in project_list:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
					print "********************************************"
			#print (contents)
			except IOError:
				print 'There are no projects. Lazy!'
			
			project_name = raw_input("What's the name of the project? " )
		
			if not os.listdir(workdir + project_name + '/Submissions'):
				print 'Submission directory is empty.'
				dec = raw_input("Try again? y or n? ")
				if dec == 'y':
					valid = False
				else:
					print 'Quitting'
					program = False
					break 	
			else:
				valid = True

		print "Type galaxies, separated by commas"
		res = raw_input('Type here: ')
		all_gals2 = res.split(',')
		all_gals = []
		for g in all_gals2:
			all_gals.append(int(g))
		all_gals = np.array(all_gals)

		prestr1 = workdir + project_name + '/Submissions/'	
		f_sub = open(prestr1 + 'Read_0.sh' , 'r')
		line = f_sub.readline()
		while line:
			input_opt = line.split()
			if input_opt[0] == 'python':
				dm_option = input_opt[-4]
				beta_option = input_opt[-3]
				plummer_option = input_opt[-1]
			line = f_sub.readline()		

		param_names = []
		nparams = 0
		if dm_option == 'PL':
			param_names.extend([r"$\rho_0$", r"$\gamma_0$", r"$\gamma_1$", r"$\gamma_2$",r"$\gamma_3$",r"$\gamma_4$"])
			nparams = nparams + 6
		else:
			param_names.extend([r"$\rho_s$", r"$r_s$", r"$\alpha$", r"$\beta$", r"$\gamma$"])
			nparams = nparams + 5
		if beta_option == 'Baes':
			nparams = nparams + 4
			param_names.extend([r"$\beta_0$", r"$\beta_{\infty}$", r"$r_a$", r"$\eta$"])
		else:
			param_names.extend([r"$\beta_0$"])
			nparams = nparams + 1

		for galaxy in all_gals:
			data = np.loadtxt(workdir + project_name + '/Galaxy_%d/' %galaxy + 'Galaxy_%d_Chains' % galaxy + project_name + '.txt')
			plot_triangle(data, param_names, nparams, workdir, project_name, galaxy)


	

	elif (option.strip() == "quit") or (option.strip() == "q"):
		program = False		


	elif option.strip() == '10':
		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ")
							
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()
		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			pass
		else:
			print 'Not valid'
			continue
			

		
		valid = False
		while valid == False:
			print "Current projects:"
			try:
				project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
				
				if np.size(project_list) == 1:
					project_list = np.array([project_list], dtype = 'str')
				for pl in project_list:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
					print "********************************************"
			#print (contents)
			except IOError:
				print 'There are no projects. Lazy!'
			
			project_name = raw_input("What's the name of the project? " )
		
			if not os.listdir(workdir + project_name + '/Submissions'):
				print 'Submission directory is empty.'
				dec = raw_input("Try again? y or n? ")
				if dec == 'y':
					valid = False
					#continue
				else:
					print 'Quitting'
					program = False
					valid = True
					#continue 	
			else:
				valid = True

		if program == False:
			break 
		
		
		print 'Which galaxies would you like to submit?'
		print '1) All 	2) Specify'
		opt = raw_input('Option: ')
		if opt == '1':
			num_allgal = glob.glob(workdir + '/KinPhotData/Galaxy*hdf5')
			all_gals = np.arange(len(num_allgal))
						
						
						
						
		elif opt == '2':
			print 'Please type galaxies, separated by commas'
			res = raw_input('Type here: ')


		
			all_gals2 = res.split(',')
			all_gals = []
			for g in all_gals2:
				all_gals.append(int(g))
			all_gals = np.array(all_gals)
		
		prestr1 = workdir + project_name + '/Submissions/'	
			
		f_sub = open(prestr1 + 'Read_%d.sh' %int(all_gals[0]) , 'r')
                line = f_sub.readline()
                while line:
                	input_opt = line.split()
                	if input_opt[0] == 'python':
                               dm_option = input_opt[-4]
                
                        line = f_sub.readline()


		for galaxy in all_gals:
			rhs = np.loadtxt(workdir + '/GalaxyData/Galaxy_%d_Rhalf.txt' %int(galaxy))
			#rsel, = np.where(rhs[:,0] == int(galaxy))
			rh = rhs

			bin_edges = np.array([0.125, 0.25, 0.50, 1, 2 , 4 , 8])*rh
			mid_bin = (np.log10(bin_edges[1:]) + np.log10(bin_edges[:-1]))/2.
			mid_bin = 10**mid_bin

			tot_bins = np.array([bin_edges[0], mid_bin[0], bin_edges[1], bin_edges[1], mid_bin[1], bin_edges[2],bin_edges[2], mid_bin[2], bin_edges[3], bin_edges[3], mid_bin[3], bin_edges[4], bin_edges[4], mid_bin[4], bin_edges[5]])

			r = np.logspace(-2,1.5, 25)

			
			
			
			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + 'Galaxy_%d_Mass'%galaxy + '.txt', invalid_raise = False)
				res = get_lims_all(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ 'Galaxy_%d_Mass'% galaxy + 'LimsAll.txt', res)
			except IOError:
				print 'Mass output not found'

			print 'Galaxy', galaxy, 'Mass done'

			
