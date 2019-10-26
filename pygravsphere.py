import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import corner
import glob
import datetime
from multiprocessing import cpu_count

cwd = str(os.getcwd())

try:
	workdir = open(cwd + '/workdir.txt', 'r')
	workdir = workdir.read()
	workdir = workdir.strip()
	
		

except IOError:
	print 'workdir.txt file does not exist. Creating file.'
	os.system("touch workdir.txt")	
	workdir = open(cwd + '/workdir.txt', 'r')
	workdir = workdir.read()
	workdir = workdir.strip()

if not workdir:
	print "workdir.txt is empty. Please update workdir.txt or pyGravSphere will operate in current directory."	
	workdir = open(cwd + '/workdir.txt', 'w')
	workdir.write(cwd + '\n')
	workdir.close()
	workdir = cwd

if workdir[-1] != '/':
		workdir  = workdir + '/'


codedir = str(os.environ["GravSpherePath"])
if codedir[-1] != '/':
		codedir  = codedir + '/'
	
sys.path.append(codedir)
from GStools import gsTools


gsTools.checkdirs(workdir, codedir)



program = True

while program == True:

	

	banner = open(codedir + '/banner.txt', 'r')
	banner = banner.read()
	print (banner)

	try:
		projects = np.loadtxt(workdir + '/projects.txt',dtype = 'str', ndmin = 0)
		updated = []
		for p in range(0, len(projects)):
			exists = os.path.isdir(workdir + '/' + projects[p])
			if exists == True:
				updated.append(str(projects[p]))
		projects = updated
		np.savetxt(workdir + '/projects.txt',updated, fmt = "%s")
				
	except IOError:	
		print "There are no current projects"
	


	print "Hello galactic dynamicist, what would you like to do?"
	print "0) Preprocess data"
	print "1) Create a new project" 
	print "2) Submit jobs"
	print "3) Analysis"
	print "4) Plot chain convergence"
	print "5) Create a corner plot"
	print "6) Get mass profile percentiles"
	print "quit/q -- a self-explanatory option"

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
		print 'Which galaxies would you like to preprocess?'
		gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 0,dtype = 'str')
		if len(gal_list) == 0:
			print "There are no galaxies, Quitting."
			break
		else:
			for it in gal_list:
				print it
		print '1) All 	2) Specify'
		opt = raw_input('Option: ')
		if opt == '1':
			all_gals = gal_list
						
						
						
						
		elif opt == '2':
			print 'Please type galaxies, separated by commas'
			res = raw_input('Type here: ')


		
			all_gals2 = res.split(',')
			all_gals = []
			for g in all_gals2:
				all_gals.append(g)

		print 'Generating input files'

		gsTools.preprocess(workdir, codedir, all_gals)		

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


		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			print 'Alright then'
				
		else:	
			print "I don't know what that means. Quitting."	
			break	

		project_name = raw_input("Please name your new project: " )
		exists = os.path.isdir(workdir + '/' + project_name)
		if exists == True:
			print 'This project already exists! Quitting!'
			break
		project_description = raw_input("Please write a short summary (or long if you want): ")		


		print 'Creating project directory'
		
		print 'Creating Submission files'
		print 'How would you like to run MCMC?'
				
		mpi_opt = raw_input("1) Serial 	2) Multiprocessing  3) MPI (if you're running on a cluster)")
		if mpi_opt == '3':
			num_cores = int(raw_input('How many cores? '))
			timevar = float(raw_input('How much time do you need (in hours)? '))
			timevar = str(datetime.timedelta(hours = timevar))
		elif mpi_opt == '2':
			ncpu = cpu_count()
			num_cores = int(raw_input('How many processes out of %d that you have?' %ncpu))
			timevar = None
			

		else:
			num_cores = None
			timevar = None

		plzh = raw_input('1) PowerLaw or 2) Zhao ? ')
		anis = raw_input('1) Baes or 2) Const ? ')
		vsps = raw_input('VSP ? y or n?')
		plummer = raw_input('1) Plummer Single  2) Plummer Components ')
		num_walkers = raw_input('How many walkers? ')
		burn_in = raw_input('Burn-in? ')
		steps = raw_input('Steps? ')
		int_points = raw_input('Integration points? ')


		param_list = []

		if plzh == '1':
			plzh = 'PL'
			param_list.extend(['rho0', 'gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4'])
		else:
			plzh = 'Zhao'
			param_list.extend(['rhos','rs','alpha','beta','gamma'])
		if anis == '1':
			anis = 'Baes'
			param_list.extend(['beta0', 'betainf', 'ra','eta'])
		else:
			anis = 'Const'
			param_list.extend(['beta0'])
		if vsps == 'y':
			vsps = 'VSP'
		else:
			vsps = 'NOVSP'
		if plummer == '1':
			plummer = 'Plummer'
		else:
			plummer = 'Plummer3'
			param_list.extend(['m1','a1','m2','a2','m3','a3'])

		param_list.extend(['mstar'])

		os.system("mkdir " + workdir + '/' + project_name) 
	


		print 'Creating submissions directory'
		os.system("mkdir " + workdir + project_name + '/Submissions')
		fopts = open(workdir + project_name + '/options.txt', 'w')
		fopts.write(plzh + '\t' + anis + '\t' + plummer + '\t' + steps + '\t' + burn_in + '\t' + num_walkers)
		fopts.close()



		print 'Copying codes'
		
		os.system("cp " + codedir + '/priors.txt ' + workdir + project_name + '/Submissions/')
		priors = np.loadtxt(workdir + project_name + '/Submissions/priors.txt', dtype = 'str', ndmin = 1)
		my_params = np.where(np.in1d(priors[:,0], np.array(param_list, dtype = 'str')))		
		priors = np.savetxt(workdir + project_name + '/Submissions/priors.txt', priors[my_params], fmt = "%s")
		np.savetxt(workdir + project_name + '/Submissions/gamsmooth.txt', np.array([1.]))
		print 'Creating OutErr directory'
		os.system("mkdir " + workdir  + project_name + '/OutErr/')
		print 'Creating output directories'
		all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		for gal in all_gals:
			out_exists = os.path.isdir(workdir + '/' + project_name + '/%s' % gal)
			if out_exists == False:
				os.system("mkdir {}/{}/{}".format(workdir, project_name, gal))

		print 'Generating submission scripts'
		gsTools.create_sub(project_name, num_cores, timevar, workdir,codedir,  anis, plzh, vsps, plummer, num_walkers, burn_in, steps, int_points, mpi_opt)		

		p_list = open(workdir + '/projects.txt', 'a+')
		p_list.write(project_name + '\n')
		p_list.close()
		d_file = open(workdir + project_name + '/ReadMe.txt', 'w')
		d_file.write(project_description)
		d_file.close()	



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
		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			print 'Alright then!'
			valid = False
		else:
			print "What was that???"
			break
		
		



		while valid == False:
			print "Current projects:"
			try:
				project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
			except IOError:
				print 'There are no projects. Lazy!'
				valid = True
				break
				
			if np.size(project_list) == 1:
				project_list = np.array([project_list], dtype = 'str')
			for pl in project_list:
				try:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
			#print (contents)
			
			
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

				
				
				sub_command = open('sub_command.txt', 'r')
				sub_com = sub_command.read()
						

				print 'Which galaxies would you like to submit?'
				gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 0,dtype = 'str')
				if len(gal_list) == 0:
					print "There are no galaxies, Quitting."
					break
				else:
					for it in gal_list:
						print it
				print '1) All 	2) Specify'
				opt = raw_input('Option: ')
				if opt == '1':
					all_gals = np.loadtxt(workdir + '/galaxy_list.txt', ndmin = 1, dtype = 'str')
					for gal in all_gals:
						exists = os.path.isdir(workdir  + project_name + '/%s/' %gal)
						if exists == False:
							os.system("mkdir " + workdir  + project_name + '/%s/' %gal)
					gsTools.submit_jobs(workdir, project_name, all_gals, sub_com)
					
					program = False
				elif opt == '2':
					print 'Please type galaxies, separated by commas'
					res = raw_input('Type here: ')
					all_gals2 = res.split(',')
					all_gals = []
					for g in all_gals2:
						all_gals.append(g)
					
					for gal in all_gals:
						exists = os.path.isdir(workdir  + project_name + '/%s/' %gal)
						if exists == False:
							os.system("mkdir " + workdir  + project_name + '/%s/' %gal)

					gsTools.submit_jobs(workdir, project_name, all_gals, sub_com)
					
					program = False
				else:
					print "What's that???"
					
					program = False
		
		
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
		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			print 'Alright then!'
			valid = False
		else:
			print "What was that??? Try again."
			break
		while valid == False:
			print "Current projects:"
			try:
				project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
			except IOError:
				print 'There are no projects. Lazy!'
				valid = True
				break
				
			if np.size(project_list) == 1:
				project_list = np.array([project_list], dtype = 'str')
			for pl in project_list:
				try:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
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

		

		if os.path.isdir(workdir+project_name + '/Analysis') == False:
			os.system("mkdir " + workdir  + project_name + '/Analysis')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Submissions')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Output')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Limits')


		

		print "Generating submission scripts"
			

		foptions = open(workdir + project_name + '/options.txt', 'r')
		input_opt = (foptions.readline()).split()
		dm_option = input_opt[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]
		steps = int(input_opt[3])
		burn_in = int(input_opt[4])
		nwalkers = int(input_opt[5])
		
		mpi_opt = raw_input("Running on a batch system y or n? " )
		if mpi_opt == 'y':
			timevar = int(raw_input("How much time do you need? (in hours)"))
			timevar = str(datetime.timedelta(hours = timevar))
		else:
			timevar = None	

		print "Do you need to cut chains? "
		print "It looks like you ran %d" %steps
		cut_off = raw_input("How many first steps fould you like to discard for each chain? ")
		print "Ignore all chains with (chi squared > N x best chi squared)? N = 10 recommended for no particular reason."
		chi_cut = raw_input("Enter N = ")
		print "For the radial profile limits, what are the radial ranges in kpc? The intervals will be in log."
		min_rad = raw_input("Minimum radius = ")
		max_rad = raw_input("Maximum radius = ")
		points = raw_input("How many log-spaced intervals? ")
		samples = raw_input("How many samples out of ~ %d that you ran? " %int((steps-burn_in -cut_off)*nwalkers))


		
		gsTools.create_ana_sub(project_name, workdir, codedir, dm_option, beta_option, plummer_option, samples, mpi_opt, cut_off, chi_cut, min_rad, max_rad, points,timevar)		

		opt_sub = raw_input('Would you like to submit a job? y or n?')
		if opt_sub.strip() == 'y':

			sub_command = open('sub_command.txt', 'r')
			sub_com = sub_command.read()


			gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 0,dtype = 'str')
			if len(gal_list) == 0:
				print "There are no galaxies, Quitting."
				break
			else:
				for it in gal_list:
					print it

			opt_gal = raw_input("Which galaxies?  1) Specify    2) All")
			if opt_gal.strip() == '1':
				print "Type galaxies, separated by commas"
				res = raw_input('Type here: ')
				all_gals2 = res.split(',')
				all_gals = []
				for g in all_gals2:
					all_gals.append(g)
				
				for galaxy in all_gals:
					if sub_com == False:
						os.system('./' + workdir + project_name + '/Analysis/Submissions/' + '%s.sh' % galaxy)
					else:
						os.system(sub_com + ' ' + workdir + project_name + '/Analysis/Submissions/' + '%s.sh' % galaxy)
			elif opt_gal.strip() == '2':
				all_gals = np.loadtxt(workdir + '/galaxy_list.txt', ndmin = 1, dtype = 'str')
				for galaxy in all_gals:
					if sub_com == False:
						os.system('./' + workdir + project_name + '/Analysis/Submissions/' + '%s.sh' % galaxy)
					else:
						os.system(sub_com + ' ' + workdir + project_name + '/Analysis/Submissions/' + '%s.sh' % galaxy)
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


	

	elif (option.strip() == '4'):
		valid = False
		while valid == False:
			print "Current projects:"
			try:
				project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
			except IOError:
				print 'There are no projects. Lazy!'
				valid = True
				break
				
			if np.size(project_list) == 1:
				project_list = np.array([project_list], dtype = 'str')
			for pl in project_list:
				try:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
			
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

		print 'Which galaxies would you like to plot?'
		gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 0,dtype = 'str')
		if len(gal_list) == 0:
			print "There are no galaxies, Quitting."
			break
		else:
			for it in gal_list:
				print it
		print '1) All 	2) Specify'
		opt = raw_input('Option: ')
		if opt == '1':
			all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1, dtype = 'str')
						
						
						
						
		elif opt == '2':
			print 'Please type galaxies, separated by commas'
			res = raw_input('Type here: ')


		
			all_gals2 = res.split(',')
			all_gals = []
			for g in all_gals2:
				all_gals.append(g)
		

		prestr1 = workdir + project_name + '/Submissions/'	
		foptions = open(workdir + project_name + '/options.txt', 'r')
		input_opt = (foptions.readline()).split()
		dm_option = input_opt[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]

		for galaxy in all_gals:
			data = np.loadtxt(workdir + project_name + '/%s/' %galaxy + '%s_Chains' % galaxy + project_name + '.txt')
			gsTools.plot_chains(data, workdir, project_name, galaxy)
			




	elif (option.strip() == '5'):
		valid = False
		while valid == False:
			print "Current projects:"
			try:
				project_list = np.genfromtxt(workdir + '/projects.txt', dtype = 'str')
			except IOError:
				print 'There are no projects. Lazy!'
				valid = True
				break
				
			if np.size(project_list) == 1:
				project_list = np.array([project_list], dtype = 'str')
			for pl in project_list:
				try:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
			
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

		print 'Which galaxies would you like to plot?'
		gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 0,dtype = 'str')
		if len(gal_list) == 0:
			print "There are no galaxies, Quitting."
			break
		else:
			for it in gal_list:
				print it
		print '1) All 	2) Specify'
		opt = raw_input('Option: ')
		if opt == '1':
			all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
			
		elif opt == '2':
			print 'Please type galaxies, separated by commas'
			res = raw_input('Type here: ')


		
			all_gals2 = res.split(',')
			all_gals = []
			for g in all_gals2:
				all_gals.append(g)

		prestr1 = workdir + project_name + '/Submissions/'	
		foptions = open(workdir + project_name + '/options.txt', 'r')
		input_opt = (foptions.readline()).split()
		dm_option = input_opt[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]
			

		
		for galaxy in all_gals:
			data = np.loadtxt(workdir + project_name + '/%s/' %galaxy + '%s_Chains' % galaxy + project_name + '.txt')
			gsTools.plot_triangle(data, workdir, project_name, galaxy)


	

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
			except IOError:
				print 'There are no projects. Lazy!'
				valid = True
				break
				
			if np.size(project_list) == 1:
				project_list = np.array([project_list], dtype = 'str')
			for pl in project_list:
				try:
					d_file = open(workdir + pl + '/ReadMe.txt', 'r')
					description = d_file.read()
					print pl
					print description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
			
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
			all_gals = np.loadtxt(workdir + '/galaxy_list.txt', ndmin = 1, dtype = 'str')
			
						
		elif opt == '2':
			print 'Please type galaxies, separated by commas'
			res = raw_input('Type here: ')


		
			all_gals2 = res.split(',')
			all_gals = []
			for g in all_gals2:
				all_gals.append(g)
			all_gals = np.array(all_gals)
		
		prestr1 = workdir + project_name + '/Submissions/'	
			
		foptions = open(workdir + project_name + '/options.txt', 'r')
		input_opt = (foptions.readline()).split()
		dm_option = input_opt[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]

		for galaxy in all_gals:
			rhs = np.loadtxt(workdir + '/GalaxyData/%s_Rhalf.txt' %int(galaxy))
			#rsel, = np.where(rhs[:,0] == int(galaxy))
			rh = rhs

			bin_edges = np.array([0.125, 0.25, 0.50, 1, 2 , 4 , 8])*rh
			mid_bin = (np.log10(bin_edges[1:]) + np.log10(bin_edges[:-1]))/2.
			mid_bin = 10**mid_bin

			tot_bins = np.array([bin_edges[0], mid_bin[0], bin_edges[1], bin_edges[1], mid_bin[1], bin_edges[2],bin_edges[2], mid_bin[2], bin_edges[3], bin_edges[3], mid_bin[3], bin_edges[4], bin_edges[4], mid_bin[4], bin_edges[5]])

			r = np.logspace(-2,1.5, 25)

			
			
			
			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_Mass'%galaxy + '.txt', invalid_raise = False)
				res = gsTools.get_lims_all(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Mass'% galaxy + 'LimsAll.txt', res)
			except IOError:
				print 'Mass output not found'

			print 'Galaxy', galaxy, 'Mass done'

			
