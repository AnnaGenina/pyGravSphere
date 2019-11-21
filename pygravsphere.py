import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import corner
import glob
import datetime
from multiprocessing import cpu_count
import warnings

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
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			projects = np.loadtxt(workdir + '/projects.txt',dtype = 'str', ndmin = 1)
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
			workdir = raw_input("Enter new name : ").strip()
			exists = os.path.isdir(workdir)
			while exists == False:
				print "Directory does not exist, please type again"
				workdir = raw_input("Enter new name : ").strip()
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
		gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(gal_list) == 0:
			print "There are no galaxies in galaxy_list, Quitting."
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
				if g in gal_list:
					all_gals.append(g)
				else:
					print g, "This galaxy does not exist. Did you add it to galaxy_list.txt?"					
					valid_galname = False
					while valid_galname == False:
						retype = raw_input("This galaxy does not exist. Retype the name of this galaxy or leave blank: ")
						retype = retype.strip()
						if (retype in gal_list) or retype == "":
							valid_galname = True
						else:
							valid_galname = False

		print 'Generating input files'

		if len(all_gals) != 0: 	
			gsTools.preprocess(workdir, codedir, all_gals)		
		else:
			print "No valid galaxies"
		print 'Done'
		program = True
		continue



	if option.strip() == "1":

		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ").strip()
			exists = os.path.isdir(workdir)
			while exists == False:
				print "Directory does not exist, please type again"
				workdir = raw_input("Enter new name : ").strip()
				exists = os.path.isdir(workdir)
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()


		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			print 'Alright then'
				
		else:	
			print "I don't know what that means. Quitting."	
			break	
		
		all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(all_gals) == 0:
			print "There are no galaxies in galaxy_list"
			sys.exit()
		gsTools.check_galdata(all_gals,workdir)

		project_name = raw_input("Please name your new project: " )
		exists = os.path.isdir(workdir + '/' + project_name)
		if exists == True:
			print 'This project already exists! Quitting!'
			break
		project_description = raw_input("Please write a short summary (or long if you want): ")		


		print 'Creating project directory'
		
		print 'Creating Submission files'
		print 'How would you like to run MCMC?'
				
		mpi_opt = raw_input("1) Serial 	2) Multiprocessing  3) MPI (if you're running on a cluster) ")
		if mpi_opt == '3':
			num_cores = int(raw_input('How many cores? '))
			timevar = float(gsTools.check_float('How much time do you need (in hours)? '))
			timevar = str(datetime.timedelta(hours = timevar))
		elif mpi_opt == '2':
			ncpu = cpu_count()
			num_cores = int(gsTools.check_float('How many processes out of %d that you have? ' %ncpu))
			timevar = None
			

		else:
			num_cores = None
			timevar = None

		sub_scripts = open(codedir + "/sub_script.txt", 'r')
		sub_script = sub_scripts.read()
		sub_script = sub_script.strip()
		sub_scripts.close()
		if (not sub_script) and (mpi_opt == '3'):
			print "You submission script template is empty. Cannot submit to the queue."
			print "Edit ", codedir + "/sub_script.txt and come back."
			print "Exiting now"
			program = True
			continue


		valid_params = False
		while valid_params == False:

			valid_params = True

			param_list = []

			plzh = raw_input('1) PowerLaw or 2) Zhao ? ')
			if (plzh == '1') or (plzh == 'PowerLaw'):
				plzh = 'PL'
				param_list.extend(['rho0', 'gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4'])
			elif (plzh == '2') or (plzh == 'Zhao'):
				plzh = 'Zhao'
				param_list.extend(['rhos','rs','alpha','beta','gamma'])
			else:
				print "Option %s does not exist" %plzh
				valid_params = False
				continue

			anis = raw_input('1) Baes or 2) Constant ? ')
			if (anis == '1') or (anis == 'Baes'):
				anis = 'Baes'
				param_list.extend(['beta0', 'betainf', 'ra','eta'])
			elif (anis == '2') or (anis == 'Constant'):
				anis = 'Const'
				param_list.extend(['beta0'])
			else:
				print "Option %s does not exist" %anis
				valid_params = False
				continue

			vsps = raw_input('VSP ? y or n? ')
			if vsps == 'y':
				vsps = 'VSP'
			elif vsps == 'n':
				vsps = 'NOVSP'
			else:
				print "Option %s for vsps does not exist" %vsps
				valid_params = False
				continue

			plummer = raw_input('1) Plummer  2) ThreePlummers ')

			if (plummer == '1') or (plummer == 'Plummer'):
				plummer = 'Plummer'
			elif (plummer == '2') or (plummer == 'ThreePlummers'):
				plummer = 'Plummer3'
				param_list.extend(['m1','a1','m2','a2','m3','a3'])
			else:
				print "Option %s does not exist" 
				valid_params = False	
				continue

			param_list.extend(['mstar'])

		correct_walkers = False
		while correct_walkers == False:
			num_walkers = (int(gsTools.check_float('How many walkers? ')))
			if len(param_list)*2 >= num_walkers:
				print "Number of walkers needs to be more than twice the number of dimensions (i.e more than %d)" %2*len(param_list)
				correct_walkers = False
				
			else:
				correct_walkers = True	
		num_walkers = str(num_walkers)

		
		burn_in = str(int(gsTools.check_float('Burn-in? ')))
		steps = str(int(gsTools.check_float('Steps? ')))
		int_points = str(int(gsTools.check_float('Integration points? ')))


		os.system("mkdir " + workdir + '/' + project_name) 
	


		print 'Creating submissions directory'
		os.system("mkdir " + workdir + project_name + '/Submissions')
		fopts = open(workdir + project_name + '/options.txt', 'w')
		fopts.write(plzh + '\t' + anis + '\t' + plummer + '\t' + (steps) + '\t' + (burn_in) + '\t' + (num_walkers))
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

		program = True
		continue

	elif option.strip() == "2":

		

		print "You are about to submit jobs."
		print "If you haven't specified the MPI option for your project previously you can only submit one job at a time to run on your desktop."
		
		quitnow = raw_input("Would you like to continue? y or n? ")
		if quitnow == 'n':
			program = True
			continue
		elif quitnow == 'y':
			pass
		else:
			print "What's that?"
			program = True
			continue
				

		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ").strip()
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()
			
		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			print 'Alright then!'
			
		else:
			print "What was that???"
			program = True
			continue
		
		
		all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(all_gals) == 0:
			print "There are no galaxies in galaxy_list"
			sys.exit()
		gsTools.check_galdata(all_gals,workdir)


		
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
				print "Project description: ", description
			except IOError:
				print "Something weird with project ", pl
			
			print "********************************************"
		#print (contents)
		
		
		
		valid_directory = False
		while valid_directory == False:
			project_name = raw_input("What's the name of the project? " )
			if not os.path.isdir(workdir + project_name + '/Submissions'):
				print "This directory does not exist"
				dec = raw_input("Try again? y or n? ")
				if dec == 'y':
					valid_directory = False
				elif dec == 'n':
					print 'Quitting'
					program = False
					sys.exit()
				else:
					print "What's that?"
					valid_directory = False	
			else:
				valid_directory = True
		valid_dirlist = False	
		while valid_dirlist == False:
			if not os.listdir(workdir + project_name + '/Submissions'):
				print 'Submission directory is empty.'
				print 'Quitting'
				program = False
				sys.exit()	
			else:
				valid_dirlist = True

			
		

		sub_command = open(codedir + '/sub_command.txt', 'r')
		sub_com = sub_command.read()
		sub_com = sub_com.strip()
		if not sub_com:
			print "You have not specified the submission command in the ", codedir,"/sub_command.txt file."
			print "Specify the command now? (e.g. sbatch)"
			sub_com = raw_input("Enter script submission command or leave blank to run a single job on desktop: ") 
				

		print 'Which galaxies would you like to submit?'
		gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(gal_list) == 0:
			print "There are no galaxies, Quitting."
			break
		else:
			print "List of galaxies:"
			for it in gal_list:
				print it
		print '1) All 	2) Specify'
		opt = raw_input('Option: ')
		if opt == '1':
	
			all_gals = np.loadtxt(workdir + '/galaxy_list.txt', ndmin = 1, dtype = 'str')
			if (not sub_com) and (len(all_gals) > 0):
				print "You have not specified the submission command."
				print "If you are not submitting scripts to the system you may only run one galaxy at a time."
				print "I will submit your first galaxy."
				all_gals = [all_gals[0]]
			for gal in all_gals:
				exists = os.path.isdir(workdir  + project_name + '/%s/' %gal)
				if exists == False:
					os.system("mkdir " + workdir  + project_name + '/%s/' %gal)
			if len(all_gals) != 0:
				gsTools.submit_jobs(workdir, project_name, all_gals, sub_com)
			
			program = True
		elif opt == '2':
			print 'Please type galaxies, separated by commas'
			res = raw_input('Type here: ')
			all_gals2 = res.split(',')
			all_gals = []
			for g in all_gals2:
				if g in gal_list:
					all_gals.append(g)
				else:
					print g, "This galaxy does not exist. Did you add it to galaxy_list.txt?"
					valid_galname = False
					while valid_galname == False:
						retype = raw_input("This galaxy does not exist. Retype the name of this galaxy or leave blank: ")
						retype = retype.strip()
						if (retype in gal_list) or retype == "":
							valid_galname = True
						else:
							valid_galname = False
			
			if (not sub_com) and (len(all_gals) > 0):
				print "You have not specified the submission command."
				print "If you are not submitting scripts to the system you may only run one galaxy at a time."
				print "I will submit your first galaxy,", all_gals[0]
				all_gals = [all_gals[0]]		
			
			for gal in all_gals:
				exists = os.path.isdir(workdir  + project_name + '/%s/' %gal)
				if exists == False:
					os.system("mkdir " + workdir  + project_name + '/%s/' %gal)

			if len(all_gals) != 0:
				gsTools.submit_jobs(workdir, project_name, all_gals, sub_com)
			
			program = True
			continue
		else:
			print "What's that??? Try again."
			
			program = True
			continue
			
		
	elif option.strip() == "3":
		
		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ").strip()
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()
			valid = False
		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			print 'Alright then!'
			valid = False
		else:
			print "What was that??? Try again."
			program = True
			continue

		all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(all_gals) == 0:
			print "There are no galaxies in galaxy_list"
			sys.exit()
		gsTools.check_galdata(all_gals,workdir)

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
					print "Project description: ", description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
			
		
			valid_directory = False
			while valid_directory == False:
				project_name = raw_input("What's the name of the project? " )
				if not os.path.isdir(workdir + project_name + '/Submissions'):
					print "This directory does not exist"
					dec = raw_input("Try again? y or n? ")
					if dec == 'y':
						valid_directory = False
					elif dec == 'n':
						print 'Quitting'
						program = False
						sys.exit()
					else:
						print "What's that?"
						valid_directory = False	
				else:
					valid_directory = True
			valid = True
		

		if os.path.isdir(workdir+project_name + '/Analysis') == False:
			os.system("mkdir " + workdir  + project_name + '/Analysis')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Submissions')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Output')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Limits')
			os.system("mkdir " + workdir  + project_name + '/Analysis/OutErr')

		
		foptions = open(workdir + project_name + '/options.txt', 'r')
		input_opt = (foptions.readline()).split()
		dm_option = input_opt[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]
		steps = int(input_opt[3])
		burn_in = int(input_opt[4])
		nwalkers = int(input_opt[5])
		
		valid_options = False
		while valid_options == False:
			print "Do you need to cut chains? "
			print "It looks like you ran %d, excluding burn-in." %(int(steps)-int(burn_in))
			cut_off = gsTools.check_float("How many first steps fould you like to discard for each chain? ")
			if (int(cut_off) >= (int(steps)-int(burn_in))):
				print "Cut off can't be more than total chains. Starting again. \n"
				valid_options = False
				continue
			
			print "Ignore all chains with (chi squared > N x best chi squared)? N = 10 recommended."
			chi_cut = gsTools.check_float("Enter N = ")
			if float(chi_cut) <= 1:
				print "Not enough points selected."
				valid_options = False
				continue
			print "For the radial profile limits, what are the radial ranges in kpc? The intervals will be log-spaced."
			min_rad = gsTools.check_float("Minimum radius = ")
			if float(min_rad) == 0:
				print "Can't do log of zero, setting minimum radius to 0.01"
				min_rad = "0.01"
			max_rad = gsTools.check_float("Maximum radius = ")
			if (float(min_rad) < 0) or (float(max_rad) < 0):
				print "Distances can't be negative. Starting again. \n"
				valid_options = False
				continue
			if (float(max_rad) < float(min_rad)):
				print "Maximum distance cannot be smaller than minimum. Starting again. \n"
				valid_options = False
				continue
			points = gsTools.check_float("How many log-spaced intervals? ")
			if int(points) <= 1:
				print "Need at least two distnaces. Starting again. \n"
				valid_options = False
				continue
			samples = gsTools.check_float("How many samples out of ~ %d that you ran? " %((int(steps)-int(burn_in) - int(cut_off))*int(nwalkers)))
			if int(samples) > (int(steps)-int(burn_in) - int(cut_off))*int(nwalkers):
				print "Can't have more analysis samples than total available. Starting again. \n"
				valid_options = False
				continue
			valid_options = True
		

		sub_scripts = open(codedir + "/sub_script.txt", 'r')
		sub_script = sub_scripts.read()
		sub_script = sub_script.strip()
		sub_scripts.close()
	
		

		valid_opt = False
		while valid_opt == False:
			mpi_opt = raw_input("Running on a batch system y or n? " )
			if mpi_opt == 'y':
				
				if not sub_script:
					print "You do not have anything in the ",codedir+"/sub_script.txt file"
					print "Will create scripts in no batch submission regime."
					mpi_opt = 'n'
					timevar = None
					valid_opt = True 
				else:
					timevar = float(raw_input("How much time do you need? (in hours) "))
					timevar = str(datetime.timedelta(hours = timevar))
					valid_opt = True
			elif mpi_opt == 'n':
				timevar = None	
				valid_opt = True
			else:
				print "What's that?"
				valid_opt = False
			

		print "Generating submission scripts"
		
		gsTools.create_ana_sub(project_name, workdir, codedir, dm_option, beta_option, plummer_option, samples, mpi_opt, cut_off, chi_cut, min_rad, max_rad, points,timevar)		


		
		print "You can submit your analysis jobs at any time by going to project_folder/Analysis/Submissions/"
		opt_sub = raw_input('Would you like to submit a job(s) now? y or n? ')
		if opt_sub.strip() == 'y':

			sub_command = open(codedir + '/sub_command.txt', 'r')
			sub_com = sub_command.read()
			sub_com = sub_com.strip()

			if (mpi_opt == 'y') and (not sub_com) and (sub_script == True):
				print "You have not specified the submission command in the ", codedir,"/sub_command.txt file."
				print "Specify the command now? (e.g. sbatch)"
				sub_com = raw_input("Enter script submission command or leave blank to run on your desktop (single galaxy only): ") 

			



			gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
			if np.size(gal_list) == 0:
				print "There are no galaxies, Quitting."
				break
			else:
				print "Galaxy list:"
				for it in gal_list:
					print it

			opt_gal = raw_input("Which galaxies?  1) Specify    2) All ")
			if opt_gal.strip() == '1':
				print "Type galaxies, separated by commas"
				res = raw_input('Type here: ')
				all_gals2 = res.split(',')
				all_gals = []
				for g in all_gals2:
					if g in gal_list:
						all_gals.append(g)
					else:
						print g, "This galaxy does not exist. Did you add it to galaxy_list.txt?"
						valid_galname = False
						while valid_galname == False:
							retype = raw_input("This galaxy does not exist. Retype the name of this galaxy or leave blank: ")
							retype = retype.strip()
							if (retype in gal_list) or retype == "":
								valid_galname = True
							else:
								valid_galname = False
				
				if (not sub_com) and (len(all_gals) > 0):
					print "You may only submit one galaxy at a time."
					print "I will run your first galaxy, ", all_gals[0]
					all_gals = [all_gals[0]]					
	
				for galaxy in all_gals:
					if not sub_com:
			
						os.system('. ' + workdir + project_name + '/Analysis/Submissions/' + '%s.sh' % galaxy)
					else:
						os.system(sub_com.strip() + ' ' + workdir + project_name + '/Analysis/Submissions/' + '%s.sh' % galaxy)
			elif opt_gal.strip() == '2':
				all_gals = np.loadtxt(workdir + '/galaxy_list.txt', ndmin = 1, dtype = 'str')
				
				if (not sub_com) and (len(all_gals) > 0):
					print "You may only submit one galaxy at a time."
					print "I will run your first galaxy, ", all_gals[0]
					all_gals = [all_gals[0]]
				for galaxy in all_gals:
					if not sub_com:
						os.system('. ' + workdir + project_name + '/Analysis/Submissions/' + '%s.sh' % galaxy)
					else:
						
						os.system(sub_com.strip() + ' ' + workdir + project_name + '/Analysis/Submissions/' + '%s.sh' % galaxy)
			else:
				print "Option doesn't exist. Quitting."	
				sys.exit()
		elif opt_sub.strip() == 'n':
			print 'Oh, I see how it is. Fine then.'	
			program = True
			continue
		else:
			print "Option doesn't exist. Quitting."	
			program = True
			continue


		program = True
		

	elif (option.strip() == '4'):
		all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(all_gals) == 0:
			print "There are no galaxies in galaxy_list"
			sys.exit()
		gsTools.check_galdata(all_gals,workdir)


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
					print "Project description: ", description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
			
			
		
			valid_directory = False
			while valid_directory == False:
				project_name = raw_input("What's the name of the project? " )
				if not os.path.isdir(workdir + project_name + '/Submissions'):
					print "This directory does not exist"
					dec = raw_input("Try again? y or n? ")
					if dec == 'y':
						valid_directory = False
					elif dec == 'n':
						print 'Quitting'
						program = False
						sys.exit()
					else:
						print "What's that?"
						valid_directory = False	
				else:
					valid_directory = True

			valid= True

		print 'Which galaxies would you like to plot?'
		gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(gal_list) == 0:
			print "There are no galaxies, Quitting."
			break
		else:
			print "Galaxy list:"
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
				if g in gal_list:
					all_gals.append(g)
				else:
					print g, "This galaxy does not exist. Did you add it to galaxy_list.txt?"
					valid_galname = False
					while valid_galname == False:
						retype = raw_input("This galaxy does not exist. Retype the name of this galaxy or leave blank: ")
						retype = retype.strip()
						if (retype in gal_list) or retype == "":
							valid_galname = True
						else:
							valid_galname = False
					
		prestr1 = workdir + project_name + '/Submissions/'	
		foptions = open(workdir + project_name + '/options.txt', 'r')
		input_opt = (foptions.readline()).split()
		dm_option = input_opt[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]

		for galaxy in all_gals:
			data = np.loadtxt(workdir + project_name + '/%s/' %galaxy + '%s_Chains' % galaxy + project_name + '.txt')
			gsTools.plot_chains(data, workdir, project_name, galaxy)
			

		program = True
		continue


	elif (option.strip() == '5'):
		all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(all_gals) == 0:
			print "There are no galaxies in galaxy_list"
			sys.exit()
		gsTools.check_galdata(all_gals,workdir)

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
					print "Project description: ", description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
			
			
		
			valid_directory = False
			while valid_directory == False:
				project_name = raw_input("What's the name of the project? " )
				if not os.path.isdir(workdir + project_name + '/Submissions'):
					print "This directory does not exist"
					dec = raw_input("Try again? y or n? ")
					if dec == 'y':
						valid_directory = False
					elif dec == 'n':
						print 'Quitting'
						program = False
						sys.exit()
					else:
						print "What's that?"
						valid_directory = False
				else:
					valid_directory = True
			valid = True

			
		print 'Which galaxies would you like to plot?'
		gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(gal_list) == 0:
			print "There are no galaxies, Quitting."
			break
		else:
			print "Galaxy list:"
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
				if g in gal_list:
					all_gals.append(g)
				else:
					print g, "This galaxy does not exist. Did you add it to galaxy_list.txt?"
					valid_galname = False
					while valid_galname == False:
						retype = raw_input("This galaxy does not exist. Retype the name of this galaxy or leave blank: ")
						retype = retype.strip()
						if (retype in gal_list) or retype == "":
							valid_galname = True
						else:
							valid_galname = False
					
		prestr1 = workdir + project_name + '/Submissions/'	
		foptions = open(workdir + project_name + '/options.txt', 'r')
		input_opt = (foptions.readline()).split()
		dm_option = input_opt[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]
			

		
		for galaxy in all_gals:
			data = np.loadtxt(workdir + project_name + '/%s/' %galaxy + '%s_Chains' % galaxy + project_name + '.txt')
			gsTools.plot_triangle(data, workdir, project_name, galaxy)


		program = True
		continue

	elif (option.strip() == "quit") or (option.strip() == "q"):
		program = False		


	elif option.strip() == '6':
		print "This will generate 0-100 percentile locations for mass, density and anisotropy at N log-spaced distance intervals."
		print "Make sure you run Analysis before using option."
		print 'The current working directory is ' + workdir
		print 'Would you like to change that?  y or n?'
		yon = raw_input("y or n : ")
		if (yon.strip() == 'y') or (yon.strip() == 'Y') or (yon.strip() == 'yes'):
			workdir = raw_input("Enter new name : ").strip()
							
			f = open(workdir + "/workdir.txt", 'w')
			f.write(workdir + '\n')
			f.close()
		elif (yon.strip() == 'n') or (yon.strip() == 'N') or (yon.strip() == 'no'):
			pass
		else:
			print 'Not valid'
			continue
			
		all_gals = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(all_gals) == 0:
			print "There are no galaxies in galaxy_list"
			sys.exit()
		gsTools.check_galdata(all_gals,workdir)
		
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
					print "Project description: ", description
				except IOError:
					print "Something weird with project ", pl
				
				print "********************************************"
			
			valid_directory = False
			while valid_directory == False:
				project_name = raw_input("What's the name of the project? " )
				if not os.path.isdir(workdir + project_name + '/Submissions'):
					print "This directory does not exist"
					dec = raw_input("Try again? y or n? ")
					if dec == 'y':
						valid_directory = False
					elif dec == 'n':
						print 'Quitting'
						program = False
						sys.exit()
					else:
						print "What's that?"
						valid_directory = False
				else:
					valid_directory = True

			valid = True

		
		print 'Which galaxies would you like to plot?'
		gal_list = np.loadtxt(workdir + '/galaxy_list.txt',  ndmin = 1,dtype = 'str')
		if np.size(gal_list) == 0:
			print "There are no galaxies, Quitting."
			break
		else:
			print "Galaxy list:"
			for it in gal_list:
				print it
		
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
				if g in gal_list:
					all_gals.append(g)
				else:
					print g, "This galaxy does not exist. Did you add it to galaxy_list.txt?"
					valid_galname = False
					while valid_galname == False:
						retype = raw_input("This galaxy does not exist. Retype the name of this galaxy or leave blank: ")
						retype = retype.strip()
						if (retype in gal_list) or retype == "":
							valid_galname = True
						else:
							valid_galname = False
					
			all_gals = np.array(all_gals)
		
		prestr1 = workdir + project_name + '/Submissions/'	
			
		min_rad = gsTools.check_float("Minimum radius = ")
		if float(min_rad) == 0:
			print "Can't do log of zero, setting minimum radius to 0.01"
			min_rad = "0.01"
		max_rad = gsTools.check_float("Maximum radius = ")
		if (float(min_rad) < 0) or (float(max_rad) < 0):
			print "Distances can't be negative. Starting again. \n"
			valid_options = False
			continue
		if (float(max_rad) < float(min_rad)):
			print "Maximum distance cannot be smaller than minimum. Starting again. \n"
			valid_options = False
			continue
		points = gsTools.check_float("How many log-spaced intervals? ")
		if int(points) <= 1:
			print "Need at least two distnaces. Starting again. \n"
			valid_options = False
			continue

		r = np.logspace(np.log10(min_rad),np.log10(max_rad), points)
		

		for galaxy in all_gals:
			
	
			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_Mass'%galaxy + '.txt', invalid_raise = False)
				res = gsTools.get_lims_all(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Mass'% galaxy + 'LimsAll.txt', res)
			except IOError:
				print 'Mass output not found'

			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_Density'%galaxy + '.txt', invalid_raise = False)
				res = gsTools.get_lims_all(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Density'% galaxy + 'LimsAll.txt', res)
			except IOError:
				print 'Density output not found'

			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_Beta'%galaxy + '.txt', invalid_raise = False)
				res = gsTools.get_lims_all(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Beta'% galaxy + 'LimsAll.txt', res)
			except IOError:
				print 'Beta output not found'

			print 'Galaxy', galaxy, 'Done'

		program = True
		continue			
