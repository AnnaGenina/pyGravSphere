import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import corner
import glob

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


codedir = str(os.environ["GravSpherePath"])
	
sys.path.append(codedir)
from GStools import gsTools


gsTools.checkdirs(workdir, codedir)



program = True

while program == True:

	banner = open(codedir + '/banner.txt', 'r')
	banner = banner.read()
	print (banner)

	print "Hello galactic dynamicist, what would you like to do?"
	print "0) Preprocess data"
	print "1) Create a new project" 
	print "2) Submit jobs"
	print "3) Analysis"
	print "4) Get limits"
	print "5) Plot chain convergence"
	print "6) Create a corner plot"
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
		mpi_opt = raw_input('Using multiple cores? y or n?')
		if mpi_opt == 'y':
			num_cores = int(raw_input('How many cores? '))
			timevar = int(raw_input('How much time do you need (hours only)? '))
			


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
		fopts.write(plzh + '\t' + anis + '\t' + plummer)
		fopts.close()



		print 'Copying codes'
		#os.system("cp " + workdir + '/gravsphere.py ' + workdir + project_name + '/Submissions/')
		#os.system("cp " + workdir + '/_gravsphere.so ' + workdir + project_name + '/Submissions/')
		os.system("cp " + codedir + '/priors.txt ' + workdir + project_name + '/Submissions/')
		#os.system("cp " + workdir + '/write_script.py ' + workdir + project_name + '/Submissions/')
		priors = np.loadtxt(workdir + project_name + '/Submissions/priors.txt', dtype = 'str')
		my_params = np.where(np.in1d(priors[:,0], np.array(param_list, dtype = 'str')))
		print priors[:,0], np.array(param_list, dtype = 'str'), my_params
		priors = np.savetxt(workdir + project_name + '/Submissions/priors.txt', priors[my_params], fmt = "%s")

		print 'Creating OutErr directory'
		os.system("mkdir " + workdir  + project_name + '/OutErr/')
		print 'Creating Triangle directory'
		os.system("mkdir " + workdir  + project_name + '/Triangle/')


		
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

				
				
				sub_command = open('sub_command.txt', 'r')
				sub_com = sub_command.read()
						

				print 'Which galaxies would you like to submit?'
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

		

		if os.path.isdir(workdir+project_name + '/Analysis') == False:
			os.system("mkdir " + workdir  + project_name + '/Analysis')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Submissions')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Output')
			os.system("mkdir " + workdir  + project_name + '/Analysis/Limits')


		samples = raw_input("How many samples?" )

		print "Generating submission scripts"
			

		foptions = open(workdir + project_name + '/options.txt', 'r')
		input_opt = line.split()
		dm_option = foptions[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]
		
		
		gsTools.create_ana_sub(project_name, workdir, codedir, dm_option, beta_option, plummer_option, samples)		

		opt_sub = raw_input('Would you like to submit a job? y or n?')
		if opt_sub.strip() == 'y':

			sub_command = open('sub_command.txt', 'r')
			sub_com = sub_command.read()

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
			break
			

		
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
		input_opt = line.split()
		dm_option = foptions[0] 
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

			
			if dm_option == 'PL':
				try:
					f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_LogLog'%galaxy + '.txt', invalid_raise = False)
					res = gsTools.get_lims_loglog(f, tot_bins)
					np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_LogLog'% galaxy + 'Lims.txt', res)
				except IOError:
					print 'Slope output not found'
			else:
				try:
					f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_LogLog'%galaxy + '.txt' , invalid_raise = False)
					res = gsTools.get_lims(f,r)
					np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_LogLog'% galaxy + 'Lims.txt', res)
				except IOError:
					print 'Slope output not found'
			print 'Galaxy', galaxy, 'Slopes done'
			
			try:	
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_Beta'%galaxy + '.txt', invalid_raise = False)
				res = gsTools.get_lims(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Beta'% galaxy +  'Lims.txt', res)
			except IOError:
				print 'Beta output not found'
			print 'Galaxy', galaxy, 'Beta done'

			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_Density'%galaxy  +'.txt' , invalid_raise = False)
					
				res = gsTools.get_lims(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/' + '%s_Density'% galaxy +  'Lims.txt', res)
			except IOError:
				print 'Density output not found'
			print 'Galaxy', galaxy, 'Density done'

			
			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_Mass'%galaxy + '.txt', invalid_raise = False)
				res = gsTools.get_lims(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Mass'% galaxy + 'Lims.txt', res)
			except IOError:
				print 'Mass output not found'

			print 'Galaxy', galaxy, 'Mass done'

			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_Plummer'%galaxy + '.txt', invalid_raise = False)
				res = gsTools.get_lims(f,r)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_Plummer'% galaxy + 'Lims.txt', res)
			except IOError:
				print 'Plummer output not found'

			print 'Galaxy', galaxy, 'Plummer done'	

			try:
				f = np.genfromtxt(workdir + project_name + '/Analysis/Output/' + '%s_SigLos'%galaxy + '.txt', invalid_raise = False)
				res = gsTools.get_lims_sig(f, workdir, galaxy)
				np.savetxt(workdir + project_name + '/Analysis/Limits/'+ '%s_SigLos'% galaxy + 'Lims.txt', res)
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

		print 'Which galaxies would you like to plot?'
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
		input_opt = line.split()
		dm_option = foptions[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]

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
			data = np.loadtxt(workdir + project_name + '/%s/' %galaxy + '%s_Chains' % galaxy + project_name + '.txt')
			gsTools.plot_chains(data, param_names, nparams, workdir, project_name, galaxy)
			




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

		print 'Which galaxies would you like to plot?'
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
		input_opt = line.split()
		dm_option = foptions[0] 
		beta_option = input_opt[1]
		plummer_option = input_opt[2]	

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
			data = np.loadtxt(workdir + project_name + '/%s/' %galaxy + '%s_Chains' % galaxy + project_name + '.txt')
			gsTools.plot_triangle(data, param_names, nparams, workdir, project_name, galaxy)


	

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
		input_opt = line.split()
		dm_option = foptions[0] 
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

			
