import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import h5py
from scipy.integrate import simps

def gs_cofficient(v1, v2):
    return np.dot(v2, v1) / np.dot(v1, v1)

def multiply(cofficient, v):
    return map((lambda x : x * cofficient), v)

def proj(v1, v2):
    return multiply(gs_cofficient(v1, v2) , v1)

def gs(X):
    Y = []
    for i in range(len(X)):
        temp_vec = X[i]
        for inY in Y :
            proj_vec = proj(inY, X[i])
            temp_vec = map(lambda x, y : x - y, temp_vec, proj_vec)
        Y.append(temp_vec)
    return Y



def plummer(x,norm,rc):
	return np.log10((4./3.)  * rc / (1 + x*x/(rc*rc))**2.) + norm

def min_plummer(args,x_data, y_data,err_data):
	norm,rc = args
	return np.sum(((plummer(x_data,norm,rc)) - np.log10(y_data))**2./err_data**2.)

def zhao(x,rhos,rs,alpha,beta,gamma):
	return (rhos) + np.log10(1./  ( (x/rs)**gamma * (1. + (x/rs)**alpha)**((beta - gamma)/alpha)   )  )

def min_zhao(args, x_data, y_data,err_data):
	rhos,rs,alpha,beta,gamma = args
	return np.sum(((zhao(x_data,rhos,10**rs,alpha,beta,gamma)) - np.log10(y_data))**2./err_data**2.)


def plummer_dens_comp(x,mass,shape):
	return 3. * (10**mass)/(4. * np.pi * (10**shape)**3.) * (1. + x**2./(10**shape)**2.)**(-5./2.)

def plummer_proj_comp(x,mass,shape):
	#return 15. * (10**mass)/(4. * np.pi * shape**5.) * (1. + x**2./shape**2.)**(-7./2.)
	return 10**mass/(np.pi * (10**shape)**2.) * (1. + x**2./(10**shape)**2)**(-2.)
def plummer_proj_sum(args, x_data, n_comp):

	sums_proj = 0
	#sums_3d = 0
	start = 0
	for i in range(0, n_comp):
		sums_proj = sums_proj + plummer_proj_comp(x_data,args[start], args[start + 1])
		#sums_3d = sums_3d + plummer_dens_comp(x_data,args[start], args[start + 1])
		start = start + 2


	
	return np.log10(sums_proj)

#R = np.logspace(-2,2, 100)
#comps = plummer_proj_sum([3,1,2,1.5, 1,2], R,3)
#plt.plot(np.log10(R), np.log10(comps[0]), label = 'analytic')
#res = []
#for val in R:
#	print val
#	x1 = np.logspace(np.log10(val),3,20000)
#	xint = (x1[1:] + x1[:-1])/2.
#	nuint = plummer_proj_sum([3,1,2,1.5, 1,2], xint,3)[1]
#	integrand = xint * nuint  * 2 / np.sqrt(xint**2. - val**2.)
#	ss = simps(integrand,xint,even = 'avg')
#	res.append(np.log10(ss))
#res = np.array(res)
#plt.plot(np.log10(R), res, label = 'projection')
#plt.show()


def min_plummer_sum(args,x_data,y_data,err_data,n_comp):
	return np.sum((plummer_proj_sum(args,x_data,n_comp) - np.log10(y_data))**2./err_data**2.)


def like_func(args, vel, vel_err):
	
	#print vel

#	return -np.sum( (-0.5 * np.log(2. * np.pi *(vel_err**2. + args[1]**2)) - 0.5 * (vel - args[0])**2./(vel_err**2. + args[1]**2.)  )   ) 
	
	return -np.sum( np.log(  (1./np.sqrt(vel_err**2. + args[1]**2))) +( -0.5 * (vel - args[0])**2./(vel_err**2. + args[1]**2.)   )   )    

	

def like_func_fixed(args,mean_vel, vel, vel_err):
	
	return -np.sum( np.log(  (1./np.sqrt(vel_err**2. + args[0]**2))) + ( -0.5 * (vel - mean_vel)**2./(vel_err**2. + args[0]**2.)   )   )    



def cov_sig(mean_vel,vel_disp,vel,vel_err):
	d = np.sum(0.5 * ( (4. * vel_disp**2)/(vel_err**2. + vel_disp**2.)**2.  - 2./(vel_err**2. + vel_disp**2.) ) - 0.5 *  (vel -mean_vel )**2. * ((8*vel_disp**2.)/(vel_err*2. + vel_disp**2.)**3. - 2./(vel_disp**2. + vel_err**2.)**2.   ))
	a = np.sum(-(1)/(vel_err*2. + vel_disp**2.))
	c = np.sum(2 * vel_disp/(vel_err*2. + vel_disp**2.)**2.)
	return a/(a*d - c*c)


def projector(volume, subhalo, observer):

	mains = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + str(volume) + 'FinalID.txt')
	main = None
	for m in mains:
		subhalos = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_Main_%d_Subhalos.txt' %m)
		if int(subhalo) in subhalos:
			main = int(m)

	com_sub = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_SubhaloCOM.txt')

	cc, = np.where(com_sub[:,0] == int(subhalo))
	centre_subhalo = com_sub[cc,1:4]

	centre_subhalo = centre_subhalo.reshape((1,3))



	if observer == 0:
		los = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_Main_' + str(int(main)) + '_COM.txt')
		observer = los.reshape((1,3))

		axis0 = observer - centre_subhalo
		
	

	else:
		los = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/KinDat/' + volume + '_' + str(int(subhalo)) + '_Axes.txt')
		observer = los[observer - 1]
		observer = observer.reshape((1,3))
		axis0 = observer 

	kindat = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/KinDat/' + volume + '_' + str(int(subhalo)) + '_KinDatFull.txt')

	
	center = axis0.reshape((3,))
	norm = center/np.sqrt(center[0]**2. + center[1]**2. + center[2]**2.)





	pos = kindat[:,0:3]
	vel = kindat[:,3:6]
	mass = kindat[:,6]

	

	com_sub = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_SubhaloCOM.txt')

	
        		
	cc, = np.where(com_sub[:,0] == int(subhalo))
	centre_subhalo = com_sub[cc,1:4]

	centre_subhalo = centre_subhalo.reshape((1,3))

	axis0 = observer - centre_subhalo
	center = axis0.reshape((3,))
	norm = center/np.sqrt(center[0]**2. + center[1]**2. + center[2]**2.)
			
	if norm[2] != 0:
		test = np.array([norm,[0., 1., 0.], [1., 0., 0.]])
	else:
		test = np.array([norm,[0., 1., 0.], [0., 0., 1.]])

	obs = np.array(gs(test))


	axis1 = obs[1,:].reshape((1,3))
	axis2 = obs[2,:].reshape((1,3))
	axis1 = axis1/ np.sqrt(axis1[0,0]**2. + axis1[0,1]**2. + axis1[0,2]**2.)
	axis2 = axis2/ np.sqrt(axis2[0,0]**2. + axis2[0,1]**2. + axis2[0,2]**2.)

	axis0 = obs[0,:].reshape((1,3))
	

	axis0 = axis0/ np.sqrt(axis0[0,0]**2. + axis0[0,1]**2. + axis0[0,2]**2.)
				
	pos_al_0 = np.zeros((len(pos), 2))
	

	pos_al_0[:,0] = (axis1 * pos ) .sum(1)
	pos_al_0[:,1] = (axis2 * pos ) .sum(1)
		
	vel_al_0 = ((vel * axis0).sum(1))
	
	vel_al_0 = vel_al_0 + np.random.normal(0,2, vel_al_0.shape)


	glob = minimize(like_func,[np.mean(vel_al_0),np.std(vel_al_0)],args = (vel_al_0,2.))
	opt_u = glob.x[0]
	vel_al_0 = vel_al_0 - opt_u



	pos2d = np.sqrt(np.sum(pos_al_0**2., axis = 1))

	
	rh = np.median(pos2d)

#	return (pos_al_0, vel_al_0)
	return (pos2d, vel_al_0, rh)



def bin_projector(volume, subhalo, observer):

	mains = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + str(volume) + 'FinalID.txt')
	main = None
	for m in mains:
		subhalos = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_Main_%d_Subhalos.txt' %m)
		if int(subhalo) in subhalos:
			main = int(m)

	los = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_' + str(int(main)) + '_10los.txt')
	observer = los[observer]
	observer = observer.reshape((1,3))

	kindat = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/KinDat/' + volume + '_' + str(int(subhalo)) + '_KinDatFull.txt')

	

	pos = kindat[:,0:3]
	vel = kindat[:,3:6]
	mass = kindat[:,6]

	

	com_sub = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_SubhaloCOM.txt')

	
        		
	cc, = np.where(com_sub[:,0] == int(subhalo))
	centre_subhalo = com_sub[cc,1:4]

	centre_subhalo = centre_subhalo.reshape((1,3))

	axis0 = observer - centre_subhalo
	center = axis0.reshape((3,))
	norm = center/np.sqrt(center[0]**2. + center[1]**2. + center[2]**2.)
			
	if norm[2] != 0:
		test = np.array([norm,[0., 1., 0.], [1., 0., 0.]])
	else:
		test = np.array([norm,[0., 1., 0.], [0., 0., 1.]])

	obs = np.array(gs(test))


	axis1 = obs[1,:].reshape((1,3))
	axis2 = obs[2,:].reshape((1,3))
	axis1 = axis1/ np.sqrt(axis1[0,0]**2. + axis1[0,1]**2. + axis1[0,2]**2.)
	axis2 = axis2/ np.sqrt(axis2[0,0]**2. + axis2[0,1]**2. + axis2[0,2]**2.)

	axis0 = obs[0,:].reshape((1,3))
	

	axis0 = axis0/ np.sqrt(axis0[0,0]**2. + axis0[0,1]**2. + axis0[0,2]**2.)
				
	pos_al_0 = np.zeros((len(pos), 2))
	

	pos_al_0[:,0] = (axis1 * pos ) .sum(1)
	pos_al_0[:,1] = (axis2 * pos ) .sum(1)
		
	vel_al_0 = ((vel * axis0).sum(1))
	vel_al_0 = vel_al_0 + np.random.normal(0,2, vel_al_0.shape)

	pos2d = np.sqrt(np.sum(pos_al_0**2., axis = 1))

	
	rh = np.median(pos2d)


	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))

	glob = minimize(like_func,[np.mean(vel_al_0),np.std(vel_al_0)],args = (vel_al_0,2.,))
	opt_u = glob.x[0]
	
	
	

	opt_disp = glob.x[1]


	veldisps = []
	veldisperrs = []

	poss = []
	for bins in binned:
		
		bin_l = minimize(like_func_fixed, [np.std(vel_al_0[bins])],args = (opt_u, vel_al_0[bins],2.,)  )
		
		vd = float(bin_l.x)
		vde = float( np.abs(cov_sig(opt_u,vd,vel_al_0[bins],2.) )) 
		veldisps.append(vd)
		veldisperrs.append( np.sqrt(vde))
		bind = float(np.median(pos2d[bins]))
		poss.append(bind)

		


#	return (pos_al_0, vel_al_0)

	poss = np.array(poss)
	veldisps = np.array(veldisps)
	veldisperrs = np.array(veldisperrs)	


	data_x = []
	data_y = []
	
	for bins in binned:
		bind = float(np.median(pos2d[bins]))
		data_x.append(float(bind))
		area = np.pi *(np.amax(pos2d[bins])**2. - np.amin(pos2d[bins])**2.)
		data_y.append(float(len(bins)) / area)

	data_x = np.array(data_x)
	data_y = np.array(data_y)

	data_y = data_y / np.amax(data_y) 

	print 'Old rh', rh
	res = minimize(min_plummer, [rh], args=(data_x, data_y), bounds = ((0.1,3),))
	rh = float(res.x)
	print 'New rh', rh


#	plt.plot(np.log10(data_x), np.log10(data_y),color = 'black')
#	plt.plot(np.log10(data_x), np.log10(plummer(data_x,rh)), color = 'blue')
#	plt.show()

#	plt.errorbar(poss, veldisps, yerr = veldisperrs)
#	plt.show()

	return (poss, veldisps, veldisperrs, rh)


#bin_projector('S5', 5,0)



def bin_projector_power(volume, subhalo, observer):

	vmax = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_Vmax.txt') 	

	v, = np.where(vmax[:,0] == int(subhalo))
	
	power = vmax[v,4] 


	mains = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + str(volume) + 'FinalID.txt')
	main = None
	for m in mains:
		subhalos = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_Main_%d_Subhalos.txt' %m)
		if int(subhalo) in subhalos:
			main = int(m)

	com_sub = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_SubhaloCOM.txt')

	cc, = np.where(com_sub[:,0] == int(subhalo))
	centre_subhalo = com_sub[cc,1:4]

	centre_subhalo = centre_subhalo.reshape((1,3))

	print observer

	if observer == 0:
		los = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/' + volume + '_Main_' + str(int(main)) + '_COM.txt')
		observer = los.reshape((1,3))
		
		axis0 = observer - centre_subhalo
		
	

	else:
		los = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/KinDat/' + volume + '_' + str(int(subhalo)) + '_Axes.txt')
		
		observer = los[observer - 1]
		
		observer = observer.reshape((1,3))
		axis0 = observer 

	kindat = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/KinDat/' + volume + '_' + str(int(subhalo)) + '_KinDatFull.txt')

	
	center = axis0.reshape((3,))
	norm = center/np.sqrt(center[0]**2. + center[1]**2. + center[2]**2.)

	print norm



	pos = kindat[:,0:3]
	vel = kindat[:,3:6]
	mass = kindat[:,6]


	R3d = np.sqrt(np.sum(pos**2., axis = 1))	


	

	
			
	if norm[2] != 0:
		test = np.array([norm,[0., 1., 0.], [1., 0., 0.]])
	else:
		test = np.array([norm,[0., 1., 0.], [0., 0., 1.]])

	obs = np.array(gs(test))


	axis1 = obs[1,:].reshape((1,3))
	axis2 = obs[2,:].reshape((1,3))
	axis1 = axis1/ np.sqrt(axis1[0,0]**2. + axis1[0,1]**2. + axis1[0,2]**2.)
	axis2 = axis2/ np.sqrt(axis2[0,0]**2. + axis2[0,1]**2. + axis2[0,2]**2.)

	axis0 = obs[0,:].reshape((1,3))
	

	axis0 = axis0/ np.sqrt(axis0[0,0]**2. + axis0[0,1]**2. + axis0[0,2]**2.)
				
	pos_al_0 = np.zeros((len(pos), 2))
	

	pos_al_0[:,0] = (axis1 * pos ) .sum(1)
	pos_al_0[:,1] = (axis2 * pos ) .sum(1)
		
	vel_al_0 = ((vel * axis0).sum(1))
	vel_al_0 = vel_al_0 + np.random.normal(0,2, vel_al_0.shape)

	pos2d = np.sqrt(np.sum(pos_al_0**2., axis = 1))

	
	rh = np.median(pos2d)


	below, = np.where(pos2d > 2.8*0.134)
	pos2d = pos2d[below]
	vel_al_0 = vel_al_0[below]


	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))

	data_x = []
	data_y = []
	
	for bins in binned:
		bind = float(np.median(pos2d[bins]))
		data_x.append(float(bind))
		area = np.pi *(np.amax(pos2d[bins])**2. - np.amin(pos2d[bins])**2.)
		data_y.append(float(len(bins)) / area)

	data_x = np.array(data_x)
	data_y = np.array(data_y)


	


	data_y = data_y 
	normal = np.log10(np.amax(data_y))

	print 'Old rh', rh
	res = minimize(min_plummer, [normal,rh], args=(data_x, data_y), bounds = ((0, normal + 5),(0.1,3),))
	normal,rh = res.x
	print 'New rh', rh



	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))

	glob = minimize(like_func,[np.mean(vel_al_0),np.std(vel_al_0)],args = (vel_al_0,2.,))
	opt_u = glob.x[0]
	
	
	

	opt_disp = glob.x[1]


	veldisps = []
	veldisperrs = []

	poss = []
	for bins in binned:
		
		bin_l = minimize(like_func_fixed, [np.std(vel_al_0[bins])],args = (opt_u, vel_al_0[bins],2.,)  )
		
		vd = float(bin_l.x)
		vde = float( np.abs(cov_sig(opt_u,vd,vel_al_0[bins],2.) )) 
		veldisps.append(vd)
		veldisperrs.append( np.sqrt(vde))
		bind = float(np.median(pos2d[bins]))
		poss.append(bind)

		




	poss = np.array(poss)
	veldisps = np.array(veldisps)
	veldisperrs = np.array(veldisperrs)	



	return (poss, veldisps, veldisperrs, rh)



def bin_projector_read(galaxy_number):

	f = h5py.File('/cosma5/data/dp004/dc-geni1/project/JR_mocks/OverallSample/Galaxy_%d' %galaxy_number + '.hdf5' , 'r')	

	
	pos = f['KinematicsPositions'].value
	pos = pos/1000.
	
	pos = np.sqrt(np.sum(pos**2,axis = 1))

	vel = f['KinematicsVelocities'].value
	mass = f['KinematicsMasses'].value

	pos_p = f['PhotometryPositions'].value
	pos_p = pos_p/1000.
	pos_p = np.sqrt(np.sum(pos_p**2,axis = 1))
	mass_p = f['PhotometryMasses'].value


	rh = np.median(pos_p)


	
	pos2d = pos_p  # photometry atm
	vel_al_0 = vel # kinematics atm


	print int(np.sqrt(len(pos2d))), 'Photometry bins', 2*int(np.sqrt(len(pos2d))), 3*int(np.sqrt(len(pos2d)))

	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))
	#binned = np.array_split(sorted_r, 279)
	
	data_x = []
	data_y = []
	data_err = []
	for bins in binned:
		bind = float(np.median(pos2d[bins]))
		data_x.append(float(bind))
		area = np.pi *(np.amax(pos2d[bins])**2. - np.amin(pos2d[bins])**2.)
		data_y.append(float(len(bins)) / area)
		data_err.append(float(np.sqrt(len(bins)))/area)
		#data_y.append(np.sum(mass_p)/area)
	data_x = np.array(data_x)
	data_y = np.array(data_y)
	#data_y2 = np.array(data_y2)
	data_err = np.array(data_err)

	data_err_log = np.zeros_like(data_err)

	for e in range(0, len(data_err)):
		distr = np.random.normal(data_y[e], data_err[e], 10000)
		data_err_log[e] = np.std(np.log10(distr))



	#data_y2 = data_y2
	normal = np.log10(np.amax(data_y))
	print 'Normal', normal
	#print 'Old rh', rh
	


	

	

#	res3 = minimize(min_plummer_sum, [normal, np.log10(rh/2), normal, np.log10(rh), normal, np.log10(rh*2)], args=(data_x, data_y,data_err_log, 3), bounds = ((0, 10),(-2,1),(0, 10),(-2,1),(0, 10),(-2,1),))
#	n1, a1, n2,a2,n3,a3 = res3.x
#	light_params = np.array([n1,a1,n2,a2,n3,a3])
#	np.savetxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/Galaxy_%d_PlumParam.txt' %galaxy_number, light_params)


	plum_param = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/Galaxy_%d_PlumParam.txt' %galaxy_number)

	n1,a1,n2,a2,n3,a3 = plum_param

	chi_plummer_sum = min_plummer_sum([n1,a1,n2,a2,n3,a3], data_x,data_y, data_err_log, 3)

	print 'Plummer sphere parameters', n1,a1,n2,a2,n3,a3

	print 'Chisq', chi_plummer_sum

	

	
	

	


	pos2d = pos   # Now kinematics

	print int(np.sqrt(len(pos2d))), 'Kinematics bins'

	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))

	glob = minimize(like_func,[np.mean(vel_al_0),np.std(vel_al_0)],args = (vel_al_0,2.,))
	opt_u = glob.x[0]
	
	
	

	opt_disp = glob.x[1]


	veldisps = []
	veldisperrs = []

	poss = []
	for bins in binned:
		
		bin_l = minimize(like_func_fixed, [np.std(vel_al_0[bins])],args = (opt_u, vel_al_0[bins],2.,)  )
		
		vd = float(bin_l.x)
		vde = float( np.abs(cov_sig(opt_u,vd,vel_al_0[bins],2.) )) 
		veldisps.append(vd)
		veldisperrs.append( np.sqrt(vde))
		bind = float(np.median(pos2d[bins]))
		poss.append(bind)

		

	poss = np.array(poss)
	veldisps = np.array(veldisps)
	veldisperrs = np.array(veldisperrs)	

	

	light_params = np.array([n1,a1,n2,a2,n3,a3])



	#return (poss, veldisps, veldisperrs, rh)
	return (poss,veldisps,veldisperrs,light_params)


def bin_projector_read_vsp(galaxy_number):

	f = h5py.File('/cosma5/data/dp004/dc-geni1/project/JR_mocks/OverallSample/Galaxy_%d' %galaxy_number + '.hdf5' , 'r')	

	
	pos = f['KinematicsPositions'].value
	pos = pos/1000.
	
	pos = np.sqrt(np.sum(pos**2,axis = 1))

	vel = f['KinematicsVelocities'].value
	mass = f['KinematicsMasses'].value

	pos_p = f['PhotometryPositions'].value
	pos_p = pos_p/1000.
	pos_p = np.sqrt(np.sum(pos_p**2,axis = 1))
	mass_p = f['PhotometryMasses'].value


	rh = np.median(pos_p)


	
	pos2d = pos_p  # photometry atm
	vel_al_0 = vel # kinematics atm


	print int(np.sqrt(len(pos2d))), 'Photometry bins', 2*int(np.sqrt(len(pos2d))), 3*int(np.sqrt(len(pos2d)))

	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))
	#binned = np.array_split(sorted_r, 279)
	
	data_x = []
	data_y = []
	data_err = []
	for bins in binned:
		bind = float(np.median(pos2d[bins]))
		data_x.append(float(bind))
		area = np.pi *(np.amax(pos2d[bins])**2. - np.amin(pos2d[bins])**2.)
		data_y.append(float(len(bins)) / area)
		data_err.append(float(np.sqrt(len(bins)))/area)
		#data_y.append(np.sum(mass_p)/area)
	data_x = np.array(data_x)
	data_y = np.array(data_y)
	#data_y2 = np.array(data_y2)
	data_err = np.array(data_err)

	data_err_log = np.zeros_like(data_err)

	for e in range(0, len(data_err)):
		distr = np.random.normal(data_y[e], data_err[e], 10000)
		data_err_log[e] = np.std(np.log10(distr))



	#data_y2 = data_y2
	normal = np.log10(np.amax(data_y))
	print 'Normal', normal
	#print 'Old rh', rh
	


	

	

#	res3 = minimize(min_plummer_sum, [normal, np.log10(rh/2), normal, np.log10(rh), normal, np.log10(rh*2)], args=(data_x, data_y,data_err_log, 3), bounds = ((0, 10),(-2,1),(0, 10),(-2,1),(0, 10),(-2,1),))
#	n1, a1, n2,a2,n3,a3 = res3.x
#	light_params = np.array([n1,a1,n2,a2,n3,a3])
#	np.savetxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/Galaxy_%d_PlumParam.txt' %galaxy_number, light_params)


	plum_param = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/Galaxy_%d_PlumParam.txt' %galaxy_number)
	vsp_param = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/Galaxy_%d_VSP.txt' %galaxy_number)


	n1,a1,n2,a2,n3,a3 = plum_param

	chi_plummer_sum = min_plummer_sum([n1,a1,n2,a2,n3,a3], data_x,data_y, data_err_log, 3)

	print 'Plummer sphere parameters', n1,a1,n2,a2,n3,a3

	print 'Chisq', chi_plummer_sum

	

	
	

	


	pos2d = pos   # Now kinematics

	print int(np.sqrt(len(pos2d))), 'Kinematics bins'

	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))

	glob = minimize(like_func,[np.mean(vel_al_0),np.std(vel_al_0)],args = (vel_al_0,2.,))
	opt_u = glob.x[0]
	
	
	

	opt_disp = glob.x[1]


	veldisps = []
	veldisperrs = []

	poss = []
	for bins in binned:
		
		bin_l = minimize(like_func_fixed, [np.std(vel_al_0[bins])],args = (opt_u, vel_al_0[bins],2.,)  )
		
		vd = float(bin_l.x)
		vde = float( np.abs(cov_sig(opt_u,vd,vel_al_0[bins],2.) )) 
		veldisps.append(vd)
		veldisperrs.append( np.sqrt(vde))
		bind = float(np.median(pos2d[bins]))
		poss.append(bind)

		

	poss = np.array(poss)
	veldisps = np.array(veldisps)
	veldisperrs = np.array(veldisperrs)	

	

	light_params = np.array([n1,a1,n2,a2,n3,a3])



	#return (poss, veldisps, veldisperrs, rh)
	return (poss,veldisps,veldisperrs,light_params,vsp_param)



def bin_projector_read_vsp_bin(galaxy_number):

	f = h5py.File('/cosma5/data/dp004/dc-geni1/project/JR_mocks/OverallSample/Galaxy_%d' %galaxy_number + '.hdf5' , 'r')	

	
	pos = f['KinematicsPositions'].value
	pos = pos/1000.
	
	pos = np.sqrt(np.sum(pos**2,axis = 1))

	vel = f['KinematicsVelocities'].value
	mass = f['KinematicsMasses'].value

	pos_p = f['PhotometryPositions'].value
	pos_p = pos_p/1000.
	pos_p = np.sqrt(np.sum(pos_p**2,axis = 1))
	mass_p = f['PhotometryMasses'].value


	rh = np.median(pos_p)


	
	pos2d = pos_p  # photometry atm
	vel_al_0 = vel # kinematics atm


	print int(np.sqrt(len(pos2d))), 'Photometry bins', 2*int(np.sqrt(len(pos2d))), 3*int(np.sqrt(len(pos2d)))

	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))
	#binned = np.array_split(sorted_r, 279)
	
	data_x = []
	data_y = []
	data_err = []
	for bins in binned:
		bind = float(np.median(pos2d[bins]))
		data_x.append(float(bind))
		area = np.pi *(np.amax(pos2d[bins])**2. - np.amin(pos2d[bins])**2.)
		data_y.append(float(len(bins)) / area)
		data_err.append(float(np.sqrt(len(bins)))/area)
		#data_y.append(np.sum(mass_p)/area)
	data_x = np.array(data_x)
	data_y = np.array(data_y)
	#data_y2 = np.array(data_y2)
	data_err = np.array(data_err)

	data_err_log = np.zeros_like(data_err)

	for e in range(0, len(data_err)):
		distr = np.random.normal(data_y[e], data_err[e], 10000)
		data_err_log[e] = np.std(np.log10(distr))



	#data_y2 = data_y2
	normal = np.log10(np.amax(data_y))
	print 'Normal', normal
	#print 'Old rh', rh
	


	

	

#	res3 = minimize(min_plummer_sum, [normal, np.log10(rh/2), normal, np.log10(rh), normal, np.log10(rh*2)], args=(data_x, data_y,data_err_log, 3), bounds = ((0, 10),(-2,1),(0, 10),(-2,1),(0, 10),(-2,1),))
#	n1, a1, n2,a2,n3,a3 = res3.x
#	light_params = np.array([n1,a1,n2,a2,n3,a3])
#	np.savetxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/Galaxy_%d_PlumParam.txt' %galaxy_number, light_params)


	plum_param = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/Galaxy_%d_PlumParam.txt' %galaxy_number)
	vsp_param = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/Galaxy_%d_VSP.txt' %galaxy_number)


	n1,a1,n2,a2,n3,a3 = plum_param

	chi_plummer_sum = min_plummer_sum([n1,a1,n2,a2,n3,a3], data_x,data_y, data_err_log, 3)

	print 'Plummer sphere parameters', n1,a1,n2,a2,n3,a3

	print 'Chisq', chi_plummer_sum

	

	
	

	


	pos2d = pos   # Now kinematics

	print int(np.sqrt(len(pos2d))), 'Kinematics bins'

	sorted_r = np.argsort(pos2d)
	binned = np.array_split(sorted_r, int(np.sqrt(len(pos2d))))
	#binned = np.array_split(sorted_r, 52)


	glob = minimize(like_func,[np.mean(vel_al_0),np.std(vel_al_0)],args = (vel_al_0,2.,))
	opt_u = glob.x[0]
	
	
	

	opt_disp = glob.x[1]


	veldisps = []
	veldisperrs = []

	poss = []
	for bins in binned:
		
		bin_l = minimize(like_func_fixed, [np.std(vel_al_0[bins])],args = (opt_u, vel_al_0[bins],2.,)  )
		
		vd = float(bin_l.x)
		vde = float( np.abs(cov_sig(opt_u,vd,vel_al_0[bins],2.) )) 
		veldisps.append(vd)
		veldisperrs.append( np.sqrt(vde))
		bind = float(np.median(pos2d[bins]))
		poss.append(bind)

		

	poss = np.array(poss)
	veldisps = np.array(veldisps)
	veldisperrs = np.array(veldisperrs)	

	

	light_params = np.array([n1,a1,n2,a2,n3,a3])



	#return (poss, veldisps, veldisperrs, rh)
	return (poss,veldisps,veldisperrs,light_params,vsp_param, rh)

def bin_projector_read_vsp_bin_copy(galaxy_number):

	f = h5py.File('/cosma5/data/dp004/dc-geni1/project/JR_mocks/OverallSample/Galaxy_%d' %galaxy_number + '.hdf5' , 'r')	

	
	pos = f['KinematicsPositions'].value
	pos = pos/1000.
	
	pos = np.sqrt(np.sum(pos**2,axis = 1))

	vel = f['KinematicsVelocities'].value
	mass = f['KinematicsMasses'].value

	pos_p = f['PhotometryPositions'].value
	pos_p = pos_p/1000.
	pos_p = np.sqrt(np.sum(pos_p**2,axis = 1))
	mass_p = f['PhotometryMasses'].value
	mass_tot = f['StellarMass3R'].value

	mass_tot = float(mass_tot)

	#rh = np.median(pos_p)  Replace with same as Justin


	


	light_params = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/ReadData/Galaxy_%d_PlumParam.txt' %galaxy_number)
	vsp_param = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/ReadData/Galaxy_%d_VSPRead.txt' %galaxy_number)
	kindat = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/ReadData/Galaxy_%d_SigP.txt' %galaxy_number)
	surfden = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/ReadData/Galaxy_%d_SurfDen.txt' %galaxy_number)
	rhalfs = np.loadtxt('/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/ReadData/GalaxiesRh.txt')


	gm, = np.where(rhalfs[:,0] == int(galaxy_number))
	rh = float(rhalfs[gm,1])


	poss = kindat[:,0]
	veldisps = kindat[:,1]
	veldisperrs = kindat[:,2]

	#return (poss, veldisps, veldisperrs, rh)
	return (poss,veldisps,veldisperrs,light_params,vsp_param, rh, surfden, mass_tot)




#a = bin_projector_read(0)
#print a

#galaxies = np.arange(31)
#for galaxy in galaxies:
#	a = bin_projector_read(galaxy)	

#poss,veldisps,veldisperrs,light_params,vsp_param, rh, surfden, mass_tot = bin_projector_read_vsp_bin_copy(0)
#print mass_tot, mass_tot.dtype
