import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt("/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/Read/ImpNum/Analysis/Limits/" + "Galaxy_0_PlummerLims.txt")

read = np.loadtxt("g0sigstar.txt")

comps = np.loadtxt("/cosma5/data/dp004/dc-geni1/project/Jfactor4/EmCee/Samples/JAna/ReadData/" + "Galaxy_0_PlumParam.txt")


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


	
	return (sums_proj)



fig,ax1 = plt.subplots(1,1, figsize = (6,6))

ax1.plot((dat[:,0]) , 10**(dat[:,1]) )
ax1.fill_between((dat[:,0]) , 10**(dat[:,2]) ,10**(dat[:,3]), color = 'grey', alpha = 0.7 )
ax1.fill_between((dat[:,0]) , 10**(dat[:,4]) ,10**(dat[:,5]), color = 'grey', alpha = 0.5 )

plt.plot(read[:,0], read[:,1], color = 'red', lw = 0.5)

plt.plot(dat[:,0],  plummer_proj_comp(dat[:,0],comps[0],comps[1]), ls = 'dashed')
plt.plot(dat[:,0],  plummer_proj_comp(dat[:,0],comps[2],comps[3]), ls = 'dashed')
plt.plot(dat[:,0],  plummer_proj_comp(dat[:,0],comps[4],comps[5]), ls = 'dashed')

plt.plot(dat[:,0], plummer_proj_sum(comps, dat[:,0], 3), color = 'blue')

surfden = np.loadtxt('read_surfden.txt')

plt.plot(surfden[:,0], surfden[:,1], color = 'orange', lw = 0.3)

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim([10,10**5])
ax1.set_xlim([10**(-2), 3])



plt.show()
