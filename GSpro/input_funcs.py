from scipy.integrate import simps
from GSpro import fitting_funcs as fits
from GSpro import profiles
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


myfontsize = 18
mylinewidth = 1
ranalmin = 0.001
ranalmax = 500.
ranalpnts = 10000
ranal = np.logspace(np.log10(ranalmin),np.log10(ranalmax),ranalpnts) #fine bins for integration

figx = 5
figy = 5

y_sigLOSmax = 45
ymin_Sigstar = 1e1
ymax_Sigstar = 1e5
yMlow = 1e6
yMhigh = 1e10



def get_surfden_bins(R,ms,Nbin,maxdatrad,maxdatfitrad,p0in,p0in_min,p0in_max, Mstar_rlim, outdir, gal_num): #photometric positions, photometric (number)masses as input, number of stars per bin, Returns bins.
    cnt = 0 #bin id
    jsum = 0.0 #stars in bin
    norm = np.zeros(len(R))
    rbin_phot_t = np.zeros(len(R))
    bin_centres = np.zeros(len(R))
    surfden_t = np.zeros(len(R)) #temporary (number)mass keeper
    index = np.argsort(R)  # sort by position
    mean_r_bin = []
    for i in range(len(R)): # for each star
        if (jsum < Nbin): #until you reach Nbin
            surfden_t[cnt] = surfden_t[cnt] + ms[index[i]] # add masses to bin
            jsum = jsum + ms[index[i]] # add masses to get one number
            rbin_phot_t[cnt] = R[index[i]] # update max position, this is bin edges
            mean_r_bin.append(R[index[i]])
        if (jsum >= Nbin): # gone above Nbin
            norm[cnt] = jsum # total number in bin
            bin_centres[cnt] = np.median(mean_r_bin) #or mean
            if (cnt == 0):
                area = np.pi*rbin_phot_t[cnt]**2.0 # if this is bin zero
            else:
                area = np.pi*(rbin_phot_t[cnt]**2.0-rbin_phot_t[cnt-1]**2.0)         #if not bin zero

            mean_r_bin = []
            surfden_t[cnt] = surfden_t[cnt]/area #divide number by area - get surface density
            jsum = 0.0
            cnt = cnt + 1 # move on to next bin

    surfdenerr_t = surfden_t / np.sqrt(norm) # Poisson error in each bin
    rbin_phot_t = bin_centres[:cnt] # there are only cnt bins
    surfden_t = surfden_t[:cnt]
    surfdenerr_t = surfdenerr_t[:cnt]
    rbin_phot = rbin_phot_t[rbin_phot_t < maxdatrad] # bins only go as far as (50 kpc)
    surfden = surfden_t[rbin_phot_t < maxdatrad]
    surfdenerr = surfdenerr_t[rbin_phot_t < maxdatrad]
    rbin_photfit = rbin_phot_t[rbin_phot_t < maxdatfitrad]
    surfdenfit = surfden_t[rbin_phot_t < maxdatfitrad]
    surfdenerrfit = surfdenerr_t[rbin_phot_t < maxdatfitrad]  # fitting purposes, exclude some really outer radius

    pfits = fits.tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_photfit,surfdenfit,surfdenerrfit) #get 3 plummer fit. p0 - priors and initial values

    Mstar_rad = rbin_phot
    norm = np.max(profiles.threeplummass(np.linspace(0,Mstar_rlim,100),\
                                pfits[0],pfits[1],\
                                pfits[2],\
                                pfits[3],pfits[4],pfits[5]))   # what is the maximum mass (at Mstar_rlim)?
    Mstar_prof = profiles.threeplummass(Mstar_rad,pfits[0],pfits[1], \
                               pfits[2],\
                               pfits[3],pfits[4],pfits[5])  / norm  # mass at each photometric bin, normalized
    Mstar_surf = profiles.threeplumsurf(ranal,pfits[0],pfits[1],\
                               pfits[2],\
                               pfits[3],pfits[4],pfits[5]) / norm   # fitted surface density at each photometric bin

    Mcum_surf = 0.0
    i = 1

    print('Norm ::', norm,pfits[0]+pfits[1]+pfits[2])

    while (Mcum_surf < (pfits[0]+pfits[1]+pfits[2])/2.0/norm): #  cum_mass/tot_mass = 0.5  compute the half-mass radius

        Mcum_surf =   Mcum_surf + \
                    2.0*np.pi*ranal[i]*Mstar_surf[i]*(ranal[i]-ranal[i-1]) #  Sum 2 pi r M dr
        #Alternatively : use Simpson's rule  Mcum_surf =  Mcum_surf + integrator(Mstar_surf[:i] * 2 * np.pi * ranal[:i],ranal[:i])
        i = i + 1
    Rhalf = ranal[i-1]
    print('Rhalf calculated: ', Rhalf)

    return rbin_phot, surfden, surfdenerr, \
        rbin_photfit, surfdenfit, surfdenerrfit, \
        Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits


def calc_virial_moments(Rhalf,nmonte,R,vz,vzerr,ms,pfits,maxdatrad,Nbin, outdir, gal_num):
    #Improve estimator using fitted surfden:
    rint = np.logspace(np.log10(Rhalf/100.0),\
                       np.log10(Rhalf*1000.0),10000)
    index = np.argsort(R) # sort all kinematics positions by distance

    #First calculate and subtract the mean vz:  ms ~ number contribution

    vzmean = np.sum(vz*ms)/np.sum(ms)   #com velocity
    vzmeanerr = 0.

    #And now the 2nd and 4th moments:
    cnt = 0
    jsum = 0.0
    norm = np.zeros(len(R))
    vlos4med = np.zeros(len(R))
    vlos2med = np.zeros(len(R))
    rbin_tmp = np.zeros(len(R))
    bin_centres = np.zeros(len(R))
    mean_r_bin = []
    for i in range(len(R)):       #fill bins with specified stars per bin
        if (jsum < Nbin):
            vlos4med[cnt] = vlos4med[cnt] + \
                (vz[index[i]]-vzmean)**4.*ms[index[i]]    # number weighted 4th moment
            vlos2med[cnt] = vlos2med[cnt] + \
                (vz[index[i]]-vzmean)**2.*ms[index[i]]    # number weighted 2nd moment
            rbin_tmp[cnt] = R[index[i]]                        # where is the furthest star in bin?
            mean_r_bin.append(R[index[i]])
            jsum = jsum + ms[index[i]]
        if (jsum >= Nbin):
            norm[cnt] = jsum
            bin_centres[cnt] = np.median(mean_r_bin) #or mean
            jsum = 0.0
            mean_r_bin = []
            cnt = cnt + 1
    vlos4med = vlos4med[:cnt]
    vlos2med = vlos2med[:cnt]
    norm = norm[:cnt]
    vlos4med = vlos4med / norm     # normalize with tot stars
    vlos2med = vlos2med / norm
    rbin_tmp = bin_centres[:cnt]

    #And Monte-Carlo the errors:
    vlos4 = np.zeros((nmonte,len(R)))    #row each iteration, col each bin
    vlos2 = np.zeros((nmonte,len(R)))
    vlos4_pureerr = np.zeros((nmonte,len(R)))
    vlos2_pureerr = np.zeros((nmonte,len(R)))
    norm = np.zeros(len(R))
    for k in range(nmonte):
        cnt = 0
        jsum = 0.0
        for i in range(len(R)):
            vz_err = (vz[index[i]]-vzmean)+\
                np.random.normal(0.0,vzerr[index[i]])   #perturb velocity with 2 km/s
            vz_pure_err = np.random.normal(0.0,vzerr[index[i]])
            if (jsum < Nbin):
                vlos4[k,cnt] = vlos4[k,cnt] + \
                    vz_err**4.*ms[index[i]]
                vlos2[k,cnt] = vlos2[k,cnt] + \
                    vz_err**2.*ms[index[i]]
                vlos4_pureerr[k,cnt] = vlos4_pureerr[k,cnt] + \
                    vz_pure_err**4.*ms[index[i]]               #4th moment of the error?
                vlos2_pureerr[k,cnt] = vlos2_pureerr[k,cnt] + \
                    vz_pure_err**2.*ms[index[i]]
                jsum = jsum + ms[index[i]]
            if (jsum >= Nbin):
                norm[cnt] = jsum
                jsum = 0.0
                cnt = cnt + 1

    vlos4tmp = np.zeros((nmonte,cnt))
    vlos4tmp = vlos4[:,:cnt]
    vlos2tmp = np.zeros((nmonte,cnt))
    vlos2tmp = vlos2[:,:cnt]
    vlos4_pe_tmp = np.zeros((nmonte,cnt))
    vlos4_pe_tmp = vlos4_pureerr[:,:cnt]
    vlos2_pe_tmp = np.zeros((nmonte,cnt))
    vlos2_pe_tmp = vlos2_pureerr[:,:cnt]
    norm = norm[:cnt]

    vlos4 = vlos4tmp / norm
    vlos2 = vlos2tmp / norm
    vlos4_pe = vlos4_pe_tmp / norm
    vlos2_pe = vlos2_pe_tmp / norm

    #And now estimate the full measurement error:
    vlos4err_meas = np.zeros(cnt)
    vlos2err_meas = np.zeros(cnt)
    vlos4_pe_meas = np.zeros(cnt)
    vlos2_pe_meas = np.zeros(cnt)

    for k in range(cnt):
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = fits.calcmedquartnine(vlos4[:,k])
        vlos4err_meas[k] = (sixhigh-sixlow)/2.0                           #mean of 16th and 84th
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = fits.calcmedquartnine(vlos2[:,k])
        vlos2err_meas[k] = (sixhigh-sixlow)/2.0
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = fits.calcmedquartnine(vlos4_pe[:,k])
        vlos4_pe_meas[k] = (sixhigh-sixlow)/2.0
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = fits.calcmedquartnine(vlos2_pe[:,k])
        vlos2_pe_meas[k] = (sixhigh-sixlow)/2.0

    #Combine with the Poisson error:
    vlos4err = np.sqrt(vlos4err_meas**2.0 + vlos4med**2.0/Nbin)
    vlos2err = np.sqrt(vlos2err_meas**2.0 + vlos2med**2.0/Nbin)

    vlos4med = vlos4med - vlos4_pe_meas  # subtract pure error
    vlos2med = vlos2med - vlos2_pe_meas

    #Demand positive:

    print('Dispersion profile:', np.sqrt(vlos2med))


    vlos2med = vlos2med[vlos4med > 0]
    vlos2err = vlos2err[vlos4med > 0]
    vlos4err = vlos4err[vlos4med > 0]
    rbin_tmp = rbin_tmp[vlos4med > 0]
    vlos4med = vlos4med[vlos4med > 0]



    p0in = np.array([1.0,0.0,1000.0])
    p0in_min = np.zeros(3)
    p0in_max = np.array([1e5,1.0,1e5])

    fitmin = int(len(rbin_tmp)/2)   # half of data range, not quite rhalf
    pfitmed = fits.v4fit(p0in,p0in_min,p0in_max,\
                        rbin_tmp[fitmin:],vlos4med[fitmin:],\
                        vlos4err[fitmin:])    #not used


    #Cut back to maxdatrad:
    rbin_tmp_full = rbin_tmp
    vlos4err_full = vlos4err
    vlos2err_full = vlos2err
    vlos4med_full = vlos4med
    vlos2med_full = vlos2med
    vzmean_full = vzmean
    vzmeanerr_full = vzmeanerr
    vlos4err = vlos4err[rbin_tmp < maxdatrad]
    vlos2err = vlos2err[rbin_tmp < maxdatrad]
    vlos4med = vlos4med[rbin_tmp < maxdatrad]
    vlos2med = vlos2med[rbin_tmp < maxdatrad]
    rbin_tmp = rbin_tmp[rbin_tmp < maxdatrad]

    #Fit straight line in log-space to vlos4med:

    pfits_powline = fits.fit_powline(rbin_tmp,vlos4med,vlos4err,Rhalf) #fit to above Rhalf
    router_min = np.max(rbin_tmp)    #  Rdata  <  Rout < 2 Rdata
    router_max = 2.0*router_min

    if (router_max < router_min):
        router_max = router_min
    gamout_min = 1.0         # Slope outside Rout
    gamout_max = 3.0
    router = (router_min+router_max)/2.0     #start in the middle
    gamout = (gamout_min+gamout_max)/2.0


    #plt.figure()
    #plt.loglog()
    #plt.errorbar(rbin_tmp_full,vlos4med_full,vlos4err_full,color='k')
    #plt.errorbar(rbin_tmp,vlos4med,vlos4err,color='b')
    #tvl4 = fits.tvl4func(rint,rbin_tmp,vlos4med,\
    #                pfits_powline[0],pfits_powline[1],\
#                    pfits_powline[2],router,gamout)
#    plt.plot(rint,tvl4,'r')
#    plt.show()


    #Plot the 1st, 2nd and 4th moments:
    if (len(rbin_tmp_full) > 1):
        #Calculate sigLOS(Rhalf) and output:
        if (np.max(rbin_tmp_full) > Rhalf):
            j=0.0
            while (rbin_tmp_full[j] < Rhalf):
                j=j+1
            print('sigLOS(Rhalf) [km/s]:', np.sqrt(vlos2med_full[j]))

        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)

        plt.loglog()
        plt.errorbar(rbin_tmp_full,vlos4med_full,vlos4err_full,color='black')
        plt.errorbar(rbin_tmp,vlos4med,vlos4err,color='blue')
        tvl4 = fits.tvl4func(rint,rbin_tmp,vlos4med,\
                        pfits_powline[0],pfits_powline[1],\
                        pfits_powline[2],router,gamout)
        plt.plot(rint,tvl4,'red')
        plt.plot(rint,pfits_powline[0]*(rint/np.max(rbin_tmp))**\
                 pfits_powline[1]+\
                 pfits_powline[2],'green')
        plt.xlim([Rhalf/10.0,Rhalf*100.0])
        plt.ylim([1.0,y_sigLOSmax**4.0*100.0])
        plt.xlabel(r'$R\,[{\rm kpc}]$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$\langle v_{\rm los}^4\rangle\,({\rm km}^4\,{\rm s}^{-4})$',\
                       fontsize=myfontsize)
        plt.savefig(outdir+'Galaxy_%s_output_vlos4.pdf' % gal_num,bbox_inches='tight')
        plt.close()
        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)

        plt.loglog()
        plt.errorbar(rbin_tmp_full,np.sqrt(vlos2med_full),\
                         vlos2err_full/2.0/np.sqrt(vlos2med_full),color='k')
        plt.errorbar(rbin_tmp,np.sqrt(vlos2med),\
                         vlos2err/2.0/np.sqrt(vlos2med),color='b')
        plt.xlim([Rhalf/10.0,Rhalf*100.0])
        plt.ylim([1.0,y_sigLOSmax])

        plt.xlabel(r'$R\,[{\rm kpc}]$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$\sigma_{\rm LOS}\,({\rm km}\,{\rm s}^{-1})$',\
                       fontsize=myfontsize)
        plt.savefig(outdir+'Galaxy_%s_output_vlos2.pdf' % gal_num,bbox_inches='tight')
        plt.close()
    #And calculate vs1 and vs2:
    tvl4 = fits.tvl4func(rint,rbin_tmp,vlos4med,\
                    pfits_powline[0],pfits_powline[1],\
                    pfits_powline[2],router,gamout)
    test_surfden = profiles.threeplumsurf(rint,pfits[0],pfits[1],pfits[2],\
                                     pfits[3],pfits[4],pfits[5])
    vs1imp = simps(tvl4*test_surfden*rint,rint)
    vs2imp = simps(tvl4*test_surfden*rint**3.0,rint)

    #Monte-Carlo to calculate the ~1sigma errors:
    vs1_samp = np.zeros(nmonte)
    vs2_samp = np.zeros(nmonte)
    for i in range(nmonte):
        vlos4_samp = vlos4med + \
            np.random.normal(0.0,vlos4err,len(vlos4med))   # sample the error

        #Fit straight line in log-space for this draw:
        pfits_powline = fits.fit_powline(rbin_tmp,vlos4_samp,vlos4err,Rhalf)

        #Draw router and gamout:
        #router = np.random.uniform(router_min, router_max)     #marginalize
        #gamout = np.random.uniform(gamout_min, gamout_max)
        router = np.random.random()*(router_max-router_min)+router_min
        gamout = np.random.random()*(gamout_max-gamout_min)+gamout_min

        tvl4 = fits.tvl4func(rint,rbin_tmp,vlos4_samp,\
                        pfits_powline[0],pfits_powline[1],\
                        pfits_powline[2],router,gamout)
        vs1_samp[i] = simps(tvl4*test_surfden*rint,rint)
        vs2_samp[i] = simps(tvl4*test_surfden*rint**3.0,rint)

    median, sixlow, sixhigh, ninelow, ninehigh,\
        nineninehigh, nineninelow = fits.calcmedquartnine(vs1_samp)
    vs1imperr = (sixhigh-sixlow)/2.0
    median, sixlow, sixhigh, ninelow, ninehigh,\
        nineninehigh, nineninelow = fits.calcmedquartnine(vs2_samp)
    vs2imperr = (sixhigh-sixlow)/2.0

    print('VirialShape vs1:', vs1imp,vs1imperr)
    print('VirialShape vs2:', vs2imp,vs2imperr)

    vs1bin = vs1imp
    vs2bin = vs2imp
    vs1err = vs1imperr
    vs2err = vs2imperr

    #Output also 2nd moment (pruning any would-be NaN values):
    rbin_kin = rbin_tmp[vlos2med > 0]
    sigpmean = np.sqrt(vlos2med[vlos2med > 0])
    sigperr = vlos2err[vlos2med > 0]/2.0/\
        np.sqrt(vlos2med[vlos2med > 0])

    #And finally calculate (for comparison only) the
    #Richardson & Fairbairn estimators:
    zeta_A = np.float(len(R))*np.sum(vz**4.0)/np.sum(vz**2.0)**2.0
    zeta_B = np.float(len(R))**2.0*np.sum(vz**4.0*R**2.0)/\
        (np.sum(vz**2.0)**2.0*np.sum(R**2.0))
    print('Richardson+Fairbairn estimators:')
    print('Nstars, zeta_A, zeta_B', len(R), zeta_A, zeta_B)

    mean_disp = np.sum(sigpmean)/np.float(len(sigpmean))
    print('Mean dispersion:', mean_disp)
    print('Mean dispersion error:', np.sqrt(np.sum((sigpmean-mean_disp)**2.0)/\
        np.float(len(sigpmean)-1))/np.sqrt(len(sigpmean)))

    #plt.errorbar(rbin_kin,sigpmean, sigperr)
    #plt.show()

    print('Started plotting')

    plt.figure(figsize = (5,5))
    plt.hist(vs1_samp,bins = 50, histtype = 'step')
    plt.axvline(vs1bin)
    plt.axvline(vs1bin + vs1err)
    plt.axvline(vs1bin - vs1err)
    plt.xlabel('v1', fontsize = 16)
    plt.savefig(outdir+'Galaxy_%s_v1hist.pdf' % gal_num,bbox_inches='tight')
    plt.close()

    plt.figure(figsize = (5,5))	
    plt.hist(vs2_samp, bins = 50, histtype = 'step')
    plt.axvline(vs2bin)
    plt.axvline(vs2bin + vs2err)
    plt.axvline(vs2bin - vs2err)
    plt.xlabel('v2', fontsize = 16)
    plt.savefig(outdir+'Galaxy_%s_v2hist.pdf' % gal_num,bbox_inches='tight')
    plt.close()
    	
    print('Finished plotting')

    return rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err
