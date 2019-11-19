import numpy as np
import input_funcs as funcs
import h5py
import matplotlib.pyplot as plt
import profiles

def galaxy_data(gal_num, outdir, data_loc):

    Nbinin_phot = 15
    Nbinin_kin = 15
    Mstar_rlim = 50.0

    
    p0in = np.array([100,100,100,0.1,0.5,0.75])
    p0in_min = np.array([10,10,10,0.05,0.05,0.05])
    p0in_max = np.array([1e4,1e4,1e4,2.0,2.0,2.0])

    tracertol = 0.5
    maxdatrad = 50.0
    maxdatfitrad = 50.0
    y_sigLOSmax = 45
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e6
    yMhigh = 1e10



    data = h5py.File(data_loc + "/%s.hdf5" %gal_num, 'r')
    pos_kin = data['KinematicsPositions'].value
    R_kin = np.sqrt(np.sum(pos_kin**2, axis = 1))/1000.0 # in kpc
    vz_kin = data['KinematicsVelocities'].value
    try:
	    vzerr_kin = data['KinVelErr'].value  # errors are all 2 km/s
    except KeyError:
            print 'Velocity error dataset does not exist :('
            print 'I will assume 2 km/s errors'
    ms_kin = data['KinematicsMasses'].value
    ms_kin = (ms_kin / np.sum(ms_kin) ) * len(ms_kin) # number relative contribution  ni/mi ~ ntot/mtot


    #Subtract space velocity:
    vz_kin = vz_kin - np.sum(ms_kin*vz_kin)/np.sum(ms_kin)  # com velocity subtracted

    pos_phot = data['PhotometryPositions'].value
    R_phot = np.sqrt(np.sum(pos_phot**2, axis = 1))/1000.0 # kpc
    ms_phot = data['PhotometryMasses'].value
    ms_phot = ms_phot / np.sum(ms_phot)*len(ms_phot)

    Mstar = data['StellarMass3R'].value  # mass within 3 Rhalf
    Mstar_err = Mstar*0.25 # assume some error ~25% of value

    if (Nbinin_phot < 0):
        Nbin = np.sqrt(np.sum(ms_phot))
    else:
        Nbin = Nbinin_phot
    print 'Number of photometric data points:', len(ms_phot)
    print 'Number of stars per photometric bin:', Nbin

    # Calculate the surface density

    rbin_phot, surfden, surfdenerr, \
            rbin_photfit, surfdenfit, surfdenerrfit, \
            Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits = \
            funcs.get_surfden_bins(R_phot,ms_phot, Nbin,maxdatrad,maxdatfitrad,p0in,p0in_min,p0in_max, Mstar_rlim, outdir, gal_num) # necessary?

 
   # And calculate the velocity disperison profile, vs1 and vs2:

    if (Nbinin_kin < 0):
        Nbin = np.sqrt(np.sum(ms_kin))
    else:
        Nbin = Nbinin_kin

    print 'Number of kinematic data points:', len(ms_kin)
    print 'Number of stars per kinematic bin:', Nbin

    nmonte = 1000  # Sampling iterations

    print 'Calculating the velocity dispersion with Nbin:', Nbin
    rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err = \
        funcs.calc_virial_moments(Rhalf,nmonte,R_kin,vz_kin,vzerr_kin,ms_kin,\
                            pfits, maxdatrad,Nbin, outdir, gal_num)

    surfden_dat = np.zeros((len(surfden), 3))
    surfden_dat[:,0] = rbin_phot
    surfden_dat[:,1] = surfden
    surfden_dat[:,2] = surfdenerr

    kin_dat = np.zeros((len(rbin_kin), 3))
    kin_dat[:,0] = rbin_kin
    kin_dat[:,1] = sigpmean
    kin_dat[:,2] = sigperr

    np.savetxt(outdir + "/%s_KinDat.txt" %gal_num, kin_dat)
    np.savetxt(outdir + "/%s_PlumParam.txt" %gal_num, pfits)
    np.savetxt(outdir + "/%s_VSPs.txt" %gal_num, np.array([vs1bin,vs1err,vs2bin, vs2err]))
    np.savetxt(outdir + "/%s_SurfDen.txt" %gal_num, surfden_dat)
    np.savetxt(outdir + "/%s_Rhalf.txt" %gal_num, [Rhalf])
    np.savetxt(outdir + "/%s_Mstar.txt" %gal_num, Mstar)

    return 0
#galaxy =0
def galaxy_data_read(gal_num, outdir):
    kin_dat = np.loadtxt(outdir + "/%s_KinDat.txt" %gal_num)
    pfits = np.loadtxt(outdir + "/%s_PlumParam.txt" %gal_num)
    vsps = np.loadtxt(outdir + "/%s_VSPs.txt" %gal_num)
    surfden_dat = np.loadtxt(outdir + "/%s_SurfDen.txt" %gal_num)
    Rhalf = np.loadtxt(outdir + "/%s_Rhalf.txt" %gal_num)
    Mstar = np.loadtxt(outdir + "/%s_Mstar.txt" %gal_num)

    Rhalf = float(Rhalf)
    Mstar = float(Mstar)


    return (kin_dat, pfits, vsps, surfden_dat, Rhalf, Mstar)
