# pyGravSphere

**pyGravSphere** is a non-parametric spherical Jeans analysis code. In essence, it is a wrapper for **EmCee**, designed for Jeans analysis, with a friendly and easy-to-use interface (I hope).

Below is my attempt at a tutorial. If any questions or suggestions arise, do feel free to contact *agenina@mpa-garching.mpg.de*.


**Note from 26/05/2020: when using the in-built preprocessing functionality pyGravsphere now defaults to <img src="https://render.githubusercontent.com/render/math?math=N/\sqrt{N}"> stars per bin.**


**A more sophisticated pyGravSphere-based binning routine for the second velocity moments has been introduced in Zoutendijk et al. 2021(https://arxiv.org/abs/2101.00253)**

**This can be found here: https://github.com/slzoutendijk/hkbin**

## Getting Started

To download **pyGravSphere**, please go to your preferred folder and type:

```
git clone git@github.com:AnnaGenina/pyGravSphere.git
```


Before you begin using **pyGravSphere** please make sure you have downloaded and installed the following dependencies ...and their dependencies.

* gcc compiler
* numpy (duh)
* scipy
* matplotlib
* EmCee https://emcee.readthedocs.io/en/latest/user/install/   (version 2.2.1 had been used)
* corner https://corner.readthedocs.io/en/latest/ (corner plot making code)
* SWIG http://www.swig.org/download.html
* h5py https://pypi.org/project/h5py/
* lmfit https://lmfit.github.io/lmfit-py/

If you will be using multiple processors

* schwimmbad https://github.com/adrn/schwimmbad  
* gnucomp & OpenMPI or similar, if you're running on a cluster

## Installing

After installing you dependencies please go to the **pyGravSphere** install folder and type:

```
python setup.py build_ext --inplace
```

If you don't get any error messages, congratulations! 

Now go to your system's *.bashrc* file, or an equivalent, and define your GravSphere path as the install directory

```
export GravSpherePath=/Users/Nerd/Desktop/pyGravSphere/
```
Success!

Now start a new terminal window for the changes to take place.

## Practical

This is a good time to decide where you want to keep all of your **GravSphere** input data and outputs. Time to pick a working directory. To select a woking directory, please update *workdir.txt* inside the install folder with your working directory path i.e. something like

```
/Users/Nerd/Desktop/MyWorkDir/
```
Do make sure that you've previously created MyWorkDir. 

Good job!

Now run **pyGravSphere** for the first time inside your install directory.

```
python pygravsphere.py
```

This will generate all the necessary folders and files you need at this step. Quit **pyGravSphere** by typing *q* or *quit*.

### Preparing a data file

**pyGravSpere** currently pre-processes input data in *.hdf5* format. In just a few short steps *.hdf5* files will become your new best friends! 

You may have two data sets: one for photometry (positions) of the stars and one for kinematics (positions and velocities). In the following we'll assume that these are the same two datasets.

So suppose you have a 5-column data file with x,y positions (***in parsecs***) on the sky, line-of-sight velocities with their errors and relative weights of each star. 

For the purposes of testing, we recommend you try out one of the Gaia Challenge datasets http://astrowiki.ph.surrey.ac.uk/dokuwiki/doku.php?id=tests:sphtri:spherical

To create a **pyGravSphere**-compatible *.hdf5* file use the following Python code:

```
import numpy as np
import h5py

f_old = np.loadtxt("gaiachallengedata.txt") # x y vz vzerr
nstars = len(f_old)

xy = f_old[:,[0,1]]*1000  #convert to parsecs if needed
vz = f_old[:,5]
vzerr = 2.
weights = np.ones((nstars,))  # Change to actual stellar masses or luminosities if needed (not magnitudes)
total_mass = 1  # Change to mass within 3 half-light radii / total mass if appropriate

f_new = h5py.File("GaiaPlumCuspIso.hdf5", 'w')

dset1 = f_new.create_dataset('PhotometryPositions',(nstars,2), dtype = xy.dtype)
dset1[...] = xy  #parsecs!!!
dset2 = f_new.create_dataset('PhotometryMasses', (nstars,), dtype = xy.dtype)
dset2[...] = weights
dset3 = f_new.create_dataset('KinematicsPositions', (nstars,2), dtype = xy.dtype)
dset3[...] = xy #parsecs!!!
dset4 = f_new.create_dataset('KinematicsVelocities', (nstars,), dtype = xy.dtype)
dset4[...] = vz
dset5 = f_new.create_dataset('KinVelErr', (nstars,), dtype = xy.dtype)
dset5[...] = vzerr
dset6 = f_new.create_dataset('KinematicsMasses', (nstars,), dtype = xy.dtype)
dset6[...] = weights
dset7 = f_new.create_dataset('StellarMass3R', (1,), dtype = xy.dtype)
dset7[...] = total_mass

f_new.close()


```
Change the contents of each data set as appropriate. The *Photometry* dataset will be used exclusively for calculating the best fit 3-component Plummer profile. The *Kinematics* dataset will be used for calculation of the binned velocity dispersion profile and the virial shape parameters. If the mass in stars is assumed to significantly contribute to your system, put in your best estimate of mass **in solar masses** in to *[StellarMass3R]* set. In the above example, each star contributes equal weight to your system. This is a good approximation for real-life galaxies. If you're applying GravSphere to a cosmological simulation and you know the masses of your stellar particles, you can use those. This will be useful in the **Pre-processing** step, where **pyGravSphere** calculates the binned velocity dispersion and stellar density profiles. If you calculate those things using your own method, then don't worry about this.

Now you have an *GaiaPlumCuspIso.hdf5* dataset. You might even have three data sets: *Galaxy_1.hdf5*, *Galaxy_2.hdf5* and *Galaxy_3.hdf5*. At this point, please open the *galaxy_list.txt* file in your working directory and type in the galaxy names: 

```
Galaxy_1
Galaxy_2
Galaxy_3
```

Please deposit you new data files into the *KinPhotDat* directory that is now inside your working directory.

### Pre-processing your data

Excellent! Now you're ready to create your MCMC input files. This includes binned velocity dispersion, surface brightness, half-light radius and two virial shape parameters. You might even have a different way of computing these, in which case you should output your pre-processed data in the same format as we do here.

*Note that, after the pre-processing, the positions of the stars are now in **kiloparsecs** .*

Go to your working directory and type:

```
python pygravsphere.py
```
You will see a list of your available options. We begin with pre-processing. Type:

```
0
```
You will be prompted to check if the working directory is correct. 

```
y
```

Now you will be asked which datesets you would like to pre-procces. You can pick *All*, which will do all three, or you can specify particular one you want.

So pick *Specify* 

```
2
```

and type your galaxy names:

```
Galaxy_3, Galaxy_1
```
The output files will make their way to the *GalaxyData* directory inside your working directory. If you would like to process your data yourself, make sure the format is the same as here.


**Note that a more sophisticated pyGravSphere-based binning routine for the second velocity moments has been introduced in Zoutendijk et al. 2021**

**This can be found here: https://github.com/slzoutendijk/hkbin**

### Preparing submission scripts

Ok, so this is the point where I emphasize that your life will be easier if you have a few cores to run **pyGravSphere** on. Like, much easier. 

#### If your institution has a cluster you can run pyGravSphere on

In your working directory, open *sub_script.txt* file. Fill the file with what your batch system submission script looks like. For example, my institution uses *slurm* and my *sub_script.txt* file looks like:

```
#!/bin/bash -l
#SBATCH --ntasks=CORENUM
#SBATCH -J GALID
#SBATCH -e GALID.err
#SBATCH -o GALID.out
#SBATCH -p cosma
#SBATCH -A durham
#SBATCH --time=TIME
module purge
module load python
module load gnu_comp/7.3.0 openmpi/3.0.1
```

For the name of the job, output file, error file, number of cores and time please keep the same names (i.e. CORENUM, GALID) so that **pyGravSphere** can replace those for you automatically. If your input scripts are vastly different to this, you might need to fiddle with the *gsTools.py* **create_sub** and **create_ana_sub** functions. Let's hope it doesn't come to that.

#### If you're running GravSphere on your laptop

Leave the *sub_script.txt* file empty.


### Creating a project

OK, we're definitely getting closer to the good stuff.

So suppose now you want to run a model with the Zhao et al. (1996) dark matter distribution, a constant anisotropy profile, 3-component Plummer fit and no virial shape parameters. And, very conveniently, you have 8 cores to run it on.

Note: if you're running on your desktop/laptop, use the multiprocessing option for faster results. Use MPI if you're running on more than one node. 


Now run **pyGravSphere** 
```
python pygravsphere.py
```
And here is your example selection:

```
Option : 1
1
The current working directory is /Users/Nerd/Desktop/MyWorkDir/
Would you like to change that?  y or n?
y or n : n
Alright then
Please name your new project: ZhaoConst3Plum
Please write a short summary (or long if you want): Tutorial
Creating project directory
Creating Submission files
How would you like to run MCMC?
1) Serial 	2) Multiprocessing  3) MPI (if you're running on a cluster) 3
How many cores? 8
How much time do you need (in hours)? 2.5
1) PowerLaw or 2) Zhao ? 2
1) Baes or 2) Const ? 2
VSP ? y or n? n
1) Plummer Single  2) Plummer Components 2
How many walkers? 500 
Burn-in? 5000
Steps? 10000
Integration points? 100
```

Think long and hard about how many *walkers* you want. About 600 is good.
*Burn-in* is the number of steps for which you want to run each of your walkers before **pyGravSphere** starts outputting the chains. This is approximately the number of steps needed for the result to converge.

*Steps* is the total steps for which you want to run your walkers. So in the example above, you will get an output of the last 5000 steps per walker.

*Integration points* sets the accuracy of your integrator: how many logarithmic bins will the integrand be divided into? 100 is fine.

Consider the amount of time you will need for the job to run. The above configuration runs in just under 1 hour on 8 cores.

## Priors

Alright! Your project has been created and your jobs are ready for submission. Now if you go to:

```
/Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum/Submissions/
```

you will find your submission scripts as well as the *priors.txt* file. This is full of default **pyGravSphere** priors. Now you might want to change those. To do so, simply edit the lower and upper priors on each parameter. If you decide you want to keep any of the parameters constant, simply edit the second last column and make sure "min" and "max" values are the same. So, for example, if I want to keep the inner slope *gamma* at 1, I will edit the priors:

```
rhos rhosmin rhosmax rhosconst 5 10. False VAR
rs rsmin rsmax rsconst -2 1.5 False VAR
alpha alphamin alphamax alphaconst 0.5 3. False VAR
beta betamin betamax betaconst 3. 7. False VAR
gamma gammamin gammamax gammaconst 1. 1. True VAR
beta0 beta0min beta0max beta0const -1. 1. False VAR
m1 m1min m1max m1const 0.5 1.5 False np.log10(VAR*lightpower[0])
a1 a1min a1max a1const 0.5 1.5 False VAR*lightpower[3]
m2 m2min m2max m2const 0.5 1.5 False np.log10(VAR*lightpower[1])
a2 a2min a2max m2const 0.5 1.5 False VAR*lightpower[4]
m3 m3min m3max m3const 0.5 1.5 False np.log10(VAR*lightpower[2])
a3 a3min a3max a3const 0.5 1.5 False VAR*lightpower[5]
mstar mstarmin mstarmax mstarconst 0.75 1.25 False np.log10(VAR*stellar_mass)

```

Excellent!

### Running the jobs

#### On a cluster
To run our imaginary *Galaxy_1* we can now go to

```
/Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum/Submissions/
```

and type, depending on you batch system, 

```
sbatch Galaxy_1.sh
```

or

```
qsub < Galaxy_1.sh
```

or 

```
bsub Galaxy_1.sh
```

Alternatively you can run **pyGravSphere** from your working directory. Type *pygravsphere.py* and select option 2 -- Submit jobs. 

#### On your laptop

Go to 
```
/Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum/Submissions/
```
and type

```
./Galaxy_1.sh

```

Now you wait.

**pyGravSphere** will deposit the output every 100 steps into the */Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum/Galaxy_1/* directory. You can monitor the progress / errors in the */Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum/OutErr/* directory.

### If your job didn't finish on time

Don't worry! Just modify the submission script *Galaxy_1.sh*, replacing *standard* to *restart* and your run will start from where it left off! Neat, eh?

### If you want to run your chains for longer

Replace *standard* in the *Galaxy_1.sh* submission script to *continue*. Change the the number of steps to whatever you like it to be. Do not change the number of walkers!


## Analysis

Now that wasn't so bad, was it? 

**pyGravSphere** will now have output your chains into the */Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum/Galaxy_X/* directory.

At this point you can do what you like with the output.

**pyGravSphere** comes with some built-in functionality to get you started.

### Get median, 68, 95 percentiles

These will be output into */Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum/Analysis/Limits/* after you've ran the Analysis option (option 3)


### Plot the evolution of your chains

This might be useful for visual evaluation of convergence. (option 4)
Note that you cannot use this when you've used multiprocessing, as the chains are not output in any particular order.


### Corner plots

You might want to output a corner plot of your posteriors. This will make a corner plot of all of your free parameters. (option 5)


### Output the percentiles of your mass/density/anisotropy models

Use option 6. For each radial distance this will output 1-100 percentile masses, densities and anisotropy. Make sure you've run Analysis previously.


## Publications using pyGravSphere 

So you got an awesome result? Great! ***You're welcome!*** Please don't foget to cite the original **GravSphere** paper:

**https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4541R/abstract**

as well as the paper that brought you here:

**https://ui.adsabs.harvard.edu/abs/2020MNRAS.tmp.2435G/abstract**

Thanks!

## Issues

Contact me at *agenina@mpa-garching.mpg.de*. If it's wrong, we'll fix it.

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details

## Future developments

The era of *Python 2* is coming to an end and we will be soon upgrading to *Python 3* as well as the latest version of **EmCee**.

In terms of the functionality, the updates will include:

* Flexible number of power law bins
* Multiple tracer populations
* Proper motions data
* Accounting for rotation

Requests? Let us know.

## Acknowledgments

**pyGravSphere** relies heavily on the **EmCee** code. Make sure you are familiar with how it works! Don't forget to cite:

**https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract**
