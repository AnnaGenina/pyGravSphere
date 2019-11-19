# pyGravSphere

**pyGravSphere** is a non-parametric spherical Jeans analysis code, ideal for applications to dwarf spheroidal kinematic data.

We've tried our best to design a user-friendly interface for **pyGravSphere** and below we include a full tutorial, however, if any questions or suggestions arise, do feel free to contact *anna.genina@durham.ac.uk*.



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
* OpenMP (if you're using a cluster)

## Installing

After installing you dependencies please go to the **pyGravSphere** install folder and type:

```
python setup.py build_ext --inplace
```

If you don't get any error messages, congratulations! 

Now go to your system's *.bashrc* file, or an equivalent, and define your GravSphere install directory

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

Good job!

Now run **pyGravSphere** for the first time inside your install directory.

```
python pygravsphere.py
```

This will generate all the necessary folders and files you need at this step. Quit **pyGravSphere** by typing *q* or *quit*.

### Preparing a data file

**pyGravSpere** currently pre-processes input data in *.hdf5* format. In just a few short steps *.hd5* files will become your new best friends! 

You may have two data sets: one for photometry (positions) of the stars and one for kinematics (positions and velocities). In the following we'll assume that these are the same two datasets.

So suppose you have a 3-column data file with x,y positions (***in parsecs***) on the sky and line-of-sight velocities with their errors. To create a **pyGravSphere**-compatible *.hdf5* file use the following Python code:

```
import numpy as np
import h5py

f_old = np.loadtxt("OldDataFile.txt") # x y vz vzerr
n_stars = len(f_old)

xy = f_old[:,[0,1]]
vz = f_old[:,2]
vzerr = f_old[:,3]
weights = np.ones((n_stars,))  # Change to actual stellar masses or luminosities if needed (not magnitudes)
total_mass = 1  # Change to mass within 3 half-light radii / total mass if appropriate

f_new = h5py.File("NewDataFile.hdf5", 'w')

dset1 = output.create_dataset('PhotometryPositions',(nstars,2), dtype = xy.dtype)
dset1[...] = xy  #parsecs!!!
dset2 = output.create_dataset('PhotometryMasses', (nstars,), dtype = xy.dtype)
dset2[...] = weights
dset3 = output.create_dataset('KinematicsPositions', (nstars,2), dtype = xy.dtype)
dset3[...] = xy #parsecs!!!
dset4 = output.create_dataset('KinematicsVelocities', (nstars,), dtype = xy.dtype)
dset4[...] = vz
dset5 = output.create_dataset('KinVelErr', (nstars,), dtype = xy.dtype)
dset5[...] = vzerr
dset6 = output.create_dataset('KinematicsMasses', (nstars,), dtype = xy.dtype)
dset6[...] = weights
dset7 = output.create_dataset('StellarMass3R', (1,), dtype = xy.dtype)
dset7[...] = total_mass

f_new.close()

```
Obviously, change the contents of each data set as appropriate. The *Photometry* dataset will be used exclusively for calculating the best fit 3-component Plummer profile. The *Kinematics* dataset will be used for calculation of the binned velocity dispersion profile and the virial shape parameters. If the mass in stars is assumed to significantly contribute to your system, put in your best estimate of mass **in solar masses** in to *[StellarMass3R]* set. In the above example, each star contributes equal weight to your system. This is a good approximation for real-life galaxies. If you're applying GravSphere to a cosmological simulation and you know the masses of your stellar particles, you can use those. This will be useful in the **Pre-processing** step, where **pyGravSphere** calculates the binned velocity dispersion and stellar density profiles. If you calculate those things using your own method, then don't worry about this.

Now you have an *NewDataFile.hdf5* dataset. Perhaps you want to name it as something more identifiable, let's say *Galaxy_1.hdf5*, or you might even have three data sets: *Galaxy_1.hdf5*, *Galaxy_2.hdf5* and *Galaxy_3.hdf5*. At this point, please open the *galaxy_list.txt* file in your working directory and type in: 

```
Galaxy_1
Galaxy_2
Galaxy_3
```

Please deposit you new data files into the *KinPhotDat* directory that is now inside your working directory.

### Pre-processing your data

Excellent! Now you're ready to create your MCMC input files. This includes binned velocity dispersion, surface brightness, half-light radius and two virial shape parameters. You might even have a different way of computing these, in which case I would advise to output your pre-processed data in the same format as we do here.

*Note that, after the pre-processing, the positions of the stars are now in **kiloparsecs** .*

Go to your working directory and type:

```
pygravsphere.py
```
You will see a list of your available options. We begin with pre-processing. Type:

```
0
```
You will be prompted to check if the working directory is correct. If you actually followed this tutorial, it should be.

```
y
```

Now you will be asked which datesets you would like to pre-procces. You can pick *All*, which will do all three, or you can specify particular one you want.

So pick *Specify* and type:

```
Galaxy_3, Galaxy_1
```
Bingo! Cogs are turning and your data will be pre-processed in no time. The output files will make their way to the *GalaxyData* directory inside your working directory.

### Preparing submission scripts

Ok, so this is the point where I re-iterate that your life will be easier if you have a few cores to run **pyGravSphere** on. Like, much easier. 

#### If your institution has a cluster you can run GravSphere on

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
module load gsl
```

For the name of the job, output file, error file, number of cores and time please keep the same names (i.e. CORENUM, GALID) so that **pyGravSphere** can replace those for you automatically. If your input scripts are vastly different to this, you might need to fiddle with the *gsTools.py* **create_sub** and **create_ana_sub** functions. Let's hope it doesn't come to that.

#### If you're running GravSphere on your laptop

Leave the *sub_script.txt* file empty.


### Creating a project

OK, we're definitely getting closer to the good stuff.

So suppose now you want to run a model with the Zhao et al. (1996) dark matter distribution, a constant anisotropy profile, 3-component Plummer fit and no virial shape parameters. And, very conveniently, you have 8 cores to run it on.

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
Please name your new project: ZhaoConst3Plum
MPI? y or n?
y or n: y
How many cores? 8
How much time do you need (hours only)? 5
1) PowerLaw or 2) Zhao ? 2
1) Baes or 2) Const ? 2
VSP ? y or n? n
1) Plummer Single  2) Plummer Components 2
How many walkers? 1000
Burn-in? 2500
Steps? 5000
Integration points? 100
Creating project directory
mkdir /Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum
Please write a short summary (or long if you want): This is an example with a Zhao, Constant Beta and three component Plummer for my tutorial. Yay!

```
Think long and hard about how many walkers you want. One thousand is good number, just saying.
*Burn-in* is the number of steps for which you want to run each of your walkers before **pyGravSphere** starts outputting the chains.

*Steps* is the total steps for which you want to run your walkers. So in the example above, you will get an output of the last 2500 steps per walker.

*Integration points* sets the accuracy of your integrator: how many logarithmic bins will the integrand be divided into? 100 is fine.

Consider the amount of time you will need for the job to run. The above configuration runs in just under 1 hour on 8 cores. It will take 7-8 hours on two cores.

## Priors

Alright! Your project has been created and your jobs are ready for submission. Now if you go to:

```
/Users/Nerd/Desktop/MyWorkDir/ZhaoConst3Plum/Submissions/
```

you will find your submission scripts as well as the *priors.txt* file. This is full of default **pyGravSphere** priors. Now you might want to change those. To do so, simply edit the lower and upper priors on each parameter. If you decide you want to keep any of the parameters constant, simply edit the last column and make sure "min" and "max" values are the same. So, for example, if I want to keep the inner slope *gamma* at 1, I will edit the priors:

```
rhomin  rhomax  rhoconst                5.	10. False
rsmin   rsmax   rsconst             np.log10(r_c)   1.5 False
alphamin  alphamax  alphaconst        0.5     3. False
betamin betamax betaconst               3 	7 False
gammamin  gammamax  gammaconst       1.	1.  True
beta0min  beta0max  beta0const      -1.	1.  False
m1min   m1max  m1const            np.log10(0.5*lightpower[0])     np.log10(1.5*lightpower[0])  False
a1min   a1max  a1const        0.5*lightpower[3]	1.5*lightpower[3] False
m2min   m2max  m2const         np.log10(0.5*lightpower[1])     np.log10(1.5*lightpower[1])  False
a2min   a2max  a2const          0.5*lightpower[4]	1.5*lightpower[4] False
m3min   m3max  m3const          np.log10(0.5*lightpower[2])     np.log10(1.5*lightpower[2]) False
a3min   a3max  a3const            0.5*lightpower[5]	1.5*lightpower[5] False
mstarmin  mstarmax mstarconst       np.log10(0.75*stellar_mass)     np.log10(1.25*stellar_mass) False

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

Alternatively you can run **pyGravSphere** from your working directory. Type *pygravsphere.py* and select option 2 -- Submit jobs. Now pick if you want to run all of your galaxies, or specify which ones you'd like to run.

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


### Corner plots

For example, you might want to output a corner plot of your posteriors.


### Plot the evolution of your chains

This might be useful for visual evaluation of convergence. 
Note that you cannot use this when you've used multiprocessing, as the chains are not output in any particular order.


### Output the percentiles of your mass/density models




## Publications using GravSphere 

So you got an awesome result? Great! ***You're welcome!*** Please don't foget to cite the original **GravSphere** paper:

**https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4541R/abstract**

as well as the paper that brought you here:

**Genina et al. (2019) in prep.**

Thanks!

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details

## Future developments

The era of *Python 2* is coming to an end and we will be soon upgrading to *Python 3* as well as the latest version of **EmCee**.

In terms of the functionality, the updates will include:

* Flexible number of power law bins
* Multiple tracer populations
* Proper motions data
* Accounting for rotation

## Acknowledgments

**pyGravSphere** relies heavily on the **EmCee** code. Make sure you are familiar with how it works! Don't forget to cite:

**https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract**
