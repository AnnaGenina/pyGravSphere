# GravSphere

**GravSphere** is a non-parametric spherical Jeans analysis code, ideal for applications to dwarf spheroidal kinematic data.

We've tried our best to design a user-friendly interface for **GravSphere** and below we include a full tutorial, however, if any questions or suggestions arise, do feel free to contact *anna.genina@durham.ac.uk*.



## Getting Started

To download **GravSphere**, please go to your preferred folder and type:

```
git install GravSphere
```


Before you begin using **GravSphere** please make sure you have downloaded and installed the following dependencies ...and their dependencies.

* gcc compiler
* numpy (duh)
* scipy
* matplotlib
* EmCee https://emcee.readthedocs.io/en/latest/user/install/
* SWIG http://www.swig.org/download.html
* mpi4py https://pypi.org/project/mpi4py/   (optional, but *STRONGLY* advised if you have some cores available)
* h5py https://pypi.org/project/h5py/


## Installing

After installing you dependencies please go to the **GravSphere** main folder and type:

```
python setup.py build
```

If you don't get any error messages, congratulations! **GravSphere** is ready to use. 


## Practical

This is a good time to decide where you want to keep all of your **GravSphere** input data and outputs. Time to pick a working directory. To select a woking directory, please update *workdir.txt* inside the main folder with your working directory path i.e. something like

```
/Users/Nerd/Desktop/MyWorkDir/
```

Good job!

### Preparing a data file

**GravSpere** currently pre-processes input data in *.hdf5* format. In just a few short steps *.hd5* files will become your new best friends! 

You may have two data sets: one for photometry (positions) of the stars and one for kinematics (positions and velocities). In the following we'll assume that these are the same two datasets.

So suppose you have a 3-column data file with x,y positions (***in parsecs***) on the sky and line-of-sight velocities. To create a **GravSphere**-compatible *.hdf5* file use the following Python code:

```
import numpy as np
import h5py

f_old = np.loadtxt("OldDataFile.txt") # x y vz
n_stars = len(f_old)

xy = f_old[:,[0,1]]
vz = f_old[:,2]
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
dset6 = output.create_dataset('KinematicsMasses', (nstars,), dtype = xy.dtype)
dset6[...] = weights
dset5 = output.create_dataset('StellarMass3R', (1,), dtype = xy.dtype)
dset5[...] = total_mass

f_new.close()

```
Obviously, change the contents of each data set as appropriate. The *Photometry* dataset will be used exclusively for calculating the best fit 3-component Plummer profile. The *Kinematics* dataset will be used for calculation of the binned velocity dispersion profile and the virial shape parameters. If the mass in stars is assumed to significantly contribute to your system, put in your best estimate of mass **in solar masses** in to *[StellarMass3R]* set. This will vary within your predefined limits in **GravSphere** MCMC run.

Now you have an *NewDataFile.hdf5* dataset. Perhaps you want to name it as something more identifiable, let's say *Galaxy_1.hdf5*, or you might even have three data sets: *Galaxy_1.hdf5*, *Galaxy_2.hdf5* and *Galaxy_3.hdf5*. At this point, please open the *galaxy_list.txt* file in your working directory and type in: 

```
Galaxy_1
Galaxy_2
Galaxy_3
```

### Pre-processing your data

Excellent! Now you're ready to to create your MCMC input files, that is binned velocity dispersion, surface brightness, half-light radius and two virial shape parameters. You might even have a different way of computing these, in which case I would advise to output your pre-processed data in the same format as we do here.

*Note that after the pre-processing the positions of the stars are now in **kiloparsecs** .*

Go to your working directory and type:

```
gravsphere.py
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
Bingo! Cogs are turning and your data will be pre-processed in no time.

### Creating a project

OK, we're definitely getting closer to the good stuff.

So suppose now you want to run a model with the Zhao et al. (1996) dark matter distribution, a constant anisotropy profile and a 3-component Plummer fit. And, very conveniently, you have 8 cores to run it on.

Now run **GravSphere** 
```
python gravsphere.py
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
Think long and hard about how many walkers you want. One thousand is good number, just sayin'.
*Burn-in* is the number of steps for which you want to run each of your walkers before **GravSphere** starts outputting the chains.
*Steps* is the total steps for which you want to run your walkers. So in the example above, you will get an output of the last 2500 steps per walker.


## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## Publications using GravSphere 

So you got an awesome result? Great! ***You're welcome!*** Please don't foget to cite the original **GravSphere** paper:

**https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4541R/abstract**

as well as the paper that brought you here:

**Genina et al. (2019) in prep.**

Thanks!

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

