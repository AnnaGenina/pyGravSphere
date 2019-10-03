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

This is a good time to decide where you want to keep all of your ***GravSphere** input data and outputs, you must pick a working directory. To select a woking directory, please edit *workdir.txt* inside the main folder with your working directory path i.e. something like

```
/Users/Nerd/Desktop/MyWorkDir/
```

Good job!




## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
