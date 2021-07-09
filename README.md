# Systematics

## Hierarchy

- *code* : all the code (algos...).
- *data* : inputs data with the *schema.json* file and others.
- *external* : other libs that are not reachable by `git clone`
- *output* : output files from executed programs (via the *code* directory)

**Note** : all temporary files or folder (like *venv*) created by libs or program are hidden by the *.gitignore* file.

## Installation

If a virtual environment already exist with all required dependencies, it's not necessary to do these following steps.

Assume that we have a terminal opened in this root repo and the latest release of *LiteBIRD_Sim* have been download from this [GitHub Repo](https://github.com/litebird/)
into the root folder as `litebird_sim-X.Y.Z.tar.gz` where *X.Y.Z* is the version.

Enter these commands :

- Create and activate a virtual environment :
```shell
virtualenv venv
source venv/bin/activate
```

- Install common libs with *pip*

```shell
pip install numpy scipy toast-cmb mpi4py ./litebird_sim-X.Y.Z.tar.gz
```

- Build and install external libs

```shell
cd external/toast-litebird
python setup.py build
python setup.py install
```

## Large File Storage (Data section)

Unfortunatly, *git* is limited by size with 100MB and our pre-computed data (synchrotron perturbation, ...) are too big (more than 200MB). The solution is to use *git-lfs*. For more informations, see [Git LFS](https://git-lfs.github.com/) docs.

**Note** : 
- The *git lfs migrate* command is not necessary when we commit new files.
- Only *fits* files are tracked and stored as large files. That means to add new heavy file with a different extension, we must proceed *git lfs track ".\<extension\>"*.

**Warning** : When we clone this repository, lfs files are referenced and not downloaded automaticaly. To force download we must execute *git lfs pull*.

## NERSC

On a local machine, in general, there is no error when we launch a *python* program. With the *NERSC* supercomputer, somes *MPI* errors (relative to *OpenMPI*) can be thrown. However, everything is fine : we must launch throught the *Slurm* script.

**Note** : most of these scipts are already written in the *code/slurm/* folder.