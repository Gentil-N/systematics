# Systematics

## Hierarchy

- *code* : all the code (algos...).
- *data* : inputs data with the *schema.json* file.
- *external* : other libs that are not reachable by `git clone`
- *output* : output files from executed programs (via the *code* directory)

Note : all temporary files or folder created by libs or program are hidden by the *.gitignore* file.

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