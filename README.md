# RIkin

Rikin implements tools for the fast computation of detailed RNA-RNA interaction kinetics
in a detailed RNAup/Intarna-inspired interaction model.

It features accurate modeling of RNA structures and interaction complexes based on their secondary structure and the full Turner nearest neighbor energy model.

## Example


## Installation

### Installation from Conda package

We recommend to install from the conda package of rikin, using mamba.

Possible create and activate an environment for rikin first:
```
conda create rikin
conda activate rikin
```

Then, install by

```mamba install -c bioconda rikin```

or alternatively conda

```conda install -c bioconda rikin```

Then, activate the environment 

```
conda activate rikin
```

and install required R pacakges

```
Rscript -e "install.packages(c('argparse','shape','RColorBrewer'))"
```


### Compilation/installation from the source repository

The tools can be compiled and installed after cloning the source repository. This requires a build toolchain with C++ compiler and autotools. We describe the installation in a conda environment (and get further specific dependencies from bioconda).

```mamba create -n rikin -c bioconda -c conda-forge viennarna locarna r-base

conda activate rikin

autoreconf -i
./configure --prefix=$CONDA_PREFIX PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
make
make install

```

Finally, see above (installation from conda package) how to install the required R packages.


## Usage


