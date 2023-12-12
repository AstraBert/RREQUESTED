#!/bin/bash

echo "Installation process started, checking if Anaconda/Miniconda and Mamba are installed..."


command -v conda >/dev/null 2>&1 || {
   echo >&2 "The installation pipeline requires Anaconda/Miniconda but it is not installed. Please check here: https://anaconda.org/ for more details. Aborting."
   exit 1
}


command -v mamba >/dev/null 2>&1 || {
   echo >&2 "The installation pipeline requires Mamba but it is not installed. Please check here: https://github.com/conda-forge/miniforge#mambaforge for more details. Aborting."
   exit 2
}

echo "Anaconda/Miniconda and Mamba are correctly installed, proceeding with the installation of necessary softwares: it may take a while..."

sleep 1

tmp=$(dirname $0)
basedir=$(realpath $tmp)

ConPath=$(which conda)
tmp=${ConPath#* }
Conda=${tmp%%/bin/co*}

# CREATE A NEW FOLDER WHERE TO STORE THE ENVIRONMENT
mkdir -p $basedir/environments

# CREATE A MAMBA/CONDA ENVIRONMENT
mamba create \
  -p $basedir/environments/rrequested  \
  -y \
  -c conda-forge \
  -c bioconda \
  python=3.10

source ${Conda}/etc/profile.d/conda.sh
conda activate \
  ${basedir}/environments/rrequested

## INSTALL necessary dependencies
mamba install \
   -c anaconda \
   -y \
   numpy=1.26.2

mamba install \
   -c bioconda \
   -y \
   edlib pandas

mamba install \
   -c conda-forge \
   -y \
   pigz

## DEACTIVATE
conda deactivate

# MAKE RREQUESTED.sh AN EXECUTABLE
echo "alias RREQ='bash ${basedir}/RREQUESTED.sh'">>~/.bash_aliases


