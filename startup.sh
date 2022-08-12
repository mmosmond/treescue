#!/bin/sh
#
# description: get everything set up after logging on
#

# load python
module load NiaEnv/2019b python/3.9.8

# make directory for virtual environments if havent already
# mkdir -p ~/.virtualenvs

# create virtual environment, if havent already
myenv="tsrescue"
# virtualenv --system-site-packages ~/.virtualenvs/$myenv

# activate virtual environment
source ~/.virtualenvs/$myenv/bin/activate 

# install snakemake in the environment if havent already
# pip install snakemake==7.8.3

#other modules
#pip install numpy
#pip install pyslim==1.0b1 #because current release has errors recapitating: https://github.com/tskit-dev/pyslim/issues/271

# for jupyter hub
#pip install ipykernel
#python -m ipykernel install --name $myenv --user
#venv2jup

# i think we might need this to prevent snakemake trying to write to a read-only file (https://github.com/snakemake/snakemake/issues/1593)
export XDG_CACHE_HOME=~/scratch

# slim
# wget https://github.com/MesserLab/SLiM/releases/download/v3.7.1/SLiM.zip -P programs/
# cd programs
# unzip SLiM.zip
# rm SLiM.zip
# module load cmake/3.21.4
# module load gcc/8.3.0
# mkdir build
# cd build
# cmake ../SLiM
# make slim
# cd ..
# mv build SLiM_3.7.1

# make DATADIR
mkdir -p data

