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

# argweaver
# module load python/2.7.15
# module load gcc
# cd programs
# git clone https://github.com/CshlSiepelLab/argweaver.git
# cd argweaver
# make
# make install prefix=$(pwd)/local
# export PATH=$PATH:$(pwd)/local
# export PYTHONPATH=$PATH:$(pwd)/local/lib/python2.7/site-packages
# now you can run arg-sim and arg-sample examples

# to use smc2bed-all you need to get samtools and bedops
# module load samtools
# wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
# tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
# mv bin bedops
# export PATH=/home/m/mmosmond/mmosmond/scratch/projects/tsrescue/programs/bedops:$PATH
# now you can run smc2bed on the examples 

# tsconvert to get ts from newicks
# git clone https://github.com/tskit-dev/tsconvert.git programs/tsconvert
# cd programs/tsconvert
# python setup.py install

# git clone https://github.com/leospeidel/relate_lib.git programs/relate_lib
# cd programs/relate_lib
# mkdir build
# cd build
# module load cmake gcc
# cmake ..
# make
