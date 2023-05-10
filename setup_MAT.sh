## Set up UPS
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh

## Set up compatible versions of ROOT and cmake
setup root v6_12_06a -q e15:prof
setup cmake v3_7_0

## Right now the installation of the MAT that has been correctly modified for uBooNE
## usage lives exclusively in Rob's user area. He plans to commit the needed change
## to a branch of the MAT repo so that it can be checked out by any future uBooNE users

## Add MAT to relevant environment variables
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/uboone/app/users/finer/MAT/opt/lib/
export PATH=$PATH:/uboone/app/users/finer/MAT/opt/bin/
export ROOT_INCLUDE_PATH=/uboone/app/users/finer/MAT/opt/include/PlotUtils:/uboone/app/users/finer/MAT/opt/include:${ROOT_INCLUDE_PATH}

