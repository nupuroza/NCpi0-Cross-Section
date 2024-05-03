## Set up UPS
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh

## Set up compatible versions of ROOT and cmake
setup root v6_12_06a -q e15:prof
setup cmake v3_7_0

## Right now the installation of the MAT that has been correctly modified for uBooNE
## usage lives exclusively in Rob's user area. He plans to commit the needed change
## to a branch of the MAT repo so that it can be checked out by any future uBooNE users.
## Source code for PlotUtils at /exp/uboone/app/users/ltong/MAT/DPUtils/PlotUtils/src

## Add MAT to relevant environment variables
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/exp/uboone/app/users/ltong/MAT/opt2/lib/
export PATH=$PATH:/exp/uboone/app/users/ltong/MAT/opt2/bin/
export ROOT_INCLUDE_PATH=/exp/uboone/app/users/ltong/MAT/opt2/include/PlotUtils:/exp/uboone/app/users/ltong/MAT/opt2/include:${ROOT_INCLUDE_PATH}

