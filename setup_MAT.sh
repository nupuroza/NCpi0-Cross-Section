## Set up UPS
source /opt/root/6.32.00_binary/root/bin/thisroot.sh

## Right now the installation of the MAT that has been correctly modified for uBooNE
## usage lives exclusively in Rob's user area. He plans to commit the needed change
## to a branch of the MAT repo so that it can be checked out by any future uBooNE users.
## Source code for PlotUtils at /exp/uboone/app/users/ltong/MAT/DPUtils/PlotUtils/src

## Add MAT to relevant environment variables
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/users/lntong/MAT/opt/lib/
export PATH=$PATH:/app/users/lntong/MAT/opt/bin/
export ROOT_INCLUDE_PATH=/app/users/lntong/MAT/opt/include/PlotUtils:/app/users/lntong/MAT/opt/include:${ROOT_INCLUDE_PATH}

