#!/bin/bash

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh

setup git v2_15_1;
setup gitflow v1_11_0;
setup mrb v1_16_02;
setup root v6_12_04e -q e15:prof
setup cmake v3_11_4;
setup eigen v3_3_4a;

./sbnfit_make_specN -x ../xmls/nuwro_ncpi0_Aug2023.xml -t NuWro_Aug2023_v8_Exclusive
./sbnfit_make_specN -x ../xmls/nuwro_ncpi0_Jun2024_Inclusive.xml -t NuWro_Jun2024_v9_Inclusive

