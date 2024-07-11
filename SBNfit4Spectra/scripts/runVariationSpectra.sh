#!/bin/bash

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup git v2_15_1;
setup gitflow v1_11_0;
setup mrb v1_16_02;
setup root v6_12_04e -q e15:prof
setup cmake v3_11_4;
setup eigen v3_3_4a;

XMLXpINC="../xmls/master_ncpi0_2gXp_inclusive_2023_samepot.xml"
XML0pEXC="../xmls/master_ncpi0_2g0p_exclusive_2023_samepot.xml"
XML1pEXC="../xmls/master_ncpi0_2g1p_exclusive_2023_samepot.xml"

TAG="v11_d25_04_24_CombInc"

./sbnfit_make_covarianceN -x $XMLXpINC -t "INC_2gXp_"$TAG -w NCPi0sig &> log.INC.$TAG".2gXp"
./sbnfit_make_covarianceN -x $XML0pEXC -t "EXC_2g0p_"$TAG -w NCPi0sig  &> log.EXC.$TAG".2g0p"
./sbnfit_make_covarianceN -x $XML1pEXC -t "EXC_2g1p_"$TAG -w NCPi0sig  &> log.EXC.$TAG".2g1p"

