#!/bin/sh

make clean
make main_profiled

pa="/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/000_Trx/"
./main_profiled $pa/000_Trx_pot_positive.pts $pa/000_Trx_mu_positive.pts $pa/000_Trx_pot_negative.pts $pa/000_Trx_mu_negative.pts $pa/000_Trx_active_center.pts $pa/000_Trx_pot_positive.pts $pa/000_Trx_mu_positive.pts $pa/000_Trx_pot_negative.pts $pa/000_Trx_mu_negative.pts $pa/000_Trx_active_center.pts 100 100000 1 0 0 0 $pa/out.txt


gprof main_profiled gmon.out > profile-data.txt
gedit profile-data.txt
