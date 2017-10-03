find sub_cycles -name "*.dat" | xargs rm

./bin/cbridge4 sub_cycles/N11t1_R82mm_vitesse_deplact5mu_par_sec
./bin/cbridge4 sub_cycles/N11t1_R82mm_vitesse_deplact5mu_par_sec shift=-25e-3 p1=-5 p2=-17 p3=9

./bin/cbridge4 sub_cycles/N11t2_R82mm_vitesse_deplact5mu_par_sec
./bin/cbridge4 sub_cycles/N11t2_R82mm_vitesse_deplact5mu_par_sec  shift=-25e-3 p1=-5 p2=-17 p3=9

./bin/cbridge4 sub_cycles/N11t3_R82mm_vitesse_deplact5mu_par_sec
./bin/cbridge4 sub_cycles/N11t3_R82mm_vitesse_deplact5mu_par_sec  shift=-25e-3 p1=-5 p2=-17 p3=9

./bin/cbridge4 sub_cycles/N11t4_R82mm_vitesse_deplact5mu_par_sec
./bin/cbridge4 sub_cycles/N11t4_R82mm_vitesse_deplact5mu_par_sec  shift=-25e-3 p1=-5 p2=-17 p3=9


