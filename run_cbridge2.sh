find sub_cycles -name "*.dat" | xargs rm

PERCENT_ENF=-5
PERCENT_PLA=-17
PERCENT_TIR=-9

./bin/cbridge4 "file='sub_cycles/N11t1_R82mm_vitesse_deplact5mu_par_sec/enfoncement_brute.txt'" main_rate=5e-3 percent=0
./bin/cbridge4 "file='sub_cycles/N11t1_R82mm_vitesse_deplact5mu_par_sec/enfoncement_brute.txt'" main_rate=5e-3 percent=$PERCENT_ENF
./bin/cbridge4 "file='sub_cycles/N11t1_R82mm_vitesse_deplact5mu_par_sec/plateau_brute.txt'" main_rate=5e-3     percent=0
./bin/cbridge4 "file='sub_cycles/N11t1_R82mm_vitesse_deplact5mu_par_sec/plateau_brute.txt'" main_rate=5e-3     percent=$PERCENT_PLA
./bin/cbridge4 "file='sub_cycles/N11t1_R82mm_vitesse_deplact5mu_par_sec/tirage_brute.txt'" main_rate=5e-3      percent=0
./bin/cbridge4 "file='sub_cycles/N11t1_R82mm_vitesse_deplact5mu_par_sec/tirage_brute.txt'" main_rate=5e-3      percent=$PERCENT_TIR

./bin/cbridge4 "file='sub_cycles/N11t2_R82mm_vitesse_deplact5mu_par_sec/enfoncement_brute.txt'" main_rate=5e-3 percent=0
./bin/cbridge4 "file='sub_cycles/N11t2_R82mm_vitesse_deplact5mu_par_sec/enfoncement_brute.txt'" main_rate=5e-3 percent=$PERCENT_ENF
./bin/cbridge4 "file='sub_cycles/N11t2_R82mm_vitesse_deplact5mu_par_sec/plateau_brute.txt'" main_rate=5e-3     percent=0
./bin/cbridge4 "file='sub_cycles/N11t2_R82mm_vitesse_deplact5mu_par_sec/plateau_brute.txt'" main_rate=5e-3     percent=$PERCENT_PLA
./bin/cbridge4 "file='sub_cycles/N11t2_R82mm_vitesse_deplact5mu_par_sec/tirage_brute.txt'" main_rate=5e-3      percent=0
./bin/cbridge4 "file='sub_cycles/N11t2_R82mm_vitesse_deplact5mu_par_sec/tirage_brute.txt'" main_rate=5e-3      percent=$PERCENT_TIR


./bin/cbridge4 "file='sub_cycles/N11t3_R82mm_vitesse_deplact5mu_par_sec/enfoncement_brute.txt'" main_rate=5e-3 percent=0
./bin/cbridge4 "file='sub_cycles/N11t3_R82mm_vitesse_deplact5mu_par_sec/enfoncement_brute.txt'" main_rate=5e-3 percent=$PERCENT_ENF
./bin/cbridge4 "file='sub_cycles/N11t3_R82mm_vitesse_deplact5mu_par_sec/plateau_brute.txt'" main_rate=5e-3     percent=0
./bin/cbridge4 "file='sub_cycles/N11t3_R82mm_vitesse_deplact5mu_par_sec/plateau_brute.txt'" main_rate=5e-3     percent=$PERCENT_PLA
./bin/cbridge4 "file='sub_cycles/N11t3_R82mm_vitesse_deplact5mu_par_sec/tirage_brute.txt'" main_rate=5e-3      percent=0
./bin/cbridge4 "file='sub_cycles/N11t3_R82mm_vitesse_deplact5mu_par_sec/tirage_brute.txt'" main_rate=5e-3      percent=$PERCENT_TIR

./bin/cbridge4 "file='sub_cycles/N11t4_R82mm_vitesse_deplact5mu_par_sec/enfoncement_brute.txt'" main_rate=5e-3 percent=0
./bin/cbridge4 "file='sub_cycles/N11t4_R82mm_vitesse_deplact5mu_par_sec/enfoncement_brute.txt'" main_rate=5e-3 percent=$PERCENT_ENF
./bin/cbridge4 "file='sub_cycles/N11t4_R82mm_vitesse_deplact5mu_par_sec/plateau_brute.txt'" main_rate=5e-3     percent=0
./bin/cbridge4 "file='sub_cycles/N11t4_R82mm_vitesse_deplact5mu_par_sec/plateau_brute.txt'" main_rate=5e-3     percent=$PERCENT_PLA
./bin/cbridge4 "file='sub_cycles/N11t4_R82mm_vitesse_deplact5mu_par_sec/tirage_brute.txt'" main_rate=5e-3      percent=0
./bin/cbridge4 "file='sub_cycles/N11t4_R82mm_vitesse_deplact5mu_par_sec/tirage_brute.txt'" main_rate=5e-3      percent=$PERCENT_TIR
