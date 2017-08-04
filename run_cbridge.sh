MAIN_RATE=5e-3
./bin/cbridge4 "file='new_data/N11t1_R82mm_vitesse_deplact5mu_par_sec/cycle_brute.txt'" main_rate=$MAIN_RATE evap_rate=50e-6 coef_push=-0.05 coef_pull=-0.25
./bin/cbridge4 "file='new_data/N11t2_R82mm_vitesse_deplact5mu_par_sec/cycle_brute.txt'" main_rate=$MAIN_RATE evap_rate=50e-6 coef_push=-0.05 coef_pull=-0.25
./bin/cbridge4 "file='new_data/N11t3_R82mm_vitesse_deplact5mu_par_sec/cycle_brute.txt'" main_rate=$MAIN_RATE evap_rate=50e-6 coef_push=-0.15 coef_pull=-0.30
./bin/cbridge4 "file='new_data/N11t4_R82mm_vitesse_deplact5mu_par_sec/cycle_brute.txt'" main_rate=$MAIN_RATE evap_rate=50e-6 coef_push=-0.05 coef_pull=-0.20
./bin/cbridge4 "file='new_data/N13t_R82mm_vitesse_deplact5mu_par_sec/cycle_brute.txt'"  main_rate=$MAIN_RATE evap_rate=50e-6 coef_push=-0.05 coef_pull=-0.25


