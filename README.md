# Whalen_et_al_2021

Simulations were run in C++ on a Linux machine. Corresponding figures were then processed in MATLAB (2020b) and Python (3.10.6). See the Figure section below for more details.

# Figures
* Figure 1
   - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -arch 2 -dt .025 -T 50 -sd 749141 -dep_off -fac_off -gstnton 0.15 0.25 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive -o >data/competitive.hst`

* Figure 2 B-D
    - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -arch 1 -dt .025 -T 50 -sd 749141 -dep_off -fac_off -gstnton 0.15 0.25 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -pg_cross 0.00 -ps_cross 1.00 -write_syn 0 data/full_seg -o >data/full_seg.hst`
    - Follow the documentation in `Whalen2021_plot_fits.m`
* Figure 2 F-G
    - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -arch 1 -dt .025 -T 50 -sd 749141 -dep_off -fac_off -gstnton 0.15 0.25 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -pg_cross 0.25 -ps_cross 0.75 -write_syn 0 data/partial_seg -o >data/partial_seg.hst`
    - Follow the documentation in `Whalen2021_plot_fits.m`

* Figure 3 B-C
    - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -arch 2 -dt .025 -T 50 -sd 749141 -dep_off -fac_off -gstnton 0.15 0.25 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive -o >data/competitive.hst`
    - Follow the documentation in `Whalen2021_plot_fits.m`
* Figure 3 D
    - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -arch 2 -dt .025 -T 50 -sd 749141 -fac_off -dep_off -gstnton 0.225 0.375 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0.2 0.2 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 0 0 -osc_mod_gpe2 0 0 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 -control 1 data/competitive_control -o >data/competitive_control.hst`
    - Follow the documentation in `Whalen2021_plot_fits.m`
* Figure 3 E-H
    - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -arch 2 -dt .025 -T 50 -sd 749141 -dep_off -fac_off -gstnton 0.15 0.25 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive -o >data/competitive.hst`
    - Follow the documentation in `Whalen2021_plot_pairwise_phase.m`

* Figure 4
    - `$ python3 plot_freq_test_jitter.py`
    - `$ python3 analyze_runs.py`
    - Run `analyze_jitter_sims.m` in MATLAB.
    - `$ python3 plot_jitter_results.py`

* Figure 5 and Supplemental 1
    - `$ python3 run_snr_stn.py`
    - Run `analyze_multi_sims.m` in MATLAB.
    - `$ python3 plot_stn_results.py`

* Figure 6 A-B
    - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -qif 1 -qif_sigma 0 -arch 2 -dt .025 -T 50 -sd 749341 -dep_off -fac_off -gstnton 0.001 0.001 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .0015 .0075 -wg .0015 .0075 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive_qif -o >data/competitive_qif.hst`
    - Follow the documentation in `Whalen2021_plot_fits.m`
* Figure 6C
    - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -qif 1 -qif_sigma 0.015 -arch 2 -dt .025 -T 50 -sd 749341 -dep_off -fac_off -gstnton 0.001 0.001 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .0015 .0075 -wg .0015 .0075 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive_qif_noise015 -o >data/competitive_qif_noise015.hst`
    - `$ ./SNr -qif 1 -qif_sigma 0.03 -arch 2 -dt .025 -T 50 -sd 749341 -dep_off -fac_off -gstnton 0.00095 0.00095 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .0015 .0075 -wg .0015 .0075 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive_qif_noise030 -o >data/competitive_qif_noise030.hst`
    - `$ ./SNr -qif 1 -qif_sigma 0.045 -arch 2 -dt .025 -T 50 -sd 749341 -dep_off -fac_off -gstnton 0.00085 0.00085 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .0015 .0075 -wg .0015 .0075 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive_qif_noise045 -o >data/competitive_qif_noise045.hst`
    - `$ ./SNr -qif 1 -qif_sigma 0.06 -arch 2 -dt .025 -T 50 -sd 749341 -dep_off -fac_off -gstnton 0.0008 0.0008 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .0015 .0075 -wg .0015 .0075 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive_qif_noise060 -o >data/competitive_qif_noise060.hst`
    - `$ ./SNr -qif 1 -qif_sigma 0.075 -arch 2 -dt .025 -T 50 -sd 749341 -dep_off -fac_off -gstnton 0.0006 0.0006 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .0015 .0075 -wg .0015 .0075 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive_qif_noise075 -o >data/competitive_qif_noise075.hst`
    - `$ ./SNr -qif 1 -qif_sigma 0.09 -arch 2 -dt .025 -T 50 -sd 749341 -dep_off -fac_off -gstnton 0.0004 0.0004 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .0015 .0075 -wg .0015 .0075 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive_qif_noise090 -o >data/competitive_qif_noise090.hst`
    - Follow the documentation in `Whalen2021_plot_qif_noise.m`

* Figure 7
    - `$ g++ -std=c++17 -O2 SNr.cc -o SNr`
    - `$ ./SNr -arch 2 -dt .025 -T 50 -sd 749141 -dep_off -fac_off -gstnton 0.15 0.25 -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 -write_syn 0 data/competitive -o >data/competitive.hst`
    - Follow the documentation in `Whalen2021_plot_balance.m`

* Supplementary Figure 2
    - Jitter with synaptic competition