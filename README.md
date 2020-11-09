# Description
This repository contains the datasets and analysis code for the 3rd Year laser spectroscopy experiment.

# Datasets

Datafiles are located in the data folder. Current files and short descriprors:
- PD1_1.csv/PD2_1.csv: Oscilloscope output of PD1/PD2 with the doppler broadened spectrum.
- PD2_2.csv: Oscilloscope output of PD2 with the pump beam introduced into the system. PD2 is the probe beam. PD1's output is lost due to an error.
- PD1_3.csv/PD2_3.csv: Oscilloscope output for PD1(reference), and PD2(probe). Setup is just as above.
- PD1_4.csv/PD2_4.csv: Oscilloscope output for PD1(reference), and PD2(probe). Setup is just as above.
- PD1_5.csv/PD2_5.csv/PD3_5.csv: Oscilloscope output for PD1(reference), PD2(probe) and PD3(etalon). Note that this dataset has a negative time-energy dependence.
- PD3_S5_S/E/F.CSV: The etalon measurements at the start middle and end of the session.
- PD1_6.csv/PD2_6.csv: Datasets taken as usual, with etalon blocked. Taken right after PD3_S5_S.
- PD1_7.csv/PD2_7.csv: Datasets taken as usual, with etalon blocked. Taken right after PD3_S5_E.
- PD1_8.csv/PD2_8.csv: Datasets taken as usual, with etalon blocked. Taken right after PD3_S5_F.
- PD2_8_B.csv: PD2 data with the pump beam blocked. Taken after PD1/2_8.
- PD1_9/PD2_9/PD2_9_B/PD3_9: Same naming convention as above. New etalon alignment.

# Analysis code

The analysis code included in the repository and their usages are:
- `plot_rough_data.py dataset [E_scaling]`: Plots the rough data with a moving average applied to it, as well as allows fitting to it using the `fit_line(xmin,xmax)` method afterwards. The `E_scaling` can be supplied to transform Time into frequency for x axis. A secondary function called `fit_splitting(amin.amax,bmin,bmax)` is included. This can be used to get the splitting between the lines characterised by bounds `[amin,amax]` and `[bmin,bmax]`.
- `plot_rough_data_etalon.py dataset [etalon]`: Same as `plot_rough_data.py`, but instead of `E_scaling` the etalon dataset can be provided. The etalon's data is plotted and found peaks are marked as a sanity check.
- `plot_hyperfine.py probe reference cull_low cull_high [E_scaling]`: Code used to remove the background of the probe beam measurement. The `cull_low` and `cull_high` parameters remove the first `cull_low` datapoints and the last `cull_high` datapoints, so they are not fitted. This is so that variation on doppler absorption depth variance is less of an issue. The `fit(xmin,xmax)` function is still present and now fits the background removed hyperfine structure. `E_scaling` can still be applied with the same effect. `fit_splitting(amin.amax,bmin,bmax)` is included with the same effect as before. `fit_structure(xmin,xmax,m1=0,m2=0,m3=0,p0=None)` is included to provide an optimisation of the structure. It fits the true and cross-over peaks. The initial values are either the true peak's rough location (m1-3) or, the optimisation parameter's initial guess (p0).
- `plot_hyperfine_etalon.py probe reference cull_low cull_high etalon`: Same as `plot_hyperfine.py`, but instead of `E_scaling` the etalon dataset can be provided. The etalon's data is plotted and found peaks are marked as a sanity check. Unlike in `plot_rough_etalon.py` here the etalon peaks are interpolated to give an accurate scaling. This also has an extra method `fit_total(amin,amax,bmin,bmax,am1,am2,am3,bm1,bm2,bm3)`, which fits a paired hyperfine structure and determines the S1/2 spliting as well. Note that while the parameters are the same as for `fit_structure`, it needs that `am1<am2<am3<bm1<bm2<bm3` to work. This file has these methods and general plotting improved.

# Setup

A schematic view of the setup used.
!(Setup)[Plots/Setup.png]
