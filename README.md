# Description
This repository contains the datasets and analysis code for the 3rd Year laser spectroscopy experiment.

# Datasets

Datafiles are located in the data folder. Current files and short descriprors:
- PD1_1.csv/PD2_1.csv: Oscilloscope output of PD1/PD2 with the doppler broadened spectrum.
- PD2_2.csv: Oscilloscope output of PD2 with the pump beam introduced into the system. PD2 is the probe beam. PD1's output is lost due to an error.
- PD1_3.csv/PD2_3.csv: Oscilloscope output for PD1(reference), and PD2(probe). Setup is just as above.
- PD1_4.csv/PD2_4.csv: Oscilloscope output for PD1(reference), and PD2(probe). Setup is just as above.

# Analysis code

The analysis code included in the repository and their usages are:
- `plot_rough_data.py dataset [E_scaling]`: Plots the rough data with a moving average applied to it, as well as allows fitting to it using the `fit_line(xmin,xmax)` method afterwards. The `E_scaling` can be supplied to transform Time into frequency for x axis. A secondary function called `fit_splitting(amin.amax,bmin,bmax)` is included. This can be used to get the splitting between the lines characterised by bounds `[amin,amax]` and `[bmin,bmax]`.
- `plot_hyperfine.py probe reference cull_low cull_high [E_scaling]`: Code used to remove the background of the probe beam measurement. The `cull_low` and `cull_high` parameters remove the first `cull_low` datapoints and the last `cull_high` datapoints, so they are not fitted. This is so that variation on doppler absorption depth variance is less of an issue. The `fit(xmin,xmax)` function is still present and now fits the background removed hyperfine structure. `E_scaling` can still be applied with the same effect. `fit_splitting(amin.amax,bmin,bmax)` is included with the same effect as before.
