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
- `plot_rough_data.py dataset`: Plots the rough data with a moving average applied to it, as well as allows fitting to it using the `fit_line(xmin,xmax)` method afterwards.
