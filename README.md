# Emissivity optimization for a micropatterned SiO<sub>2</sub> layer

This repository contains MATLAB/Octave scripts to perform rigorous coupled-wave analysis (RCWA) simulations for a SiO<sub>2</sub> layer decorated with micropillars.

Running the code requires the RETICOLO v9 RCWA code:
> Jean-Paul Hugonin, & Philippe Lalanne. (2021). Light-in-complex-nanostructures/RETICOLO: V9. Zenodo. [https://doi.org/10.5281/zenodo.4419063](https://doi.org/10.5281/zenodo.4419063)

The repository contains:
* Exhaustive_parameter_search.m: main script setting up and running all simulations
* Input_param.csv: the explored parameter space
* SiO2_J. Kischkat, S. Peters.txt: complex refractive index of SiO2 taken from the refractiveindex.info [database](https://refractiveindex.info/?shelf=main&book=SiO2&page=Kischkat)
* Abs, Rup, Tup csv files containing the output of the exhausive search
* Exhaustive_parameter_search_plot_results.m: script interpolating and plotting the output. Running this script should generate the following figure


