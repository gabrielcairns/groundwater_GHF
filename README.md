# groundwater_GHF
Code to accompany the manuscript "Modification of Antarctic geothermal heat flux by groundwater flow" by GJ Cairns, GP Benham and IJ Hewitt. Code was written by GJ Cairns, with the variable meshgrid in 'variable_grid_300.mat' produced from code by IJ Hewitt.

This code can be used to generate all Figures apart from Figure 2 in the manuscript as follows:
make_antarctic_figures.m – Figures 1, 5, 6, 8, 11 and 12
plot_test_case.m – Figures 3, 4, 9 and 10
plot_CS_test_case.m – Figures 7

All other scripts are auxiliary functions used in the above.

NOTE: This code requires several packages and datasets in order to run, which are listed as follows:

Packages:
Antarctic Mapping Tools (https://mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools)
MEaSUREs (https://mathworks.com/matlabcentral/fileexchange/47329-measures)
Antarctic boundaries, grounding line and masks from InSAR (https://mathworks.com/matlabcentral/fileexchange/60246-antarctic-boundaries-grounding-line-and-masks-from-insar)
Hatchfill2 (https://mathworks.com/matlabcentral/fileexchange/53593-hatchfill2)
cmocean perceptually-uniform colormaps (https://mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)

Data:
"/Bedmap3/" 
Link: https://data.bas.ac.uk/full-record.php?id=GB/NERC/BAS/PDC/01615
Pritchard et al. (2025), Bedmap3 updated ice bed, surface and thickness gridded datasets for Antarctica. Scientific data, 12(1), 414.

"AIS_BaTh_v1.nc"
Link: https://zenodo.org/records/15556691
Seiner, O. et al. (2025). A synthesis of the basal thermal state of the Antarctic ice sheet. Journal of Glaciology, 1-30.

"Ant_Crust.mat"
Link: https://zenodo.org/records/10242299
Li, L., & Aitken, A. R. A. (2024). Crustal heterogeneity of Antarctica signals spatially variable radiogenic heat production. Geophysical Research Letters, 51(2), e2023GL106201.

"aq1_01_20.nc" 
Link: https://doi.pangaea.de/10.1594/PANGAEA.924857
Stål, T., et al. (2021). Antarctic geothermal heat flow model: Aq1. Geochemistry, Geophysics, Geosystems, 22(2), e2020GC009428.

"/ICESat1_ICESat2_mass_change_updated_2_2021/" 
Link: http://hdl.handle.net/1773/45388
Smith, B., et al. (2020). Pervasive ice sheet mass loss reflects competing ocean and atmosphere processes. Science, 368(6496), 1239-1242.
