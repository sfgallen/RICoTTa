# RICoTTa
 
### ***R***iver ***I***nversion Model using ***Co***smogenic Radionuclide, Marine ***T***errace, and Low-Temperature ***T***hermochronological  Dat***a***

[![DOI](https://zenodo.org/badge/601847671.svg)](https://zenodo.org/badge/latestdoi/601847671)

This series of codes takes a suite of data sets (i.e., longitudinal river profiles, cosmogenic radionuclide concentrations, marine terrace rock uplift calculations, and low-temperature thermochronometric ages) and uses a Bayesian inversion to constrain rock uplift histories and parameters in the detachment-limited stream power model.

These codes are modified versions of the algorithm described in Glotzbach (2015). The version of the codes in this repository was applied in Gallen et al. (2023). Please see those references for details on the model.

The three example scripts provided will run an abbreviated simulation using the datasets in Gallen et al. (2023). Note that the data provided with the examples was prepared using SRTM digital elevation models and TopoToolbox v2 (Schwanghart and Scherler, 2014) and ChiProfiler (Gallen and Wegmann, 2017).

## Disclaimer:

We _strongly_ recommend that only expert Matlab users deeply familiar with the forward and inverse models used in this software modify the codes provided. We make this disclaimer because it is fairly easy to make mistakes that will inhibit the usefulness of the model results.

### Furthermore, before attempting to apply these codes to another study area, it is important to make sure that a number of assumptions in the natural setting are met. These assumptions include, but are not limited to:

-	The catchment is an erosional, fluvially-dominated, bedrock system (i.e., no deposition, no evidence of past glaciation, and no alluvial rivers).
-	The catchment network has no obvious evidence of recent reorganization or unstable drainage divides (e.g., no evidence of recent river capture or strongly asymmetric drainage divides).
-	Erodibility is relatively uniform throughout the catchment (i.e., roughly uniform rock type, no strong precipitation gradients).
-	Rock uplift appears spatially uniform (i.e., no evidence of local or regional tilting and no active faults cutting the catchment).
-	The uplift and exhumation history can reasonably be approximated as one-dimentionsal in the vertical direction (i.e., no significantly horizontal tectonic velocities).
-	The catchments are not too large (no firm size threshold, but the larger the basin, the more likely the assumptions will be violated).
-	The outlet location is fixed with time.

Note that this is not a comprehensive list, but these are some of the assumptions that must be met for this code to be meaningfully applied to a given catchment.

For an example of working through a list of assumptions in conducting similar catchment-scale inversions, please see Gallen and Fernández-Blanco (2021).

### If you use these codes of modified versions of these codes for scientific research, cite Glotzbach (2015) and Gallen et al. (2023).

### References:

Gallen, S.F., Seymour, N.M., Glotzbach, C., Stockli, D.F., O'Sullivan, P., _accepted_, Calabrian forearc uplift paced by slab-mantle interactions during subduction retreat: Nature Geoscience.

Gallen, S.F. and Fernández-Blanco, D., 2021, A New Data-driven Bayesian Inversion of Fluvial Topography Clarifies the Tectonic History of the Corinth Rift and Reveals a Channel Steepness Threshold: JGR-Earth Surface, v. 126, p. e2020JF005651. https://doi.org/10.1029/2020JF005651.

Gallen, S.F. and Wegmann, K.W., 2017, River profile response to fault growth and linkage: An example from the Hellenic forearc, south-central Crete, Greece: Earth Surface Dynamics, v. 5, p. 161-186. https://doi.org/10.5194/esurf-5-161-2017/

Schwanghart, W. and Scherler, D., 2014. TopoToolbox 2–MATLAB-based software for topographic analysis and modeling in Earth surface sciences. Earth Surface Dynamics, 2(1), pp.1-7. https://doi.org/10.5194/esurf-2-1-2014

Glotzbach, C., 2015. Deriving rock uplift histories from data-driven inversion of river profiles: Geology, 43(6), pp.467-470. https://doi.org/10.1130/G36702.1

