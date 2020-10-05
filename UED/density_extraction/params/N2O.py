import os
import numpy as np

data_parameters = {
    "molecule"           : "N2O",
    "experiment"         : "UED",
    "simulate_data"      : True,
    "simulate_data_scale": 1e2/14,
    "scale"              : 1.,
    "data_fileName"      : os.path.join(
          "/reg/ued/ana/scratch/N2O/timeBases",
          "UED_fit_results_Temp-30.0-150.0-121_Ints-1.0-10.0-91.h5"),
    #"perturb_range"     : np.array([0.01, 0.01, 0.01]),  # angs
    #"perturb_range"     : np.array([0.02, 0.02, 0.01]),  # angs
    "perturb_range"     : np.array([0.0075, 0.0075, 0.0075]),  # angs
    #"data_projections"  : np.arange(6)*2,
    "fit_bases"         : np.array([2]),
    #"bases_address"     : "bases/UED_best_fit_Nbases-2/",
    "isMS"              : True,#False,
    #"binning"           : "log",
    #"min_q"             : 0.0,
    #"max_q"             : 14.12,
    "fit_range"         : [3,11.5],#[2.5, 8],
    "elEnergy"          : 3.7e6,
    "init_geo_xyz"      : "XYZ/N2O.xyz", 
    "scat_amps_path"    : "/reg/neh/home5/khegazy/baseTools/simulation/scatteringAmplitudes/3.7MeV/"
}
