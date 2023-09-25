To use it as an add-on module, follow these steps:

1) Ensure all code paths are updated to reflect the correct file locations in the Filepath (e.g., Filepath.GCM_AMIPW), like this:

../GRIST_increm/amipw_physics/ozone/

../src/atmosphere/gcm/
../src/atmosphere/gcm/io_amipw_template/
../src/atmosphere/dynamics/dycore/
../src/atmosphere/dynamics/dycore/hpe/
../src/atmosphere/dynamics/dycore/nhd/
......

2) Compile the model code with the additional “-DCLIMATE_O3” flag to enable the module as an add-on.
