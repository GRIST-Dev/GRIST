To use this MPS as an add-on module for GRIST_kernel, follow these steps:

1) Ensure all code paths are updated to reflect the correct file locations in the Filepath (e.g., Filepath.GCM_AMIPW), like this:

../GRIST_increm/mlPhy/src/atmosphere/physics/amipw_physics/
../GRIST_increm/mlPhy/src/infrastructure/namelist/
../GRIST_increm/mlPhy/src/infrastructure/io/

../src/atmosphere/gcm/
../src/atmosphere/gcm/io_amipw_template/
../src/atmosphere/dynamics/dycore/
../src/atmosphere/dynamics/dycore/hpe/
../src/atmosphere/dynamics/dycore/nhd/
......


2) Compile the model code with the additional “-DResNet” flag to enable the module as an add-on.


