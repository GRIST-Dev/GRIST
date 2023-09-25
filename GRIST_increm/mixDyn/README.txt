The mixDyn code is a specific implementation of the mixed-precision dynamics module in GRIST, as described in the following paper:

Chen, S., Y. Zhang, Y. Wang, Z. Liu, X. Li, and W. Xue, (2024), Mixed-precision computing in the GRIST dynamical core for weather and climate modelling. Geosci. Model Dev., 17(16), 6301-6318.doi:10.5194/gmd-17-6301-2024.

The add-on module remains consistent with the description in the paper but has been updated to conform with the latest kernel version. It produces bit-identical results when configured with the DBL setup.

To use it as an add-on module, follow these steps:

1) Ensure all code paths are updated to reflect the correct file locations in the Filepath (e.g., Filepath.GCM_AMIPW), like this:

../mixDyn/src/atmosphere/gcm/
../mixDyn/src/atmosphere/dtp/
../mixDyn/src/atmosphere/dynamics/dycore/
../mixDyn/src/atmosphere/dynamics/tracer_transport/
../mixDyn/src/atmosphere/dynamics_mixed/dycore/
../mixDyn/src/atmosphere/dynamics_mixed/dycore/hpe/
../mixDyn/src/atmosphere/dynamics_mixed/dycore/nhd/
../mixDyn/src/atmosphere/dynamics_mixed/tracer_transport/

../src/atmosphere/gcm/
../src/atmosphere/gcm/io_amipw_template/
../src/atmosphere/dynamics/dycore/
../src/atmosphere/dynamics/dycore/hpe/
../src/atmosphere/dynamics/dycore/nhd/
......


2) Compile the model code with the additional “-DMIXCODE” flag to enable the module as an add-on.
Note: Using the “-DMIXCHECK -DREG_KESI" flag compiles the MIXCODE in a DBL mode, producing bit-identical results with the kernel’s DBL code. This serves as a regression test to verify that the code behavior remains unchanged.
