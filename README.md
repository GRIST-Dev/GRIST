# Global-Regional Integrated Forecast System (GRIST)

[**Online Tutorial**](https://grist-tutorial.readthedocs.io/en/latest/)  

- **Model Version:** 1.2  
- **Latest Code Version:** A25.1.6, Fortran

---

## ✨ Features

-  **Global-Regional Integrated Dynamical Core**: The concept of Global-Regional Integration can be represented in various ways, including global uniform high-resolution (km-scale) models, global models with localized high-resolution areas, and global and regional models developed within the same codebase, such that the regional model is a truncated version of the global model. All of these approaches are pursued by GRIST in an integrated manner, optimizing resource input.
-  **Two Integrated Physics Suites**:  
  - **PhysW / AMIPW_Physics**  
  - **PhysC / AMIPC_Physics**  
-  **Simplified Models** for testing physics and numerical operators.
-  **Extensive Testing** across multiple configurations:  
   Single-column mode  
   Shallow-water mode  
   3D dry/moist idealized dynamics  
   AMIP, NWP-style, and GSRM  
-  **Custom Diagnostic Templates** for various applications *(see the [GRIST tutorial](https://grist-tutorial.readthedocs.io/en/latest/) for details)*.
-  **Cost effective computational performance**: week-scale SDPD for GSRM-scale applications on normal CPU machines.
---

## General Usage

1. **Modify Environment Setup**  
   - Edit `build.sh` in the "build" directory on your local machine.  
   - Choose a model type (e.g., `amipw` or `amipc`).  
   - Run:  
     ```bash
     ./build.sh amipw  # or ./build.sh amipc
     ```
   
2. **Compilation**  
   - After successful compilation, the following files are generated:  
     - `ParGRIST_amipw.exe / ParGRIST_amipc.exe` (Model executable)  
     - `partition.exe`  
   - Place the model executable in the **"run"** directory for each case.  

3. **Running the Model**  
   - Prepare **namelist files** and necessary **input data**.  
   - Ensure the required **runtime environment** is set up.  
   - Submit the job for execution.  

---

## GRIST Kernel

The **GRIST_kernel** is the **minimum codebase** supporting all GRIST v1.0 functions, including results from published papers *(bit regression is ensured for identical setups and environments)*.  

- **GRIST_kernel** can be **easily extended** with **GRIST_increm** (add-on modules) with minimal effort.  

---

## 📌 Quick Info

A comprehensive description of the development and evaluation of the dynamical core framework is given in [Zhang et al. (2019)](https://doi.org/10.1029/2018MS001539), [(2020)](https://doi.org/10.1175/MWR-D-19-0305.1) and [Zhang et al. (2024)](https://doi.org/10.1002/qj.4804). The two baseline physics suites, PhysW and PhysC, specifically tailored for GRIST, are described and evaluated based on single column modeling ([Li et al. 2023](https://doi.org/10.5194/gmd-16-2975-2023)).   


The full-model studies involving GRIST-PhysW, GRIST-PhysC, and their AMIP simulations have been discussed in [Zhang et al. (2021)](https://doi.org/10.1029/2021MS002592) and [Li et al. (2022)](https://doi.org/10.1029/2021JD036069), focusing on model performance and climate simulation analysis. [Fu et al. (2024)](https://doi.org/10.1007/s00382-024-07205-2) [(2025)](https://doi.org/10.1007/s00382-024-07527-1) further contributed to the intercomparison of the long-term state of two model climates and the simulation of the Madden-Julian Oscillation (MJO).  

The model has found applications across several domains, including global variable-resolution modeling, global storm-resolving modeling, and climate change modeling (e.g., [Zhou et al. 2020](https://doi.org/10.5194/gmd-13-6325-2020); [Zhang et al. 2022](https://doi.org/10.1029/2022EA002401); [Sun et al. 2024](https://doi.org/10.1016/j.scib.2023.11.013)). A Chinese-language brief introduction to the model framework  was provided by [Wang et al. (2024)](http://www.cmalibrary.cn/amst/2024/202404/fmbd/202409/t20240927_165283.htm). Some background of the model development research can be found in [Yu et al. (2019)](https://doi.org/10.1007/s00376-019-8203-1) and [Zhang et al. (2023)](https://doi.org/10.1007/978-3-031-40567-9_1). More published references can be found at [**Online Tutorial**](https://grist-tutorial.readthedocs.io/en/latest/references.html).  

---
