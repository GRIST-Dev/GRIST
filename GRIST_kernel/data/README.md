Global-to-Regional Integrated forecast SysTem (GRIST)
================================================================================
GRIST is a forecast model designed for flexible and general-purpose usages in global weather and climate modeling research and applications. It is developed in a hierarchical manner, covering spectra of experimental protocols for global atmospheric modeling. It has been used for experimental earth weather and climate modeling applications.  

Table of Contents 
--------------------------------------------------------------------------------
- [Data information](#data_information)
- [Horizontal interpolation](#horizontal_interpolation)
- [Vertical interpolation](#vertical_interpolation)
- [Time slice](#time_slice)
- [Processing tools](#processing_tools)

Data_information
--------------------------------------------------------------------------------
The online data diagnostics workflow is designed to be as "light" as possible to minimize uncessary enginering complexity. Only necessary information is generated. Fortunately, powerful tools exist for data post-processing. The data are fully compatible with popular tools including NCO, CDO, and NCL. Processing data output from GRIST is very simple and highly efficient, for both weather and climate applications.  

Raw data are organized in 1d (2d sphere), 2d (2d sphere+vertical) and 3d (2d sphere+vertical+ntracer) arrays. Most variables are defined at the node point/primal cell, except:    
(i) normal velocity at edge point; unnecessary for most situations;  
(ii) vorticity at dual-cell/triangle point. 

Horizontal_interpolation
--------------------------------------------------------------------------------
1. One may display the unstructured data directly using NCL.  
2. One may use NCO and CDO to post-interpolate the unstructured data to a regular latitude-longitude mesh. To do so, see post.  
3. One may directly fast-check the raw unstructured data using panoply on a map ('georeferenced Feature-Type'); for this, please first use 'nccopy' to convert the data format from CDF5 to e.g., NetCDF classic.  

Vertical_interpolation
--------------------------------------------------------------------------------
GRIST uses a Dry-mass vertical coordinate, this implies that only the dry hydrostatic surface pressure satisfies the relation set by the vertical coordinate coefficients. Thus, typical interpolation functions for hybrid coordinates (e.g., vinth2p) in NCL are not appropriate for GRIST. To interpolate the model data to a regular pressure-level (implying moist hydrostatic pressure), recommend to use a simple linear interpolation (such as wrf_interp_1d/3d of NCL) for vertical interpolation. To do so, one needs the moist hydrostatic pressure (mpressure) from the model. see post.

Time_slice
--------------------------------------------------------------------------------
GRIST data do not have a time dimension for accumulation. Two time sampling intervals for history files are supported:  
(i) A user-defined frequency that generates data at every N steps (including average and instantaneous state; h1 history);  
(ii) Montly-mean output (h0 history).  
Sampling intervals beyond h0 and h1 at one runtime (e.g., h2, h3) requires code extension.    

One may use NCO for post-processing the data to desired time resolution: to add an unlimited time dimension and concatenate them together.

Processing_tools
--------------------------------------------------------------------------------
1. https://code.mpimet.mpg.de/projects/cdo/
2. http://nco.sourceforge.net/
3. http://www.ncl.ucar.edu/
4. https://www.giss.nasa.gov/tools/panoply/
