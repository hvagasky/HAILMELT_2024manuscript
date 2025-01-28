Code Associated with Vagasky et al. (2025, MWR)

HAILMELT is a Python-based one-dimensional hall fall trajectory model that is based on HAILCAST (Poolman 1992, Brimelow et al. 2002, 2006; Jewell and Brimelow 2009, Adams-Selin and Zeigler 2016, Adams-Selin et al. 2019, Adams-Selin 2023, and Pounds et al. 2024). HAILMELT is unique since it is a simplified hail trajectory model that only simulates falling melting hail and includes several melting and terminal velocity parameterizations. 

Model details:
- Hail in HAILMELT is not allowed to grow via collection or deposition and the model stops if the hail is swept upward by an updraft. 
- Hail is assumed spherical and sheds water when the water mass exceeds 2.0E-4 kg
- HAILMELT's time resolution is currently set to 1 second
- Included parameterizations: 
    - Melting: Heymsfield 1987, Wang and Chueh 2020, and Zheng and List 1994
    - Terminal Velocity: Rasmussen and Heymsfield 1987, diameter – velocity relationship from Heymsfield et al. 2018

Model Input:
- Sounding text file (usually obtained from https://weather.uwyo.edu/upperair/sounding.html)
- Assumed initial hail size distribution (Can be read in from an input file or manually defined within the code)
- Assumed vertical velocity profile (Currently this is assumed to be 0 m/s below the LCL and 5 m/s above the LCL).

Versions of the model 
- HAILMELT.ipynb - Basic verison
- HAILMELT_EnvironmentalSensitivity.ipynb - Runs a series of basic simulations, each with slightly perturbed environment conditions, vertical velocity is assumed to be 0 m/s in each of these simulations. 
- HAILMELT_VertVelSensitivity.ipynb - Runs a series of basic simulations, each with different vertical velocities

How to run the basic model: 
- Prior to running HAILMELT the f2py Python module must be installed (https://numpy.org/doc/stable/f2py/)
    - f2py used to compile each of the fortran scripts so that they can be used within Python.
    - Once f2py is installed and is activated within a python environment Fortran code can be compiled from the command line using the following syntax: python -m numpy.f2py -c -m FortranFileName FortranFileName.f90
- HAILMELT.ipynb contains the basic HAILMELT code and is used to run the model.
    - All inputs that must be made are within the first 11 sections of the code. Variables to update include:
          - Directory and file names of the input data is located
          - Directory to store figures
          - Set flag to determine if shedding is included and plots showing the time evolution of shedding
          - Set flad to determine if the environmental and hail soundings are plotted
          - Set the model time step and maximum duration of the simulations
          - Set the height profile
          - Portions of the code related to reading in the environmental soundings may need to be slightly modified in order to properly read in exact text file.
    - Once all the inputs are set the code can be run.

Basic model output:
- Plots:
    - Time series of hail height and hail diameter
    - Hail diameter and height distribution
    - Parameterization comparison matrix plots of hail diameter, height, and duration
- Netcdf File
    - Temperature and dew point temperature profiles of the hail environment
    - Time series of hail diameter, height, and terminal velocity.
        - The arrays have dimensions of (# melt parameterizations, # terminal parameterizations, # hail stones). The indices in the melting parameterization dimension correspond to: 0=RH87, 1=WC2020, 2=ZL1994. The indices in the terminal velocity parameterizations correspond to: 0=H2018, 1=RH87
          
The versions of the codes that run the environmental and vertical velocity senstivity are run similarly and have similar outputs. All necessary changes are made in the first 11 sections of code. 

Questions or comments can be listed in the discussion section of emailed to hvagasky@aer.com




