# Read APC propeller performance data
This repository includes functions and test scripts to import and process APC propeller performance data file in MATLAB. See the test script file to explore how it works. Examples includes:
- plotting all performance data (propeller coefficients and forces)
- interpolating the available data at specific RPM or airspeed

![example](https://user-images.githubusercontent.com/52099779/232206371-5da17179-9c1c-43dc-ae23-cafc77d81a24.png)

### Notes
Developed and tested with MATLAB R2021b.

APC propellers performance data can be found on the website: [https://www.apcprop.com/technical-information/performance-data/](https://www.apcprop.com/technical-information/performance-data/)

All the interpolations in the test script files are made with the native ``scatteredInterpolant`` MATLAB function. Given the available performance data, where datasets at lower RPM provide results in an airspeed range narrower than data at higher RPM, it may be wise to disable the extrapolation of the ``scatteredInterpolant`` function. For instance, the ``testFunction.m`` uses the ``scatteredInterpolant`` function with default methods and may provide bumpy plots at the highest velocities, while the ``testPerfo1.m`` and the ``testPerfo2.m`` script files are more advanced, providing data normalization before interpolation, and avoiding jumps in the plots.
