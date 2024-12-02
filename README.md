# Dynamical-S-gate-decoding
Codes and circuits for "Transversal Logical Clifford gates on rotated surface codes with reconfigurable neutral atom arrays"

## Required python packages:

Stim, Pymatching (usually included in stim), Sinter, joblib (required for multi-processing), scipy

## Decoding experiments
The following decoding experiments are under circuit noise $p=0.001$. (The circuit noise parameter can be easily changed
in each of the following files.) 

Run 'exp_X_memory.py' -> decode the X-memory circuits

Run 'exp_S_2_distance.py' -> use the plain decoder to decode S-2 circuits with $n_{pad}=\left\lfloor \frac{d}{2} \right \rfloor+1$ and 
$n_{m}=d+1$. The parameters $n_{pad}$ and $n_{m}$ can be adjusted in 'exp.py'.

Run 'exp_S_2_near_time_bd.py' -> use the VTB (or other VTB-based) decoder to decode S-2 circuits with $n_{pad}=2$ and $n_{m}=d$. 
Different decoder choices are also provided in this file. 

Running the above experiments would generate data files containing shot and error information. 


## Reading from data files (in the data_files folder) and plotting
In the data_files folder:

Run 'readnplot_threshold.py' to obtain Fig. 4c.

Run 'readnplot_distance.py' to obtain Fig. 4d.

Run 'readnplot_near_timebd.py' to obtain Fig. 4e. 

