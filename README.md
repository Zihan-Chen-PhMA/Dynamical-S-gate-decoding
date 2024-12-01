# Dynamical-S-gate-decoding
Codes and circuits for "Transversal Logical Clifford gates on rotated surface codes with reconfigurable neutral atom arrays"

## Required python packages:
Running our codes on python=3.10 with the newest versions of the following packages is OK. 

Stim, Pymatching (usually included in stim), Sinter, joblib (required for multi-processing), scipy

## Decoding experiments
The following decoding experiments are under circuit noise $p=0.001$. (The circuit noise parameter can be easily changed
in each of the following files.) 

Run 'exp.py' -> decode the X-memory circuits

Run 'exp.py' -> use the plain decoder to decode S-2 circuits with $n_{pad}=\left\lfloor \frac{d}{2} \right \rfloor+1$ and 
$n_{m}=d+1$. The parameters $n_{pad}$ and $n_{m}$ can be adjusted in 'exp.py'.

Run 'exp.py' -> use the VTB-PR (or plain/VTB/VTB-FR/FR) decoder to decode S-2 circuits with $n_{pad}=2$ and $n_{m}=d$. 
Different decoder choices are provided in ''. 

Running the above experiments would generate data files containing shot and error information. 


## Reading from data files (in the data folder) and plotting
In the data folder:

Run 'readnplot_threshold.py' to obtain Fig. 4c.

Run 'readnplot_distance.py' to obtain Fig.4d.

Run 'readnplot_near_time_bd.py' to obtain Fig.4e. 

