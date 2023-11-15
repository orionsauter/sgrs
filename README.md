# Release Simulation
Simulation of the S-GRS test mass release process. The dynamics are handled by `TMrel_HG.m` and `TMsys_HG.m`.

## Parametric Study
To run a series of simulations with varying parameters, use `DesignStudy.m` and `AnalyzeStudy.m`. These are set up for HiPerGator. From the login node, the following example runs a series of y-impulse values:
```
module load matlab
matlab -r "DesignStudy('ImpulseStudy','iy',[50e-6;60e-6;70e-6]);exit;"
sbatch Run_ImpulseStudy.sh
mv ImpulseStudy* ImpulseStudy/
matlab -r "AnalyzeStudy('ImpulseStudy','y imp. [kg m/s]');exit;"
``````
Before trying to run any of this code, be sure to change the Slurm inputs to your own information.
