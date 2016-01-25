# motionEnergy

This code can be run on LISA, using the stopos system to submit jobs efficiently. 

First, compile matlab code:
module load matlab
module load mcr
mcc -m -R -singleCompThread MotionEnergyWrap2.m (to disable multithreading)

(see MotionEnergyWrap2.m for the use of scratch space)

Then, submit jobs using stopos
module load stopos
stopos create -p pool
stopos add paramfile -p pool
qsub -t 1-50 stoposjob (submit this njobs/8 times, every node has minimum 8 cores)

check progress with: stopos status