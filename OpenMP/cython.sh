#!/bin/bash 
#PBS -l nodes=1:ppn=16,walltime=10:00:00
# Define the working directory 
export MYDIR="path_to_your_directory"
cd $PBS_O_WORKDIR 

export JOBNO="`echo $PBS_JOBID | sed s/.master.cm.cluster//`"
export CONF="$MYDIR/machines.$JOBNO" 
for i in `cat $PBS_NODEFILE`; 
do echo $i >> $CONF 
done
export NUMPROC=`cat $PBS_NODEFILE|wc -l`
python ./run_script.py 16
