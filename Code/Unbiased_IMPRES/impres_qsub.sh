#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=2G
#$ -N Auslander
#$ -t 1-500:1
#$ -tc 200
#$ -e impres_e.txt
#$ -o impres_o.txt

#Runs Auslander.m X times as an array job
# -t determines the number of independent feature sets to be generated

export PATH=$(echo $PATH | sed -e 's;R2013a;R2018a;')

iteration=${SGE_TASK_ID}
directory="/mnt/grid/atwal/hpc/data/data/jacarter/Auslander/Run"

cd $directory

matlab -nodisplay -nosplash -nodesktop -r "run('Auslander.m');exit;"
