#author: Samuel Ahuno
#purpose: run pycoQC on sequencing_summary.txt files from `guppy`. sequencing metrics was generated by megaloodon pipeline  

import sys
import os
import glob
import subprocess
import pycoQC
#!pip install pycoQC

#how to run 
## $ conda activate pycoQC
## $ python pycoQC_sequencing.py /juno/work/greenbaum/projects/TRI_EPIGENETIC/megalodon/
#python /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/methyl_longRead_wf/scripts/python/pycoQC_sequencing.py /work/greenbaum/projects/ont_pipeline/projects/SPECTRUM_MTHY/results/methylation/megalodon/
#python /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/methyl_longRead_wf/scripts/python/pycoQC_sequencing.py /work/greenbaum/projects/ont_pipeline/projects/DIVYA_BRCA/results/methylation/megalodon/
# '/juno/work/greenbaum/projects/TRI_EPIGENETIC/megalodon/'.split("/")[-2]


#set paths to meagalodon data
#dir_parent = '/juno/work/greenbaum/projects/TRI_EPIGENETIC/megalodon/'
dir_parent = sys.argv[1] #get path from command line
proj_dir = dir_parent.split("/")[-3] #get project name from path
dirs=glob.glob(f'{dir_parent}*/sequencing_summary.txt')#.split("/")[-1]
#dirs.split("/")#[-1]

sampleNames=[] #define empty list to store sample data
for dir in dirs:
    filename=dir.split("/")[-2]
    sampleNames.append(filename)
print(sampleNames)


#subprocess.run(["pycoQC", "-h"]) 
#subprocess.run(["mkdir", "-p", "{proj_dir}/pycoQC/json", "{proj_dir}/pycoQC/html"]) #create dir to keep results
os.system(f'mkdir -p {proj_dir}/pycoQC/json {proj_dir}/pycoQC/html')


from multiprocessing import Pool
def run_pycoQC(sampleNames): #define function to run pycoQC
    os.system(f'pycoQC -f {dir_parent}{sampleNames}/sequencing_summary.txt -o {proj_dir}/pycoQC/html/{sampleNames}_pycoQC_sequencing_summary.html -j {proj_dir}/pycoQC/json/{sampleNames}_pycoQC_sequencing_summary.json --quiet')

#run pycoQC in parallel
pool = Pool(5)
pool.map(run_pycoQC, sampleNames)

#dirs