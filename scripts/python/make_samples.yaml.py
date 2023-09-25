#name:Samuel ahuno
#date: Sept 23rd 2023
#create yaml files

import glob as glob

with open("/juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/config/samples_automatic.yaml", "w") as f:
    f.write("samples:\n")
    for file_path in glob.glob("/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/*/*_5mC.bed", recursive=True):
        sample_name = file_path.split("/")[-1].split(".")[0]
        f.write(f"    {sample_name}: {file_path}\n")