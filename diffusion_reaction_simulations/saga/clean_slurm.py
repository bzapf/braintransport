import os

p = "/cluster/home/bazapf/nobackup/sleepCode/scripts/saga/slurmouts/"
for file in os.listdir(p):
    filepath = p + file
    if filepath.endswith(".out") and "slurm" in filepath:
        os.system("rm " + filepath)
print("Cleaning userwork, this might take a while...")
for file in os.listdir("/cluster/work/users/bazapf"):
    if len(file) == len("4791857"):
        os.system("rm -r /cluster/work/users/bazapf/" + file)
