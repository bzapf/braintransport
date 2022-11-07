import os
import pathlib
import time
jobpath = "/cluster/home/bazapf/nobackup/sleepCode/scripts/saga/slurmouts/"

while True:
    for job in os.listdir(jobpath):
        jobfile = jobpath + job

        jobid = str(pathlib.Path(job).stem)

        file1 = open(jobfile, 'r')


        try:
            Lines = file1.readlines()
        except UnicodeDecodeError:
            print("UnicodeDecodeError at job", job, "will continue to next job")
            continue

        for line in Lines:
            # print(line)

            if "Out Of Memory" in line:
                os.system("scancel " + jobid)
                print("Cancelled job", jobid)


    print("Sleeping for one hour")
    time.sleep(60 * 60)