import os
import json
import pathlib

intervals = [
    # (0, 0.001),
    (1.2, 2.6), 
    (4.5, 7.4), (20, 30), (40, 56)
]

labels = [# r"0", 
            r"$\sim$2", r"$\sim$6", r"$\sim$24", r"$\sim$48"]

ylabels = {"avg": "Brain-wide (mmol)", "white": "Subcortical white matter (mmol)", "gray": "Cerebral cortex (mmol)"}

groups = {}

sleep = ["105", "175", "176", "178", "183", "190", "191", "205", "215", "218", "228", "240"]
sleepdep = ["199", "227", "230", "235", "236", "241"]

sleep += ['002', '078', '091', "127", "172"]
sleepdep += ["249"]

groups["sleep"] = sleep
groups["sleepdep"] = sleepdep

nodti = ["002", "078", "091", "127", "172", "249"]

groups["nodti"] = nodti

groups["t1map"] = ["105", "175", "178", "183", "190", "191", "205", "215", 
                    "228", "240", "199", "227", "235", "236", "241"]

groups["all"] =sleep + sleepdep

groups["nodtinot1"] = ['002', '078', '091', '127', '172', '176', '218', '230', '249']

groups["dti"] = [x for x in groups["all"] if x not in groups["nodti"]]

for key, group in groups.items():
    groups[key] = sorted(group)

pats = groups["all"]



pats = ['002','078','091','105','127','172','175','176','178','183','190','191','199','205','215','218','227','228','230','235','236','240','241','249']

# print(pats)




nonconverged = ["091", "105", "172", "175", "176", "178", "190", "191", "235", "249"]



# def resultfoldername(pat):
#     return "alphatests"



datafolder = "/home/basti/Dropbox (UiO)/Sleep/"



def make_best_subfolder(pat, resultfolder, subfoldername):
    
    pat = str(pat)
    subfolder = os.path.join(pat, resultfolder, subfoldername, "")

    if (not os.path.isdir(subfolder)):
        # print(pat, "has no iterfolder, subfolder=", subfolder)
        return None

    subfolders = sorted(os.listdir(subfolder), key=lambda x: int(x[1:]), reverse=True)

    if len(subfolders) == 0:
        # print("Empty iterfolder for", pat, "continue")
        return None
    else:
        # print(pat, subfolders[0])
    
        subfolder = os.path.join(subfolder, subfolders[0], "")

    if (not os.path.isdir(subfolder)) or (not os.path.isfile(subfolder + "/params")) or (len(os.listdir(subfolder)) == 0):
        print(pat, "is missing, CONTINUE")

        return None

    return os.path.join(datafolder, subfolder)



def reaction_resultfolder(pat, best=True, k=None, subfoldername="avgDTIavgT1"):

    return_folder = None

    if best:
        
        return make_best_subfolder(pat, resultfolder="diffusion_reaction", subfoldername=subfoldername)

    else:
        assert k is not None
        assert "k" in k

        return_folder = os.path.join(datafolder, pat, "diffusion_reaction", "avgDTIavgT1", k, "")
        
        if not os.path.isdir(return_folder):
            # print(k)
            return None

    try:
        params = json.load(open(return_folder + "params"))
    except json.decoder.JSONDecodeError:
        print("vim " + return_folder + "params")
    
    if not params["concfolder"] == "FIGURES_CONC":
        assert pat in groups["nodtinot1"]

    return return_folder






