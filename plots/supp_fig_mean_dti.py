import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from pathlib import Path
from scipy import stats
from matplotlib.markers import TICKDOWN
import matplotlib

"""

Script to generate bar plot that compares mean diffusivity in white / gray between groups

"""

try:
    from dolfin import *
except:
    pass

# matplotlib.rcParams["lines.linewidth"] = 3
# matplotlib.rcParams["axes.linewidth"] = 3
# matplotlib.rcParams["axes.labelsize"] = "xx-large"
# matplotlib.rcParams["grid.linewidth"] = 1
# matplotlib.rcParams["xtick.labelsize"] = "xx-large"
# matplotlib.rcParams["ytick.labelsize"] = "xx-large"
# matplotlib.rcParams["legend.fontsize"] = "xx-large"
# matplotlib.rcParams["font.size"] = 14

params = {'font.size': 9,}
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams.update(params) 

from matplotlib.markers import TICKDOWN
def significance_bar(start, end, height, fontsize, displaystring, text_dh, linewidth=1., 
                    markersize=5, boxpad=0.2, color='k'):
    assert start != end
    # draw a line with downticks at the ends
    plt.plot([start, end], [height, height], '-', color=color, lw=linewidth,
             marker=TICKDOWN, markeredgewidth=linewidth, markersize=markersize)
    # draw the text with a bounding box covering up the line
    plt.text(0.5 * (start + end), text_dh + height, displaystring, ha='center', va='center',
             # bbox=dict(facecolor='1.', edgecolor='none', # boxstyle='Square,pad=' + str(boxpad)), 
             size=fontsize)


path_to_files = "/home/basti/Dropbox (UiO)/Sleep/"

from definitions import groups

group = "all"
pats = groups[group]

# fig1 = plt.figure()
# ax1 = fig1.gca()

# fig2 = plt.figure()
# ax2 = fig2.gca()

# fig3 = plt.figure()
# ax3 = fig3.gca()

data_dict = {}

for pat in pats:

    break

    dt = 1800
    res = 32
    if pat == 241:
        res = 64

    mesh = Mesh()

    global dest

    hdf = HDF5File(mesh.mpi_comm(), "%s/%.3d/mesh/parenchyma%d.h5" %
            (path_to_files, pat, res), "r")

    hdf.read(mesh, "/mesh", False)
    SD = MeshFunction("size_t", mesh, mesh.topology().dim())
    hdf.read(SD, "/subdomains")
    
    bnd = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    hdf.read(bnd, "/boundaries")
    lookup_table = MeshFunction("size_t", mesh, mesh.topology().dim())
    hdf.read(lookup_table, '/lookup_table')
    lookup_boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    hdf.read(lookup_boundary, "/lookup_boundary")

    TensorSpace = TensorFunctionSpace(mesh, 'DG', 0)
    MDSpace = FunctionSpace(mesh, 'DG', 0)
    MD = Function(MDSpace)
    Kt = Function(TensorSpace)
    hdf.read(MD, '/MD')
    hdf.read(Kt, '/DTI')

    V = FunctionSpace(mesh, "CG", 1)

    write_regions = np.unique(lookup_table.array()[:])

    # GRAY = 1. WHITE = 2. BRAIN STEM = 3.
    dx_SD = Measure('dx')(domain=mesh, subdomain_data=SD)
    dx_lookup = Measure('dx')(domain=mesh, subdomain_data=lookup_table)
    ds_lookup = Measure('ds')(domain=mesh, subdomain_data=lookup_boundary)

    mean_d_gray = assemble(MD * dx_SD(1)) / assemble(1 * dx_SD(1))

    mean_d_white = assemble(tr(Kt) * dx_SD(2)) / (3 * assemble(1 * dx_SD(2)))

    breakpoint()

    # breakpoint()
    # form = 3 - tr((Kt / tr(Kt)) ** 2) ** (-1)
    # fa2 = assemble(sqrt(0.5 * form) * dx_SD(1))

    data_dict[pat] = [mean_d_gray, mean_d_white]

    # breakpoint()

    if pat in sleepers:
        data_dict[pat] = data_dict[pat] + ["sleep"]
    else:
        data_dict[pat] = data_dict[pat] + ["no sleep"]

    # ax2.plot(ts, avgs_ds_at_ts, linewidth=0, color=c, marker="o")
    
    # ax1.plot(t, avgds, linewidth=1, marker="x", color=c)

    # print(len(t))

def dtibarplot(fs, dpi, size_inches):

    statistics_path = "/home/basti/Dropbox (UiO)/Sleep/PLOTS/concentration_statistics/"

    patdf = pd.read_csv(statistics_path + "dti_stats_all.csv", engine='python')
    print(patdf)

    meand_white_nosleep = patdf.iloc[:6, 2] * 1e3
    meand_gray_nosleep = patdf.iloc[:6, 1] * 1e3

    meand_white_sleep = patdf.iloc[6:, 2] * 1e3
    meand_gray_sleep = patdf.iloc[6:, 1] * 1e3

    sleep_means = np.array([np.mean(meand_gray_sleep), np.mean(meand_white_sleep)])
    nosleep_means = np.array([np.mean(meand_gray_nosleep), np.mean(meand_white_nosleep)])

    sleep_std = stats.sem(np.array([np.array(meand_gray_sleep).tolist(), np.array(meand_white_sleep).tolist()]), axis=1, nan_policy="omit")
    nosleep_std = stats.sem(np.array([np.array(meand_gray_nosleep).tolist(), np.array(meand_white_nosleep).tolist()]), axis=1, nan_policy="omit")

    # breakpoint()

    labels = ["gray", "white"]

    # breakpoint()
    x = np.arange(len(labels))  # the label locations
    width = 0.3 # the width of the bars
    capsize_ = 2
    fontsize = fs

    
    cs = "tab:blue" #"royalblue"
    cns = "tab:red"

    dtifig, ax = plt.subplots(dpi=dpi)
    dtifig.set_size_inches(size_inches[0], size_inches[1]) 

    rects1 = ax.bar(x - width / 2, sleep_means, width, color=cs, label='Sleep', 
            yerr=sleep_std, capsize=capsize_,
            )
    rects2 = ax.bar(x + width / 2, nosleep_means, width, color=cns, label='Sleep deprivation',
            yerr=nosleep_std, capsize=capsize_,
            )
    for i in range(2):
        bar_centers = x[i] + np.array([-width / 2., width / 2])
        h = 0.03 + max(nosleep_means[i] + nosleep_std[i] / 2, sleep_means[i] + sleep_std[i] / 2)
        text_dh = 0.01
        # significance_bar(bar_centers[0], bar_centers[1], fontsize=fontsize, height=h, displaystring="n.s.", text_dh=text_dh)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('mean diffusivity $(10^{-3}\,\mathrm{mm}^2\, \mathrm{s}^{-1})$', fontsize=fontsize)
    # ax.set_title('Scores by group and gender')
    ax.set_xticks(x, labels, )
    # ax.tick_params(axis='y', which='major', labelsize=fontsize)
    ax.legend(fontsize=fontsize)

    ax.set_ylim(0.8, 1.02)
    plt.locator_params(axis='y', nbins=6)

    # ax.bar_label(rects1, padding=3)
    # ax.bar_label(rects2, padding=3)
    # dtifig.tight_layout()

    

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    # dtifig.savefig(plotpath + "avg_mean_d.pdf", dpi=3*dpi)
    dtifig.savefig(plotpath + "avg-mean-d.png", dpi=3*dpi)

    # fig.savefig("/home/basti/Dropbox (UiO)/Sleep/PLOTS/concentration_statistics/avg_mean_d.pdf")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    plotpath = "/home/basti/Dropbox (UiO)/Sleep/PLOTS/concentration_statistics/"
    plotpath = "/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/Brain influx and clearance during sleep and sleep deprivation/figures/data/"

    latex_textwidth = 7.02352416667 # inches
    figwidth = latex_textwidth / 3
    figheight = figwidth * 1

    size_inches = (figwidth, figheight)

    # from scripts.plots.generic_plot import figsize, dpi
    dtibarplot(fs=None, dpi=300, size_inches=size_inches)


    MD_water = 3e-3
    MD_gadovist = 3.8e-4 # (taken from Gd-DPTA)
    scale_diffusion_gad = MD_gadovist / MD_water

    MD_gray = 0.00092582074369026 # mean from all pats with DTI in gray
    MD_white = 0.000867643624860488 # mean from all pats with DTI in white
    MD_brainstem = 0.0010602528519216016 # mean from all pats with DTI in brainstem
    
    print("Mean gadobutrol diffusivities")
    for d in [MD_gray, MD_white, MD_brainstem]:
        print(d * scale_diffusion_gad * 1e4)
    
    exit()
