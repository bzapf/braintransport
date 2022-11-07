import numpy as np
import os
import matplotlib.pyplot as plt
import scipy


def correlationplot(x_data, y_data, ylabel, xlabel, roi, colors, fontsize=18, labelsize=18, legendsize=18, dpi=300,
                    figsize=None,
                    aspect=1, title="", labels=None, figurepath=None, xlim=None, ylim=None):

    plt.close("all")

    if labels is None:
        labels = ["" for _ in range(y_data.shape[-1])]

    ms = 12

    assert len(y_data.shape) == 3
    assert len(x_data.shape) == 3

    assert x_data.shape[-1] == 1

    assert np.sum((np.isnan(x_data))) <= 1
    assert np.sum((np.isnan(y_data))) <= y_data.shape[-1] # only one (for 241)

    x_data = x_data[..., 0]

    for tid, tlabel in zip(range(y_data.shape[0]), [r"$\sim 2\,$h", r"$\sim 6\,$h", r"$\sim 24\,$h", r"$\sim 48\,$h"]):


        xd_t = x_data[tid, :]
        xd_t = xd_t[~np.isnan(xd_t)]


        tests = []

        cs = []

        for j in range(y_data.shape[-1]):

            conc = y_data[tid, :, j]
            conc = conc[~np.isnan(conc)]

            # breakpoint() 

            test = scipy.stats.pearsonr(xd_t, conc)
            
            print("tid", tid, "pearson r:", format(test[0], ".2f"))

            tests.append(test)
            cs.append(conc)

        fig1, ax1 = plt.subplots(dpi=dpi, figsize=figsize)
        
        sortarr = np.argsort(xd_t)

        loc = "center"

        # if len(title) > 0:
        #     loc = "left"

        ax1.set_title(title + tlabel, fontsize=fontsize, loc=loc)

        count = 0

        coloriter = iter(colors)
        labeliter = iter(labels)

        markers = iter(["o", "*", "s"])


        for conc, test in zip(cs, tests):

            label2 = next(labeliter) + "r=" + format(test[0], ".2f")

            # breakpoint()
            ax1.plot(xd_t[sortarr], conc[sortarr], c=next(coloriter), label=label2,
                    markersize=ms,
                    linewidth=0, marker=next(markers))

        if y_data.shape[-1] == 3:
            if "ds" in roi:
                ax1.set_ylim(0, 0.1)
                ax1.set_xlim(0, 0.1)
            else:
                ax1.set_ylim(0, 0.25)
                ax1.set_xlim(0, 0.25)

        if xlim is not None:
            ax1.set_xlim(xlim)
            assert ylim is not None
            ax1.set_ylim(ylim)
            
            xticks = np.linspace(xlim[0], xlim[1], 4)
            yticks = np.linspace(ylim[0], ylim[1], 4)
            plt.xticks(xticks, [format(x, ".1f") for x in xticks])
            plt.yticks(yticks, [format(x, ".1f") for x in yticks])

        ax1.set_box_aspect(aspect=aspect)

        plt.ylabel(ylabel, fontsize=labelsize)
        plt.xlabel(xlabel, fontsize=labelsize)

        ax1.tick_params(axis='x', which='major',  labelsize=labelsize)        
        ax1.tick_params(axis='y', which='major', labelsize=labelsize)

        # plt.locator_params(axis='y', nbins=5)
        # plt.locator_params(axis='x', nbins=5)



        plt.legend(fontsize=legendsize)
        plt.tight_layout()

        if figurepath is not None:
            assert ".png" not in figurepath
            plt.savefig(figurepath + str(tid) + ".png", dpi=600)

        else:
            plt.show()
            exit()

