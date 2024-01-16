import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from astropy.io import fits

plt.rcParams.update({"text.usetex": True})

cmapI = cm.gist_heat
cmapP = cm.PRGn


def fig2(res_path='constantCylinder/example1/data/polaris_detector_nr0001.fits.gz'):

    # setup fig 1
    fig = plt.figure(figsize=(9, 1.2))
    gs = fig.add_gridspec(2, 4,
                          height_ratios=[0.1, 1], hspace=0.03, wspace=0.22)

    axs = [fig.add_subplot(gs[1, i]) for i in range(4)]
    caxI = fig.add_subplot(gs[0, 0])
    caxP = fig.add_subplot(gs[0, 1:])

    # total and polarized intensity maps for detector I
    # shape: (type,Nlambda,Npixel,Npixel)
    imaps = fits.getdata(res_path)

    # account for detector rotation (x->-z, y->y):
    imaps = np.swapaxes(imaps, -1, -2)
    imaps = np.flip(imaps, axis=-2)

    # plot total intensity map
    im = axs[0].imshow(imaps[0][0],
                       cmap=cmapI,
                       norm=colors.LogNorm(vmin=np.amax(imaps[0][0])*1e-7,
                                           vmax=np.amax(imaps[0][0])),
                       extent=[-3, 3, -3, 3],
                       origin="lower"
                       )

    # polarized intensity Ip = sqrt(Q²+U²+V²)
    ip = np.sqrt(np.sum([imaps[i+1][0]**2 for i in range(3)], axis=0))
    for i in range(3):
        imp = axs[i+1].imshow(imaps[i+1][0] / ip,
                              cmap=cmapP,
                              norm=colors.Normalize(vmin=-1, vmax=1),
                              extent=[-3, 3, -3, 3],
                              origin="lower"
                              )

    # axis
    for ax, l in zip(axs, ["$I$", "$Q/I_P$", "$U/I_P$", "$V/I_P$"]):
        ax.set_xlabel("$y$ in cm")
        if l == "$I$":
            textc = "w"
        else:
            textc = "k"
        ax.text(0.05, 0.95, l, c = textc, va="top", ha="left", transform=ax.transAxes)
        ax.set_ylim(-1.5, 1.5)
    axs[0].set_ylabel("$z$ in cm")

    # layout colorbar I
    cbarI = fig.colorbar(im, cax=caxI, orientation="horizontal")
    cbarI.ax.set_xlabel("$I$ in Wm$^{-1}$m$^{-2}$")
    cbarI.ax.xaxis.set_ticks_position("top")
    cbarI.ax.xaxis.set_label_position("top")

    # layout Q, U, V
    cbarP = fig.colorbar(imp, cax=caxP, orientation="horizontal")
    cbarP.ax.set_xlabel("normalized polarized intensities")
    cbarP.ax.xaxis.set_ticks_position("top")
    cbarP.ax.xaxis.set_label_position("top")

    return fig


def fig3(res_path='constantCylinder/example1/data/polaris_detector_nr0002.fits.gz'):

    # setup fig 1
    fig = plt.figure(figsize=(9, 1.2))
    gs = fig.add_gridspec(2, 4,
                          height_ratios=[0.1, 1], hspace=0.03, wspace=0.22)

    axs = [fig.add_subplot(gs[1, i]) for i in range(4)]
    caxI = fig.add_subplot(gs[0, 0])
    caxP = fig.add_subplot(gs[0, 1:])

    # total and polarized intensity maps for detector I
    # shape: (type,Nlambda,Npixel,Npixel)
    imaps = fits.getdata(res_path)

    # account for detector rotation (x->x, y->-z):
    imaps = np.flip(imaps, axis=-2)

    # plot total intensity map
    im = axs[0].imshow(imaps[0][0],
                       cmap=cmapI,
                       norm=colors.LogNorm(vmin=np.amax(imaps[0][0])*1e-3,
                                           vmax=np.amax(imaps[0][0])),
                       extent=[-3, 3, -3, 3],
                       origin="lower"
                       )

    # polarized intensity Ip = sqrt(Q²+U²+V²)
    ip = np.sqrt(np.sum([imaps[i+1][0]**2 for i in range(3)], axis=0))
    for i in range(3):
        imp = axs[i+1].imshow(imaps[i+1][0] / ip,
                              cmap=cmapP,
                              norm=colors.Normalize(vmin=-1, vmax=1),
                              extent=[-3, 3, -3, 3],
                              origin="lower"
                              )

    # axis
    for ax, l in zip(axs, ["$I$", "$Q/I_P$", "$U/I_P$", "$V/I_P$"]):
        ax.set_xlabel("$x$ in cm")
        if l == "$I$":
            textc = "w"
        else:
            textc = "k"
        ax.text(0.05, 0.95, l, c = textc, va="top", ha="left", transform=ax.transAxes)
        ax.set_ylim(-1.5, 1.5)
    axs[0].set_ylabel("$z$ in cm")

    # layout colorbar I
    cbarI = fig.colorbar(im, cax=caxI, orientation="horizontal")
    cbarI.ax.set_xlabel("$I$ in Wm$^{-1}$m$^{-2}$")
    cbarI.ax.xaxis.set_ticks_position("top")
    cbarI.ax.xaxis.set_label_position("top")

    # layout Q, U, V
    cbarP = fig.colorbar(imp, cax=caxP, orientation="horizontal")
    cbarP.ax.set_xlabel("normalized polarized intensities")
    cbarP.ax.xaxis.set_ticks_position("top")
    cbarP.ax.xaxis.set_label_position("top")

    return fig


def fig4(res_path='constantCylinder/example2/data/polaris_detector_nr0002.fits.gz'):

    # setup fig 1
    fig = plt.figure(figsize=(9, 1.2))
    gs = fig.add_gridspec(2, 4,
                          height_ratios=[0.1, 1], hspace=0.03, wspace=0.22)

    axs = [fig.add_subplot(gs[1, i]) for i in range(4)]
    caxI = fig.add_subplot(gs[0, 0])
    caxP = fig.add_subplot(gs[0, 1:])

    # total and polarized intensity maps for detector I
    # shape: (type,Nlambda,Npixel,Npixel)
    imaps = fits.getdata(res_path)

    # account for detector rotation (x->x, y->-z):
    imaps = np.flip(imaps, axis=-2)

    # plot total intensity map
    im = axs[0].imshow(imaps[0][0],
                       cmap=cmapI,
                       norm=colors.LogNorm(vmin=np.amax(imaps[0][0])*1e-3,
                                           vmax=np.amax(imaps[0][0])),
                       extent=[-3, 3, -3, 3],
                       origin="lower"
                       )

    # polarized intensity Ip = sqrt(Q²+U²+V²)
    ip = np.sqrt(np.sum([imaps[i+1][0]**2 for i in range(3)], axis=0))
    for i in range(3):
        imp = axs[i+1].imshow(imaps[i+1][0] / ip,
                              cmap=cmapP,
                              norm=colors.Normalize(vmin=-1, vmax=1),
                              extent=[-3, 3, -3, 3],
                              origin="lower"
                              )

    # axis
    for ax, l in zip(axs, ["$I$", "$Q/I_P$", "$U/I_P$", "$V/I_P$"]):
        ax.set_xlabel("$x$ in cm")
        if l == "$I$":
            textc = "w"
        else:
            textc = "k"
        ax.text(0.05, 0.95, l, c = textc, va="top", ha="left", transform=ax.transAxes)
        ax.set_ylim(-1.5, 1.5)
    axs[0].set_ylabel("$z$ in cm")

    # layout colorbar I
    cbarI = fig.colorbar(im, cax=caxI, orientation="horizontal")
    cbarI.ax.set_xlabel("$I$ in Wm$^{-1}$m$^{-2}$")
    cbarI.ax.xaxis.set_ticks_position("top")
    cbarI.ax.xaxis.set_label_position("top")

    # layout Q, U, V
    cbarP = fig.colorbar(imp, cax=caxP, orientation="horizontal")
    cbarP.ax.set_xlabel("normalized polarized intensities")
    cbarP.ax.xaxis.set_ticks_position("top")
    cbarP.ax.xaxis.set_label_position("top")

    return fig


if __name__ == '__main__':

    fig2 = fig2()
    fig3 = fig3()
    fig4 = fig4()

    plt.show()