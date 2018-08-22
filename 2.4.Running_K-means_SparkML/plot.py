import numpy as np
import sys
import os.path
from subprocess import call
import matplotlib.colors as cls
import matplotlib.pyplot as plt

# if len(sys.argv) < 3:
#     sys.exit("Need 2 inputs: km and sid")
# km = int(sys.argv[1])  # Number of Clusters
km = 10
# sid = int(sys.argv[2])  # id
sid = 1
# dir1 = './CTD/'
# mdnm = 'Aqua_b42_TR'
# ctdfnm = 'MODIS_{}.cent_km{:02d}_sid{:02d}'.format(mdnm,km,sid)
# fnm=dir1+ctdfnm+'.dpdat'
ncl = km
nelem = 42

obsctd = centers
obsctd = obsctd.reshape([ncl, nelem])

obscf = np.sum(obsctd, axis=1) * 100.
np.set_printoptions(precision=3, suppress=True)
print(obscf)


def cent_show(ax1, i, ctd):

    global ncl
    nx = 6  # TAU (Optical Thickness)
    ny = 7  # CTP
    if len(ctd.reshape(-1)) != nx * ny:
        print("Error: centroid data size is bad:", ctd.shape)
        sys.exit()

    cm = plt.cm.get_cmap('jet', 512)
    cmnew = cm(np.arange(512))
    cmnew = cmnew[72:, :]
    newcm = cls.LinearSegmentedColormap.from_list("newJET", cmnew)
    newcm.set_under('white')

    props = dict(norm=cls.LogNorm(vmin=0.1, vmax=30), cmap=newcm, alpha=0.9)
    pic1 = ax1.imshow(ctd, interpolation='nearest', aspect=0.84, **props)

    # Axis Control
    xlabs = ['0', '1.3', '3.6', '9.4', '23', '60', '379']
    ylabs = [1000, 800, 680, 560, 440, 310, 180, 10]

    ax1.set_xlim(-0.5, 5.5)
    ax1.set_ylim(-0.5, 6.5)
    ax1.set_xticks(np.ones(nx + 1) * range(nx + 1) - 0.5)
    if i < ncl - 2:
        ax1.set_xticklabels([])
    else:
        ax1.set_xticklabels(xlabs)

    ax1.set_yticks(np.ones(ny + 1) * range(ny + 1) - 0.5)
    if i % 3 == 1:
        ax1.set_yticklabels(ylabs)
    elif i % 3 == 2:
        ax1.set_yticklabels([])
    else:  # i%2 == 0:
        ax1.set_yticklabels(ylabs)
        ax1.yaxis.tick_right()

    # Add colorbar
    if i == 3:
        tt = [0.1, 0.3, 1, 3, 10, 30]
        # tt=[0.2,0.5,1,2,5,10,30]
        tt2 = [str(x) + '%' for x in tt]

    for j in range(7):
        for i in range(6):
            if abs(ctd[j, i]) > 5.0:
                # ax1.annotate(str(ctd[j,i]),xy=(ix,iy))
                ax1.annotate("%4.1f" % (ctd[j, i]), xy=(i, j), ha='center', va='center', stretch='semi-condensed',
                             fontsize=10)
    return pic1


def add_colorbar_horizontal(ax1, pic1, tt, tt2=None):
    # Add colorbar
    if tt2 == None:
        tt2 = tt
    pos1 = ax1.get_position().bounds  ##<= (left,bottom,width,height)
    cb_ax = fig.add_axes([0.1, pos1[1] - 0.07, 0.8, 0.015])
    # cb_ax = fig.add_axes([1.0,0.2,0.03,0.6])  ##<= (left,bottom,width,height)
    cb = fig.colorbar(pic1, cax=cb_ax, orientation='horizontal', ticks=tt, extend='both')
    cb.ax.set_xticklabels(tt2, size=12, stretch='condensed')
    return cb


def cent_show_common(ax1, i, cf):

    # add a title.
    subtit = "CR" + str(i) + ". CF=%4.1f%%" % (cf)
    print(subtit)
    ax1.set_title(subtit, x=0.0, ha='left', fontsize=12, stretch='semi-condensed')

    # Draw Guide Line
    ax1.axvline(x=1.5, linewidth=0.7, color='k', linestyle=':')
    ax1.axvline(x=3.5, linewidth=0.7, color='k', linestyle=':')
    ax1.axhline(y=1.5, linewidth=0.7, color='k', linestyle=':')
    ax1.axhline(y=3.5, linewidth=0.7, color='k', linestyle=':')

    # Ticks
    ax1.tick_params(axis='both', which='major', labelsize=10)

    return

# plotting basics
fig, axs = plt.subplots(4, 3)  ## (ny,nx)
fig.set_size_inches(7.5, 9.6)  ## (lx,ly)
plt.suptitle("MODIS Aqua Tropical CRs, K={}, id={}".format(km, sid), fontsize=18, y=0.98)
lf = 0.09;
rf = 0.91
bf = 0.06;
tf = 0.92
fig.subplots_adjust(hspace=0.18, wspace=0.06, left=lf, right=rf, top=tf, bottom=bf)

for ii, ax in enumerate(axs.flat):
    if ii < km:
        vv = np.copy(obsctd[ii, :].reshape([7, 6])) * 100.
        vv = vv[::-1, :]
        pic1 = cent_show(ax, ii + 1, vv)
        cent_show_common(ax, ii + 1, obscf[ii])
        if ii == int((ncl - 1) / 3) * 3:
            ax.set_xlabel('Optical Thickness', fontsize=13)
            ax.set_ylabel('Pressure (hPa)', fontsize=13, labelpad=0)
        if ii == km - 1:
            # Add Color Bar
            tt = [0.1, 0.3, 1, 3, 10, 30]
            tt2 = [str(x) + '%' for x in tt]
            cb = add_colorbar_horizontal(ax, pic1, tt, tt2)
    else:
        ax.set_visible(False)

# -----------------------------------

# plt.tight_layout()

outdir = "./"
fnout = ctdfnm + ".png"

# Show or Save
# plt.show()
plt.savefig(outdir + fnout, bbox_inches='tight', dpi=175)
