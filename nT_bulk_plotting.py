#!/usr/bin/python
# Filename: nT_bulk.py
#
# Author: Jacob Vander Vliet
# Version: 2
# Date: 4/12/2013
# Last Modified: 4/12/13
#
# Description:
#   Creates the phase diagram of the gas for cells that
#   contribute to absorption
#
# Instructions:
#   python nT_bulk.py <normalize?> <plot together?>
#
# If normalize is high, the ALL1 and ALL8 are normalized to the SN counts

import numpy as np
import matplotlib.pyplot as plt
from os import system
#import scipy.interpolate
import math
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#import matplotlib as mpl
import sys
from barymp import *

verbose = 0
normalize = int(sys.argv[1])

if len(sys.argv)!=4:
    bulk = 0
else:
    bulk = int(sys.argv[2])

galID_list = ['D9o2', 'D9q', 'D9m4a']
feedback_list = ['dwSN', 'dwALL_1', 'dwALL_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']


file_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/abscells/'
numbins = 50
i = -1
histos = []
snHistos = []
all1Histos = []
all8Hisots = []
xed = []
yed = []
for galID in galID_list:
    i+=1
    for ion in ion_list:
        print galID, '\t', ion

        # Open the data 
        abs_file = galID+'.'+ion+'.bulk_abscells.dat'
        lognH, logT = np.loadtxt(file_loc+abs_file, skiprows=1, usecols=(7, 8), unpack=True)
        print lognH[0], logT[0]

        # Bin the data
        H, xedges, yedges = np.histogram2d( lognH, logT, bins=numbins)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        
        # Normalize each cell by the total number of cells 
        # if normalize is high
        if normalize==1:
            with np.errstate(divide='ignore'):
                H = H / Hfull
                H[Hfull==0] = 0.
            # Mask the bins where the count is zero
            Hmasked = np.ma.masked_where(H==0,H)
        else:
            # Mask the bins where the count is zero
            # Take the log of the count
            Hmasked = np.ma.masked_where(H==0,H)
            Hmasked = np.log10(Hmasked)


        # Rotate and filp the histogram
        Hmasked = np.rot90(Hmasked)
        Hmasked = np.flipud(Hmasked)        
        xed.append(xedges)
        yed.append(yedges)
        histos.append(Hmasked)
        
        if feedback_list[i]=='dwSN':
            snHistos.append(Hmasked)
        elif feedback_list[i]=='dwALL_1':
            all1Histos.append(Hmasked)
        elif feedback_list[i]=='dwALL_8':
            all8Histos.append(Hmasked)
        else:
            print 'Unknown feedback name: ', feedback_list[i]
            sys.exit()


# Normalize dwALL_1 and dwALL_8 to the dwSN counts if
# the normalize flag is high
if normalize==1:
    for i in range(0,len(all1Histos)):
        allH = all1Histos[i]
        snH = snHistos[i]
        allH = allH[i]/snH[i]
        allH[snH==0] = 0.
        all1Histos[i] = allH
    for i in range(0,len(all1Histos)):
        allH = all1Histos[i]
        snH = snHistos[i]
        allH = allH[i]/snH[i]
        allH[snH==0] = 0.
        all1Histos[i] = allH
    
# Create the subplots
fig,((p11,p12,p13), (p21,p22,p23), (p31,p32,p33), (p41,p42,p43)) = plt.subplots(4,3,figsize=(10.2,10.8))
plot_list = [p11, p21, p31, p41, p12, p22, p32, p42, p13, p23, p33, p43]

labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']

labels = ['(a)', '(d)', '(g)', '(j)', 
          '(b)', '(e)', '(h)', '(k)', 
          '(c)', '(f)', '(i)', '(l)']
for i in range(0,len(plot_list)):
    ax = plot_list[i]
    H = histos[i]
    if i<4:     
        # Sn galaxy
        H = snHistos[i]
    elif i>=4 and i<8:
        # All1 galaxy
        H = all1Histos[i-4]
    else:
        # ALL8 galaxy
        H = all8Histos[i-8]
    xedges = xed[i]
    yedges = yed[i]
    mesh = ax.pcolormesh(xedges, yedges, H)
#        mesh = ax.pcolormesh(xedges,yedges,H, vmin=0.0, vmax=1.0)
    ax.set_xlim([-8, 1])
    ax.set_ylim([2,8])
    ax.text(-1, 7, labels[i])
    if i==0:
        ax.set_title('dwSN', size=12)
        ax.set_ylabel('HI \n $\log$ (T) [K] ', multialignment='center')
    elif i==1:
        ax.set_ylabel('MgII \n $\log$ (T) [K] ', multialignment='center')
    elif i==2:
        ax.set_ylabel('CIV \n $\log$ (T) [K] ', multialignment='center')
    elif i==3:
        ax.set_ylabel('OVI \n $\log$ (T) [K] ', multialignment='center')
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==4:
        ax.set_title('dwALL_1', size=12)
    elif i==7:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==8:
        ax.set_title('dwALL_8', size=12)
    elif i==11:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')

    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.set_label('$\log$ (Counts)')
    
dropx = [p12, p22, p32, p13, p23, p33, p11, p21, p31]
plt.setp([a.get_xticklabels() for a in dropx],visible=False)

dropy = [p12, p22, p32, p13, p23, p33, p42, p43]
plt.setp([a.get_yticklabels() for a in dropy],visible=False)

plt.tight_layout()
if normalize==0:
    s = file_loc+'abscell_phase.master_bulk.eps'
else:
    s = file_loc+'abscell_phase.master_bulk_normed.eps'

plt.savefig(s, bbox_inches='tight')
