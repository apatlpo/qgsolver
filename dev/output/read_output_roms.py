#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Read result from a PV inversion
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, MFDataset
import re

# maybe temporary
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

d2r = np.pi/180.


class qgrun():
    ''' Container for now of useful info
    '''
    def __init__(self, path):
        
        #
        self.path = path
        print 'Load run from following path: '+path
        
        # load log and get useful info
        self.read_log()
        
        # load grid
        nc = Dataset(path+'input/roms_metrics.nc','r')
        self.x = nc.variables['lon'][:,:]
        self.y = nc.variables['lat'][:,:]
        #self.x = nc.variables['lon'][self.jstart:self.jend,self.istart:self.iend]
        #self.y = nc.variables['lat'][self.jstart:self.jend,self.istart:self.iend]
        nc.close()
        
        # load data
        self.nc = MFDataset(path+'output/output*.nc')
        
    def read_log(self):
        # Read output.mpi to extract parameters
        logfile = os.path.join(self.path, 'qgsolver/output.mpi')
        f = open(logfile)
        for line in iter(f):
            if 'istart' in line:
                i = re.findall('\d+',line)
                self.istart = int(i[0])
                self.iend = int(i[1])
                self.jstart = int(i[2])
                self.jend = int(i[3])                 
            #print line
        



if __name__ == "__main__":
    
    # data path
    # datadir='/home/caparmor-work/aponte/qg_nemo/dev/data/'
    datadir='/home1/datawork/aponte/roms_qg0/'
    pref = 'roms_qg0_'

    # read files
    d = qgrun(datadir)
    #print d.x.shape
    #print d.nc.variables['psi'][:].shape
    #print d.nc.variables['t'][:]
    
    # manual input for now
    dt = 1.
    ti = dt*np.arange(1,d.nc.variables['psi'][:].shape[0]+1)
    
    # make of movie of psi    
    tmax=150.
    tsize=10
    lx=d.x/1e3
    ly=d.y/1e3
    x_tks=np.arange(0,1200,300)
    itg=0
    for it,lti in enumerate(ti[np.where(ti<tmax)]):
    
        # clear the figure
        plt.close()
        f, ax = plt.subplots(1,1, figsize=(5,7))
    
        print "time index=", it, "/", ti.size
    
        ### psi
        psi=d.nc.variables['psi'][it,-1,:,:]
    
        lax=ax
        lax.cla()
        toplt=psi[:]
        lvls=20
        #dl=2.5
        #lvls=np.arange(-45,45.+dl,dl)
        #dl=5
        #lvls=np.arange(-100.,100.+dl,dl)
        im = lax.contourf(lx,ly,toplt, lvls, cmap='PuOr')
        lax.set_aspect('equal', 'box','C')
        divider = make_axes_locatable(lax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax, format="%.1f")
        lax.set_title('psi [m2/s], ' + str(lti) +"d",fontsize=tsize)
        lax.grid(True)
        lax.set_ylabel('y [km]')
        lax.set_xlabel('x [km]')
        lax.set_xticks(x_tks)

        #time.sleep(.5) 
    
        # print the figure
        #if lti<tmax[i]:
        plt.savefig('img/fig%04d.tiff'%itg, dpi=300)
        itg+=1
    
    movie=pref+"psi"
    # clean pre-existing movies
    os.system("rm -rf  movies/"+movie+".gif  movies/"+movie+".mp4 >& /dev/null ")
    #os.system("convert -delay 20 img/*.tiff "+pref+"ssh_xi.gif")
    os.system("convert -delay 100 img/*.tiff movies/"+movie+".gif")
    # converts to mp4 for vlc playing
    os.system("ffmpeg -f gif -i movies/"+movie+".gif  movies/"+movie+".mp4")
    



    # make of movie of pv
    k=45
    #tmax=5.
    itg=0
    for it,lti in enumerate(ti[np.where(ti<tmax)]):
    
        # clear the figure
        plt.close()
        f, ax = plt.subplots(1,1, figsize=(5,7))
    
        print "time index=", it, "/", ti.size
    
        ### psi
        q=d.nc.variables['q'][it,k,:,:]
    
        lax=ax
        lax.cla()
        toplt=q[:]
        lvls=20
        #dl=2.5
        #lvls=np.arange(-45,45.+dl,dl)
        #dl=5
        #lvls=np.arange(-100.,100.+dl,dl)
        im = lax.contourf(lx,ly,toplt, lvls, cmap='PuOr')
        lax.set_aspect('equal', 'box','C')
        divider = make_axes_locatable(lax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax, format="%.1f")
        lax.set_title('q [1/s], ' + str(lti) +"d",fontsize=tsize)
        lax.grid(True)
        lax.set_ylabel('y [km]')
        lax.set_xlabel('x [km]')
        lax.set_xticks(x_tks)

        #time.sleep(.5) 
    
        # print the figure
        #if lti<tmax[i]:
        plt.savefig('img/fig%04d.tiff'%itg, dpi=300)
        itg+=1
    
    movie=pref+"q"
    # clean pre-existing movies
    os.system("rm -rf  movies/"+movie+".gif  movies/"+movie+".mp4 >& /dev/null ")
    #os.system("convert -delay 20 img/*.tiff "+pref+"ssh_xi.gif")
    os.system("convert -delay 100 img/*.tiff movies/"+movie+".gif")
    # converts to mp4 for vlc playing
    os.system("ffmpeg -f gif -i movies/"+movie+".gif  movies/"+movie+".mp4")
    




