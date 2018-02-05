#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Read result from a PV inversion and time stepping
"""

import sys, os
import xarray as xr
import matplotlib.pyplot as plt
import re

g=9.1

class qgrun():
    ''' Container for now of useful info
    '''
    def __init__(self, path):
        
        #
        self.path = path
        print('Load run from following path: '+path)
        
        # load log and get useful info
        #self.read_log()
        
        # load grid
        self.m = xr.open_dataset(path+'input/roms_metrics.nc')
        _ds = xr.open_dataset(path+'input/roms_pv.nc')
        self.f0 = _ds['f0']
       
        # load data
        #self.nc = MFDataset(path+'output/output*.nc')
        self.bg = xr.open_dataset(path+'input/roms_bg.nc')
        self.bg.rename({'zt': 'z'}, inplace=True)
        self.ds = xr.open_mfdataset(path+'output/output_*.nc', concat_dim='t', compat='identical')
        self.t = self.ds['t']
        (self.Nt, self.N, self.M, self.L) = self.ds['psi'].shape
        
        # add back background state
        self.ds['psi'] += self.bg['psi']
        self.ds['q'] += self.bg['q']

    def __getitem__(self,key):
        if key in self.ds:
            return self.ds[key]
        elif key in self.m:
            return self.m[key]
        
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
    rdir = '/home1/datawork/aponte/qgsolver/'
    datadir = rdir+'roms_out_200/'
    pref = datadir.split('/')[-2]+'_'

    # read files
    d = qgrun(datadir)
    d.ds = d.ds.isel(z=d.N-1)
    #print(d.ds)
    #print(d.bg)
    #print(d.m)
    
    for it,t in enumerate(d.t):
        # clear the figure
        plt.close()
        f, ax = plt.subplots(1,1, figsize=(5,7))
        toplt = d.ds['psi'].sel(t=t)*d.f0/g
        toplt.plot(ax=ax,vmin=-.5,vmax=.5)
        toplt.plot.contour(ax=ax,levels=[0.],color='w')
        ax.set_aspect('equal')
        ax.set_title('t = %.1f' %t)
        plt.tight_layout()
        plt.savefig('img/fig%04d.png'%it, dpi=300)
        print('%d/%d'%(it,d.Nt))
    
    movie=pref+"psi"
    # clean pre-existing movies
    os.system("rm -rf  movies/"+movie+".gif  movies/"+movie+".mp4 >& /dev/null ")
    com = "ffmpeg -y -r 4 -i img/fig*.png  movies/"+movie+".mp4"
    print(com)
    os.system("ffmpeg -y -r 4 -i img/fig*.png  movies/"+movie+".mp4")

