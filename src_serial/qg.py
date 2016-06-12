#!/usr/bin/python
# -*- encoding: utf8 -*-

from .grid import *
from .solver import *


class qg():
    """ QG object
    """
    
    def __init__(self,
                 grid_uniform_input = {'Lx':1.e7, 'Ly':1.5e7, 'H':4.e3, 
                                       'Nx':100, 'Ny':150, 'Nz':10},
                 grid_lonlat_file = None,
                 N2 = 1e-4
                 ):
        
        # determine if run is parallel
        print 'The libraries assumes the code is run in parallel'
        self.parallel=True
        
        # horizontal and vertical grids
        if grid_lonlat_file is None:
            print 'The grid is uniform'
            self.grid = grid_uniform(**grid_uniform_input)
        elif self.parallel:
            print 'The grid is in lon/lat space, parallel run'
            print 'The grid file is: '+grid_lonlat_file
            self.grid = grid_lonlat_parallel(grid_lonlat_file)
        else: 
            print 'The grid is in lon/lat space, serial run'
            print 'The grid file is: '+grid_lonlat_file
            self.grid = grid_lonlat_serial(grid_lonlat_file)
            
        # vertical stratification
        self.N2 = N2
                
        # initiate pv inversion solver
        if self.parallel:
            self.inv=pvinversion_parallel()
        else:
            self.inv=pvinversion_serial()

        
    def set_psi(self, analytical_psi=True, file_psi=None):
        """ Set psi to a given value
        """
        print 'Set psi to ...'
    
    
    def setup_time_stepping(self):
        """ Setup the time stepping
        Create additional vectors
        """
        pass
        



