#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

import time
import sys

#try:
#    from qgsolver.qg import qg_model
#    from qgsolver.io import write_nc
#except:
#    print 'qgsolver not yet in path'
#    sys.exit

from qgsolver.qg import qg_model
from qgsolver.io import write_nc


#
#==================== Synthetic curvilinear case ============================================
#



def curvilinear_runs(ncores_x=8, ncores_y=8, ping_mpi_cfg=False):
    ''' Tests with curvilinear grid
    ''' 
    
    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        # proceeds with computations            
    
        qg = qg_model(hgrid = 'curv_metrics.nc', vgrid = 'curv_metrics.nc',
                      f0N2_file = 'curv_pv.nc',
                      K = 1.e3, dt = 0.5*86400.e0)
        qg.case='curv'
        #
        qg.set_q(file_q='curv_pv.nc')
        qg.invert_pv()
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg)
        
        test=1
        if test==0:
            # one time step and store
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg, create=False)
        elif test==1:
            while qg.tstepper.t/86400. < 200 :
                qg.tstep(1)
                write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg, create=False)
    
        return qg

def main(ping_mpi_cfg=False):    

    print('This script needs to updated')
    sys.exit()

    # qg = uniform_grid_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    #qg = roms_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    qg = nemo_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    # qg = test_L(ping_mpi_cfg=ping_mpi_cfg)
    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print 'Test done \n'


if __name__ == "__main__":
    main()
    
    
