#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

#import qgsolver.qg as qg
from qgsolver.qg import qg_model
from qgsolver.io import write_nc


def uniform_grid_runs():
    ''' Tests with uniform grid, closed domains
    '''
    qg = qg_model(hgrid = {'Nx':150, 'Ny':100}, vgrid = {'Nz':3 },
            K = 0.e0, dt = 0.5*86400.e0)
    #
    qg.set_q()
    qg.invert_pv()
    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg)
    #
    test=2
    if test==0:
        # one time step and store
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)
    elif test==1:
        # write/read/write
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output1.nc', qg, create=True)
        qg.set_q(file_q='output.nc')
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output1.nc', qg, create=False)
    else:
        while qg.tstepper.t/86400. < 200 :
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)
    
    return qg


def curvilinear_runs():
    ''' Tests with curvilinear grid
    ''' 
    
    qg = qg_model(hgrid = 'curv_metrics.nc', vgrid = 'curv_metrics.nc',
                  f0N2_file = 'curv_pv.nc',
                  K = 0.e0, dt = 0.5*86400.e0)
    #
    qg.set_q(file_q='curv_pv.nc')
    qg.invert_pv()
    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg)
    
    #
    #===========================================================================
    # test=2
    # if test==0:
    #     # one time step and store
    #     qg.tstep(1)
    #     write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)
    # elif test==1:
    #     # write/read/write
    #     qg.tstep(1)
    #     write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output1.nc', qg, create=True)
    #     qg.set_q(file_q='output.nc')
    #     qg.tstep(1)
    #     write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output1.nc', qg, create=False)
    # else:
    #     while qg.tstepper.t/86400. < 200 :
    #         qg.tstep(1)
    #         write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)
    #===========================================================================

    return qg
    

if __name__ == "__main__":
    
    #qg = uniform_grid_runs()
    
    qg = curvilinear_runs()
    

    if qg._verbose:
        print 'Test done \n'
