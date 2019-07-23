#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Run qgsolver with outputs from idealized ROMS numerical simulations

python run_datarmor.py seamount/pvfpsi_tD_bD roms_seamount.py /home1/datahome/aponte/qgsolver/input/seamount/
#python run_datarmor.py seamount roms_seamount.py /home1/datahome/aponte/qgsolver/input/seamount/

"""

import time
import sys

from qgsolver.qg import qg_model
from qgsolver.state import add

#
#==================== ROMS case ============================================
#

d2s = 86400.

def roms_input_runs(ncores_x=4, ncores_y=6, ping_mpi_cfg=False):
    ''' Tests with roms configuration (spatially uniform, vertically stretched)
    '''

    #ncores_x=2; ncores_y=4; # desktop
    #ncores_x=32; ncores_y=12; # datarmor

    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        
        # proceeds with computations
        start_time = time.time()
        cur_time = start_time
        
        # Top and Bottom boundary condition type: 'N' for Neumann, 'D' for Dirichlet
        #bdy_type = {'top': 'N_PSI', 'bottom': 'N_PSI'}
        #bdy_type = {'top': 'N_PSI', 'bottom': 'D'}
        #bdy_type = {'top': 'D', 'bottom': 'N_PSI'} # dirichlet on top
        bdy_type = {'top': 'D', 'bottom': 'D'} # dirichlet on top/bottom        

        # vertical subdomain
        vdom = {'Nz': 48}

        # horizontal subdomain
        hdom = {'Nx': 100, 'Ny': 150}

        datapath = '../input/'
        outdir = '../output/'
        hgrid = datapath+'roms_metrics.nc'
        vgrid = datapath+'roms_metrics.nc'
        file_q = datapath+'roms_pv.nc'
        file_psi = datapath+'roms_psi.nc'
        file_bg = datapath+'roms_bg.nc'

        params = {'ncores_x': ncores_x, 'ncores_y' :ncores_y,
                  'hgrid': hgrid, 'vgrid': vgrid, 'vdom': vdom, 'hdom': hdom, 
                  'mask': True, 'f0N2_file': file_q,
                  'K': 2000.e0, 'dt': 0.02 * d2s,
                  'verbose': 1}

        ##               
        qg = qg_model(boundary_types=bdy_type, **params)

        # start filling in variables
        qg.set_q(file=file_q)
        qg.set_psi(file=file_psi)   
        qg.write_state(filename=outdir+'input.nc')
        
        # substract background state
        bstate = qg.set_bstate(file=file_bg)
        add(qg.state, bstate, da=None, a2=-1.)
        qg.write_state(filename=outdir+'input.nc', append=True)

        # reset PV from psi
        #qg.pvinv.q_from_psi(qg.state.Q, qg.state.PSI)

        # make copy of the original state
        state0 = qg.state.copy(qg.da)
        
        ## PV inversion

        qg.invert_pv()
        qg.write_state(filename=outdir+'output_full.nc')
        
        # set lateral boundary conditions to 0
        # this does not work as top and bottom boundary conditions are also using the PSI provided
        #qg.state = state0.copy(qg.da)
        #qg.invert_pv(PSI=qg.state.PSI*0)
        #qg.write_state(filename=outdir+'output_full_lat0.nc')
        
        ## PV inversion without PV
    
        qg.state = state0.copy(qg.da)
        qg.state.Q = qg.state.Q*0
        qg.invert_pv()
        qg.write_state(filename=outdir+'output_bsqg.nc')

        # set lateral boundary conditions to 0
        # this does not work as top and bottom boundary conditions are also using the PSI provided
        #qg.state = state0.copy(qg.da)
        #qg.state.Q = qg.state.Q*0
        #qg.invert_pv(PSI=qg.state.PSI*0)
        #qg.write_state(filename=outdir+'output_bsqg_lat0.nc')
        
        if qg._verbose>0:
            print('----------------------------------------------------')
            print('Elapsed time for all ',str(time.time() - cur_time))
                         
        return qg




def main(ping_mpi_cfg=False):    
    
    #
    qg = roms_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    # 
    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print('All done \n')


if __name__ == "__main__":
    main()
    
    
