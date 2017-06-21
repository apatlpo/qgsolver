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
from qgsolver.inout import write_nc

import cProfile, pstats, StringIO

#
#==================== ROMS case ============================================
#



def roms_input_runs(ncores_x=2, ncores_y=4, ping_mpi_cfg=False):
    ''' Tests with roms configuration (spatially uniform, vertically stretched)
    '''

    # ncores_x=2; ncores_y=4; # desktop
    # ncores_x=4; ncores_y=3; # desktop
    # ncores_x=4; ncores_y=12; # datarmor
    # ncores_x=8; ncores_y=12; # datarmor
    ncores_x=32; ncores_y=12; # datarmor

    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        # proceeds with computations
    
        # start_time = time.time()
        # cur_time = start_time
    
        # MPI decomposition of the domain
        # case must be defined before ncores for run_caparmor.py
        casename='roms'
        
        # Top and Bottom boundary condition type: 'N' for Neumann, 'D' for Dirichlet
        #bdy_type = {'top':'N', 'bottom':'N'}
        # bdy_type = {'top':'D', 'bottom':'D', 'periodic':True}
        # bdy_type = {'top':'D', 'bottom':'N', 'periodic':True}
        bdy_type = {'top':'N', 'bottom':'N', 'periodic':True}

        # vertical subdomain
        #vdom = {'kdown': 0, 'kup': 49, 'k0': 0 }
        # vdom = {'kdown': 0, 'kup': 30, 'k0': 10 }
        #vdom = {'Nz': 30}
        vdom = {'Nz': 50}

        # horizontal subdomain
        # hdom = {'istart': 0, 'iend': 255, 'i0': 0, 'jstart':0, 'jend':721,  'j0': 0}
        # hdom = {'Nx': 256, 'j0':300, 'Ny':200}
        hdom = {'Nx': 256, 'Ny': 720}
        # 256 = 2^8
        # 720 = 2^4 x 3^2 x 5


        datapath = '../input/'
        outdir = '../output/'
        hgrid = datapath+'roms_metrics.nc'
        vgrid = datapath+'roms_metrics.nc'
        file_q = datapath+'roms_pv.nc'
        file_psi = datapath+'roms_psi.nc'
        file_rho = datapath+'roms_rho.nc'
        qg = qg_model(hgrid=hgrid, vgrid=vgrid, f0N2_file=file_q, K=1.e0, dt=0.01 * 86400.e0,
                      vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y, 
                      bdy_type_in=bdy_type, substract_fprime=True, verbose=2)        
        qg.case=casename
    
        qg.set_q(file_q=file_q)
    
        qg.set_psi(file_psi=file_psi)
    
        qg.set_rho(file_rho=file_rho)
   
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], outdir+'input.nc', qg)
    
        qg.invert_pv()

        qg.get_uv()
        write_nc([qg._U, qg._V, qg.RHO, qg.PSI, qg.Q], ['u', 'v', 'rho', 'psi', 'q'], outdir+'output0.nc', qg)

        # compute CFL
        CFL = qg.compute_CFL()
        if qg._verbose>0: print "CFL="+str(CFL)

        test=0
        if test==0:
            #for it in range(100):
            for it in range(10):
                qg.tstep(10)
                #qg.tstep(5)
                qg.update_rho()
                #qg.get_uv()
                write_nc([qg.RHO, qg.PSI, qg.Q], ['rho', 'psi', 'q'], outdir+'output'+str(it+1)+'.nc',qg)
                #write_nc([qg._U, qg._V, qg.RHO, qg.PSI, qg.Q], ['u', 'v', 'rho', 'psi', 'q'], outdir+'output'+str(it+1)+'.nc',qg)
        elif test==1:
            # write/read/write
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], outdir+'output1.nc', qg, create=True)
            qg.set_q(file_q='data/output.nc')
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], outdir+'output1.nc', qg, create=False)
        elif test==2:
            while qg.tstepper.t/86400. < 200 :
                qg.tstep(1)
                write_nc([qg.PSI, qg.Q], ['psi', 'q'], outdir+'output.nc', qg, create=False)
             

        # if qg.rank == 0: print '----------------------------------------------------'
        # if qg.rank == 0: print 'Elapsed time  ', str(cur_time - start_time)
    
        return qg




def main(ping_mpi_cfg=False):    
    
    # turn profiling on
    flagProfile=True
    
    flagProfile = flagProfile and not ping_mpi_cfg
    if flagProfile:
        pr = cProfile.Profile()
        pr.enable()
        # stream = open('profile.log', 'w');
    #
    qg = roms_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    # 
        
    if flagProfile and qg._verbose>0:
        pr.disable()
        s = StringIO.StringIO()
        # s = stream
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        # ps.print_stats()
        #print s.getvalue()        
        ps.dump_stats('profile.log')

    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print 'Test done \n'


if __name__ == "__main__":
    main()
    
    
