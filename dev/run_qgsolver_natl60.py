#!/usr/bin/python
# -*- encoding: utf8 -*-

"""  Run QG solver for NATL60 input data
"""

import time
import sys
import copy

from qgsolver.qg import qg_model
from qgsolver.io import write_nc, read_nc_petsc
from qgsolver.solver import pvinversion


# Parameters

datapath = 'data/'
file_q = datapath+'nemo_pvregfrom_sossheig.nc'
file_psi = datapath+'nemo_psi_true.nc'
file_rho = datapath+'nemo_rho_true.nc'
file_psi_bg = datapath+'nemo_psi_tmean.nc' # optional
file_psi_ot = datapath+'nemo_psi0regfrom_sossheig.nc'  # optional

# Boundary condition type: 
#   true boundary conditions: 
#    'N': Neumann (density-like),
#    'D': Dirichlet (streamfunction),
#   other boundary conditions: 
#    'NBG': Neumann-background,
#    'DBG': Dirichlet-background,
#    'NOT': Neumann-other
#    'DOT': Dirichlet-other

bdy_type = {'top':'D', 'bottom':'NBG', 'lateral': 'DOT'}

# LMX domain: Nx=1032, Ny=756, Nz=300

# vertical subdomain
# vdom = {'kdown': 0, 'kup': 50-1, 'k0': 100 }     # small 
vdom = {'kdown': 0, 'kup': 150-1, 'k0': 100 }     # large

# horizontal subdomain
# hdom = {'istart': 0, 'iend': 200-1, 'i0': 50, 'jstart': 0, 'jend': 200-1,  'j0': 50}   # small 
hdom = {'istart': 0, 'iend': 944-1, 'i0': 50, 'jstart': 0, 'jend': 672-1,  'j0': 50}   # large 

# 448=8x56
# 512=8x64

ncores_x=16; ncores_y=14  # large datarmor  

#================================================================

# Boundary conditions

def set_rhs_bdy(qg,fpsi_bg,fpsi_ot):
    """
    Set South/North, East/West, Bottom/Top boundary conditions
    Set RHS along boundaries for inversion, may be an issue
    for time stepping
    :param da: abstract distributed memory object of the domain
    :param qg: qg_model instance
    :return:
    """

    rhs = qg.da.getVecArray(qg.pvinv._RHS)
    mx, my, mz = qg.da.getSizes()
    (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

    istart = qg.grid.istart
    iend = qg.grid.iend
    jstart = qg.grid.jstart
    jend = qg.grid.jend
    kdown = qg.grid.kdown
    kup = qg.grid.kup

    if qg.case == "roms" or 'nemo' or 'uniform':

        psi = qg.da.getVecArray(qg.PSI)
        if fpsi_bg:
            psi_bg = qg.da.getVecArray(qg.PSI_BG)
        if fpsi_ot:
            psi_ot = qg.da.getVecArray(qg.PSI_OT)
        rho = qg.da.getVecArray(qg.RHO)

        if qg.bdy_type['top'] in ['D']:
            psi_up=psi
        if qg.bdy_type['bottom'] in ['D']:
            psi_down=psi
        if qg.bdy_type['lateral'] in ['D']:
            psi_lat=psi

        if qg.bdy_type['top'] in ['N']:
            rho_up=rho
        if qg.bdy_type['bottom'] in ['N']:
            rho_down=rho

        if qg.bdy_type['top'] in ['NBG','DBG']:
            psi_up=psi_bg
        if qg.bdy_type['bottom'] in ['NBG','DBG']:
            psi_down=psi_bg
        if qg.bdy_type['lateral'] in ['DBG']:
            psi_lat=psi_bg

        if qg.bdy_type['top'] in ['NOT','DOT']:
            psi_up=psi_ot
        if qg.bdy_type['bottom'] in ['NOT','DOT']:
            psi_down=psi_ot
        if qg.bdy_type['lateral'] in ['DOT']:
            psi_lat=psi_ot

        # lower ghost area
        if zs < kdown:
            for k in range(zs,kdown):
                for j in range(ys, ye):
                    for i in range(xs, xe):                    
                        # rhs[i,j,k]=sys.float_info.epsilon
                        rhs[i,j,k]=psi_down[i, j, k]
        # bottom bdy
        k = kdown
        if qg.bdy_type['bottom'] in ['N']: 
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = - qg.g*0.5*(rho_down[i, j, k]+rho_down[i, j, k+1])/(qg.rho0*qg.f0)
        elif qg.bdy_type['bottom'] in ['NBG','NOT']: 
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = (psi_down[i,j,k+1]-psi_down[i,j,k])/qg.grid.dzw[k] 
        elif qg.bdy_type['bottom'] in ['D','DBG','DOT']:
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = psi_down[i,j,k]

                    
        else:
            print qg.bdy_type['bottom']+" unknown bottom boundary condition"
            sys.exit()

        # debug: computes vertical bdy from psi
        # for j in range(ys, ye):
        #     for i in range(xs, xe):
        #         rhs[i, j, k] = (psi[i,j,k+1]-psi[i,j,k])/qg.grid.dzw[k] 


        if ze > kup+1:
            for k in range(kup+1,ze):
                for j in range(ys, ye):
                    for i in range(xs, xe):
                        # rhs[i,j,k]=sys.float_info.epsilon   
                        rhs[i,j,k]= psi_up[i,j,k]
                                    
        # upper bdy
        k = kup
        if qg.bdy_type['top'] in ['N']: 
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = - qg.g*0.5*(rho_up[i, j, k]+rho_up[i, j, k-1])/(qg.rho0*qg.f0)
        elif qg.bdy_type['top'] in ['NBG','NOT']: 
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = (psi_up[i,j,k]-psi_up[i,j,k-1])/qg.grid.dzw[k-1] 
        elif qg.bdy_type['top'] in ['D','DBG','DOT']:
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = psi_up[i,j,k]

        else:
            print qg.bdy_type['top']+" unknown top boundary condition"
            sys.exit()

        # south bdy
        if ys <= jstart:
            #j = 0
            for k in range(zs, ze):
                for j in range(ys,min(ye,jstart+1)):
                    for i in range(xs, xe):
                        rhs[i, j, k] = psi_lat[i, j, k]
        # north bdy
        if ye >= jend:
            #j = my - 1
            for k in range(zs, ze):
                for j in range(max(ys,jend),ye):
                    for i in range(xs, xe):
                        rhs[i, j, k] = psi_lat[i, j, k]
        # west bdy
        if xs <= istart:
            #i = 0
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(xs,min(xe,istart+1)):
                        rhs[i, j, k] = psi_lat[i, j, k]
        # east bdy
        if xe >= iend:
            #i = mx - 1
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(max(xs,iend),xe):
                        rhs[i, j, k] = psi_lat[i, j, k]

        # debug: computes vertical bdy from psi    
        # for j in range(ys, ye):
        #     for i in range(xs, xe):
        #         rhs[i, j, k] = (psi[i,j,k]-psi[i,j,k-1])/qg.grid.dzw[k-1]

    else:

        # bottom bdy
        if zs == 0:
            k = 0
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = 0.
        # upper bdy
        if ze == mz :
            k = mz-1
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = 0.

    if qg.pvinv._verbose>0:
        print 'set RHS along boudaries for inversion '

#================================================================

def set_rhs_mask(qg):
        """
        Set mask on rhs: where mask=0 (land) rhs=psi
        - param da: abstract distributed memory object of the domain
        - param qg: qg_model instance
             qg.grid.D[qg.grid._k_mask]: mask
        - self.rhs : vorticity whith boundary conditions
        return: masked rhs
        """

        rhs = qg.da.getVecArray(qg.pvinv._RHS)
        mask = qg.da.getVecArray(qg.grid.D)
        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        istart = qg.grid.istart
        iend = qg.grid.iend
        jstart = qg.grid.jstart
        jend = qg.grid.jend
        kdown = qg.grid.kdown
        kup = qg.grid.kup
        kmask = qg.grid._k_mask

        if qg.case == "roms" or 'nemo' or 'uniform':

            if qg.bdy_type['lateral'] == 'D':
                psi = qg.da.getVecArray(qg.PSI)
            elif qg.bdy_type['lateral'] == 'DBG':
                psi = qg.da.getVecArray(qg.PSI_BG)
            elif qg.bdy_type['lateral'] == 'DOT':
                psi = qg.da.getVecArray(qg.PSI_OT)

            if qg.rank==0: print qg.bdy_type
            
            # interior
            for k in range(zs,ze):
                for j in range(ys, ye):
                    for i in range(xs, xe):
                        if mask[i,j,kmask]==0.:
                            rhs[i, j, k] = psi[i,j,k]

        if qg.pvinv._verbose>0:
            print 'set RHS mask for inversion '

#================================================================

def read_petsc(qg,fpsi_bg,fpsi_ot,
               file_q,file_rho,file_psi,file_psi_bg=None,file_psi_ot=None):

    cur_time = time.time()
    # set from files
    read_nc_petsc(qg.Q, 'q', file_q, qg, fillmask=0.)
    if qg.rank == 0:
        print '----------------------------------------------------'
        print 'Elapsed time setting Q ', str(time.time() - cur_time)
    cur_time = time.time()
 
    read_nc_petsc(qg.PSI, 'psi', file_psi, qg, fillmask=0.)
    if qg.rank == 0:
        print '----------------------------------------------------'
        print 'Elapsed time setting PSI ', str(time.time() - cur_time)
    cur_time = time.time()
 
    read_nc_petsc(qg.RHO, 'rho', file_rho, qg, fillmask=0.)
    if qg.rank == 0:
        print '----------------------------------------------------'
        print 'Elapsed time setting RHO ', str(time.time() - cur_time)
    cur_time = time.time()

    if fpsi_bg:
        qg.PSI_BG = qg.da.createGlobalVec() 
        read_nc_petsc(qg.PSI_BG, 'psi', file_psi_bg, qg, fillmask=0.)
        if qg.rank == 0:
            print '----------------------------------------------------'
            print 'Elapsed time setting PSI_BG ', str(time.time() - cur_time)
        cur_time = time.time()

    if fpsi_ot:
        qg.PSI_OT = qg.da.createGlobalVec()          
        read_nc_petsc(qg.PSI_OT, 'psi', file_psi_ot, qg, fillmask=0.)
        if qg.rank == 0:
            print '----------------------------------------------------'
            print 'Elapsed time setting PSI_OT ', str(time.time() - cur_time)
        cur_time = time.time()


def pvinv_solver(qg,fpsi_bg,fpsi_ot):
    
    # ONE = qg.set_identity()
    # qg.pvinv.L.mult(ONE,self._RHS)
    # write_nc([self._RHS], ['id'], 'data/identity.nc', qg)

    # qg.PSI.set(0)
    
    # compute L*PSI and store in self._RHS
    qg.pvinv.L.mult(qg.PSI,qg.pvinv._RHS)
    # store L*PSI in netcdf file lpsi.nc
    write_nc([qg.pvinv._RHS], ['rhs'], 'data/lpsiin.nc', qg)
    # copy Q into RHS
    qg.Q.copy(qg.pvinv._RHS)
    if qg.pvinv._substract_fprime:
        # substract f-f0 from PV
        qg.pvinv.substract_fprime_from_rhs(qg)
        if qg.pvinv._verbose>0:
            print 'Substract fprime from pv prior to inversion'
    # fix boundaries
    #self.set_rhs_bdy(qg)
    set_rhs_bdy(qg, fpsi_bg, fpsi_ot)
    # mask rhs 
    set_rhs_mask(qg)
    # store RHS in netcdf file rhs.nc
    write_nc([qg.pvinv._RHS], ['rhs'], 'data/rhs.nc', qg)
    # actually solves the pb
    qg.pvinv.ksp.solve(qg.pvinv._RHS, qg.PSI)
    # compute L*PSI and store in self._RHS
    qg.pvinv.L.mult(qg.PSI,qg.pvinv._RHS)
    # store L*PSI in netcdf file lpsi.nc
    write_nc([qg.pvinv._RHS], ['Lpsi'], 'data/lpsiout.nc', qg)

def nemo_input_runs(ncores_x=ncores_x, ncores_y=ncores_y, ping_mpi_cfg=False):
    
    """ Set up processing parameters and run qg solver.
    """

    # activate bg or other field
    fpsi_bg=False
    fpsi_ot=False
    for key, value in bdy_type.items():
        if value in ['NBG','DBG']:
            fpsi_bg=True
        if value in ['NOT','DOT']:
            fpsi_ot=True
    
    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        # proceeds with computations
        start_time = time.time()
        cur_time = start_time
    
        # MPI decomposition of the domain
        # case must be defined before ncores for run_caparmor.py
        casename = 'nemo'
    
        hgrid = datapath+'nemo_metrics.nc'
        vgrid = datapath+'nemo_metrics.nc'

        # equivalent boundary condition (to initialize qg solver)
        bdy_type_tmp=copy.copy(bdy_type)
        for key, value in bdy_type_tmp.items():
            if value in ['N','NBG','NOT']:
                bdy_type_tmp[key]='N'
            if value in ['D','DBG','DOT']:
                bdy_type_tmp[key]='D'

        qg = qg_model(hgrid=hgrid, vgrid=vgrid, f0N2_file=file_q, K=1.e0, dt=0.5 * 86400.e0,
                      vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y, 
                      bdy_type_in=bdy_type_tmp, substract_fprime=True)
        
        # reset boundary condition to prescribed value
        qg.bdy_type.update(bdy_type)
    
        qg.case=casename

        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for qg_model ', str(time.time() - cur_time)
        cur_time = time.time()

       
#         # read 3D variables (variable is set to None if optional and non-existing)
#         read_nc_3D(qg, ['PSI', 'PSI_BG', 'PSI_OT', 'Q', 'RHO'], 
#                    [file_psi, file_psi_bg, file_psi_ot, file_q, file_rho])

        
        # initialize 3D-variables
        
        # set analytically
#         qg.set_q_analytically()
#         qg.set_rho_analytically()
#         qg.set_psi_analytically()

        read_petsc(qg,fpsi_bg,fpsi_ot,file_q,file_rho,file_psi,file_psi_bg,file_psi_ot)        

        # build the list of variables to write in input.nc
        if not fpsi_bg and not fpsi_ot:
            petsc_writein=[qg.PSI, qg.Q]
            vname_writein=['psi', 'q']            
        elif fpsi_bg and not fpsi_ot:
            petsc_writein=[qg.PSI, qg.PSI_BG, qg.Q]
            vname_writein=['psi', 'psi_bg', 'q']
        elif not fpsi_bg and fpsi_ot:
            petsc_writein=[qg.PSI, qg.PSI_OT, qg.Q]
            vname_writein=['psi', 'psi_ot', 'q']
        elif fpsi_bg and fpsi_ot:
            petsc_writein=[qg.PSI, qg.PSI_BG, qg.PSI_OT, qg.Q]
            vname_writein=['psi', 'psi_bg', 'psi_ot', 'q']
        
        write_nc(petsc_writein, vname_writein, 'data/input.nc', qg)
            
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for write_nc ', str(time.time() - cur_time)
        cur_time = time.time()

        if qg._verbose>1:
            print 'Inversion done'

        #qg.pvinv.solve(qg)
        pvinv_solver(qg,fpsi_bg,fpsi_ot)
        
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for invert_pv ', str(time.time() - cur_time)
        cur_time = time.time()

        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg)
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for write_nc ', str(time.time() - cur_time)
        cur_time = time.time()
    
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time  ', str(cur_time - start_time)
    
        return qg


def main(ping_mpi_cfg=False):    
    
    qg = nemo_input_runs(ping_mpi_cfg=ping_mpi_cfg)
        
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print 'Test done \n'


if __name__ == "__main__":
    main()   
