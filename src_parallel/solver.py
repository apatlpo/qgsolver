#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys

from .set_L import *
from .io import write_nc

#
#==================== Serial solver ============================================
#

class pvinversion():
    """ PV inversion, parallel
    """
    
    def __init__(self, qg):
        """ Setup the PV inversion solver
        """
                
        self._verbose = qg._verbose
                        
        # create the operator
        self.L = qg.da.createMat()
        #
        if self._verbose>0:
            print 'Operator L declared'

        # Fill in operator values
        if qg.grid._flag_hgrid_uniform and qg.grid._flag_vgrid_uniform:
            set_L(self.L, qg)
        else:
            set_L_curv(self.L, qg)

        #
        if self._verbose>0:
            print 'Operator L filled'

        # global vector for PV inversion
        self._RHS = qg.da.createGlobalVec()

        # local vectors
        #self._localRHS  = qg.da.createLocalVec()
        #self._localPSI  = qg.da.createLocalVec()

        # create solver
        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_WORLD)
        self.ksp.setOperators(self.L)
        # self.ksp.setType('cg')
        # self.ksp.setType('gmres')
        self.ksp.setType('bicg')
        self.ksp.setInitialGuessNonzero(True)
        # and incomplete Cholesky for preconditionning
        # self.ksp.getPC().setType('icc')
        # self.ksp.getPC().setType('bjacobi')
        # self.ksp.getPC().setType('asm')
        #self.ksp.getPC().setType('none')
        # set tolerances
        # self.ksp.setTolerances(rtol=1e-10) # nope
        self.ksp.setTolerances(max_it=1000)
        #
        #
        for opt in sys.argv[1:]:
	        PETSc.Options().setValue(opt, None)
        self.ksp.setFromOptions()
        
        if self._verbose>0:
            print 'PV inversion is set up'
            
            

    def solve(self, qg):
        """ Compute the PV inversion
        """
        # qg.pvinv.L.mult(qg.PSI,self._RHS)
        # write_nc([self._RHS], ['Lpsi'], 'data/lpsi.nc', qg)
        # copy Q into RHS
        qg.Q.copy(self._RHS)
        # substract f-f0 from PV
        self.substract_fprime_from_rhs(qg)
        # fix boundaries
        self.set_rhs_bdy(qg)
        #write_nc([self._RHS], ['rhs'], 'data/rhs.nc', qg)
        # qg.PSI.set(0)
        # actually solves the pb
        self.ksp.solve(self._RHS, qg.PSI)
        #
        # debug test:
        #self.L.mult(qg.PSI, self._RHS)
        #
        if self._verbose>1:
            print 'Inversion done'
            
    def substract_fprime_from_rhs(self, qg):
        """
        Substract f'=f-f0 from the rhs used in PV inversion
        """
        rhs = qg.da.getVecArray(self._RHS)
        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        D = qg.da.getVecArray(qg.grid.D)

        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):                    
                    rhs[i,j,k] -= D[i,j,qg.grid._k_f] - qg.f0
        
        if self._verbose>0:
            print 'Substract f-f0 from pv prior to inversion'

        

    def set_rhs_bdy(self, qg):
        """
        Set South/North, East/West, Bottom/Top boundary conditions
        Set RHS along boundaries for inversion, may be an issue
        for time stepping
        :param da: abstract distributed memory object of the domain
        :param qg: qg_model instance
        :return:
        """

        rhs = qg.da.getVecArray(self._RHS)
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
            rho = qg.da.getVecArray(qg.RHO)

            # lower ghost area
            if zs < kdown:
                for k in range(zs,kdown):
                    for j in range(ys, ye):
                        for i in range(xs, xe):                    
                            rhs[i,j,k]=sys.float_info.epsilon
            # bottom bdy
            k = kdown
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = - qg.g*rho[i, j, k]/(qg.rho0*qg.f0)
            
            # upper ghost area (!!! ze=Nz and not Nz-1)
            if ze > kup+1:
                for k in range(kup+1,ze):
                    for j in range(ys, ye):
                        for i in range(xs, xe):
                            rhs[i,j,k]=sys.float_info.epsilon                
            # upper bdy
            k = kup
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = - qg.g*rho[i, j, k]/(qg.rho0*qg.f0)
            
            # debug: computes vertical bdy from psi
            # bottom bdy
            # if zs == 0:
            #     k = 0
            #     for j in range(ys, ye):
            #         for i in range(xs, xe):
            #             rhs[i, j, k] = (psi[i,j,k+1]-psi[i,j,k])/qg.grid.dzc[k]
            # # upper bdy
            # if ze == mz:
            #     k = mz-1
            #     for j in range(ys, ye):
            #         for i in range(xs, xe):
            #             rhs[i, j, k] = (psi[i,j,k]-psi[i,j,k-1])/qg.grid.dzc[k-1]

            # south bdy
            if ys <= jstart:
                #j = 0
                for k in range(zs, ze):
                    for j in range(ys,min(ye,jstart+1)):
                        for i in range(xs, xe):
                            rhs[i, j, k] = psi[i, j, k]
            # north bdy
            if ye >= jend:
                #j = my - 1
                for k in range(zs, ze):
                    for j in range(max(ys,jend),ye):
                        for i in range(xs, xe):
                            rhs[i, j, k] = psi[i, j, k]
            # west bdy
            if xs <= istart:
                #i = 0
                for k in range(zs, ze):
                    for j in range(ys, ye):
                        for i in range(xs,min(xe,istart+1)):
                            rhs[i, j, k] = psi[i, j, k]
            # east bdy
            if xe >= iend:
                #i = mx - 1
                for k in range(zs, ze):
                    for j in range(ys, ye):
                        for i in range(max(xs,iend),xe):
                            rhs[i, j, k] = psi[i, j, k]

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

        if self._verbose>0:
            print 'set RHS along boudaries for inversion '


#
#==================== Time stepper ============================================
#


class time_stepper():
    """ Time stepper, parallel with petsc4py
    4 steps explicit RungeKutta
    """
    
    def __init__(self, qg, dt, t0 = 0.):
        
        self._verbose = qg._verbose
        
        ### physical parameters
        # laplacian parameter, should move outside of here
        self.K = qg.K
        
        ### time variables
        self.dt = dt
        self._t0 = t0
        self.t = t0
        #print 't = %e d' % (self.t/86400.)
        
        ### 4 steps explicit RungeKutta parameters
        self._a = [1./6., 1./3., 1./3., 1./6.]
        self._b = [0.5, 0.5, 1.]

        ### additional global vectors
        self._RHS0 = qg.da.createGlobalVec()
        self._RHS1 = qg.da.createGlobalVec()
        self._dRHS = qg.da.createGlobalVec()
        
        if self._verbose>0:
            print 'PV time stepper is set up'


    def go(self, qg, nt):
        """ Carry out the time stepping
        """
        _tstep=0
        #for i in xrange(nt):
        while _tstep < nt:
            # update time parameters and indexes
            self.t += self.dt
            _tstep += 1
            #
            qg.Q.copy(self._RHS0) # copies Q into RHS0
            qg.Q.copy(self._RHS1) # copies Q into RHS1
            for rk in range(4):
                if qg.grid._flag_hgrid_uniform and qg.grid._flag_vgrid_uniform:
                    self._computeRHS(qg)
                else:
                    self._computeRHS_curv(qg)
                if rk < 3: qg.Q.waxpy(self._b[rk]*self.dt, self._dRHS, self._RHS0)
                self._RHS1.axpy(self._a[rk]*self.dt, self._dRHS)
            self._RHS1.copy(qg.Q) # copies RHS1 into Q
            # reset q at boundaries
            self.set_rhs_bdy(qg)
            if self._verbose>0:
                print 't = %f d' % (self.t/86400.)
#         if self._verbose>0:
#             print 'Time stepping done'


    
    def _computeRHS(self,qg):
        """ Compute the RHS of the pv evolution equation i.e: J(psi,q)
        Jacobian 9 points (from Q-GCM):
        Arakawa and Lamb 1981:
        DOI: http://dx.doi.org/10.1175/1520-0493(1981)109<0018:APEAEC>2.0.CO;2
        """
    
        ### compute PV inversion to streamfunction
        qg.invert_pv()
        
        ### declare local vectors
        local_RHS  = qg.da.createLocalVec()
        #local_dRHS  = qg.da.createLocalVec()
        local_PSI  = qg.da.createLocalVec()

        ###
        qg.da.globalToLocal(qg.Q, local_RHS)
        qg.da.globalToLocal(qg.PSI, local_PSI)
        #qg.da.globalToLocal(self._dRHS, local_dRHS)
        #
        q = qg.da.getVecArray(local_RHS)
        psi = qg.da.getVecArray(local_PSI)
        dq = qg.da.getVecArray(self._dRHS)
        #dq = qg.da.getVecArray(local_dRHS)
        #
        mx, my, mz = qg.da.getSizes()
        dx, dy, dz = qg.grid.dx, qg.grid.dy, qg.grid.dz
        idx, idy, idz = [1.0/dl for dl in [dx, dy, dz]]
        idx2, idy2, idz2 = [1.0/dl**2 for dl in [dx, dy, dz]]
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
        #
        ### advect PV:
        # RHS= -u x dq/dx - v x dq/dy = -J(psi,q) = - (-dpsi/dy x dq/dx + dpsi/dx x dq/dy)
        #
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if (i==0    or j==0 or
                        i==mx-1 or j==my-1):
                        # lateral boundaries
                        dq[i, j, k] = 0.
                    else:
                        ### Jacobian
                        #
                        # naive approach leads to noodling (see Arakawa 1966, J_pp)
                        # dq[i, j, k] = - ( -dpsidy * dqdx + dpsidx * dqdy)
                        # 
                        # Arakawa Jacobian
                        #
                        J_pp =   ( q[i+1,j,k] - q[i-1,j,k] ) * ( psi[i,j+1,k] - psi[i,j-1,k] ) \
                               - ( q[i,j+1,k] - q[i,j-1,k] ) * ( psi[i+1,j,k] - psi[i-1,j,k] )
                        J_pp *= idx*idy*0.25
                        #
                        J_pc =   q[i+1,j,k] * (psi[i+1,j+1,k]-psi[i+1,j-1,k]) \
                               - q[i-1,j,k] * (psi[i-1,j+1,k]-psi[i-1,j-1,k]) \
                               - q[i,j+1,k] * (psi[i+1,j+1,k]-psi[i-1,j+1,k]) \
                               + q[i,j-1,k] * (psi[i+1,j-1,k]-psi[i-1,j-1,k])
                        J_pc *= idx*idy*0.25
                        #
                        J_cp =   q[i+1,j+1,k] * (psi[i,j+1,k]-psi[i+1,j,k]) \
                               - q[i-1,j-1,k] * (psi[i-1,j,k]-psi[i,j-1,k]) \
                               - q[i-1,j+1,k] * (psi[i,j+1,k]-psi[i-1,j,k]) \
                               + q[i+1,j-1,k] * (psi[i+1,j,k]-psi[i,j-1,k])
                        J_cp *= idx*idy*0.25
                        #
                        dq[i, j, k] = ( J_pp + J_pc + J_cp )/3.
                        #
                        ### Dissipation
                        dq[i, j, k] +=   self.K*(q[i+1,j,k]-2.*q[i,j,k]+q[i-1,j,k])*idx2 \
                                       + self.K*(q[i,j+1,k]-2.*q[i,j,k]+q[i,j-1,k])*idy2  



   
    def _computeRHS_curv(self,qg):
        """ Compute the RHS of the pv evolution equation i.e: J(psi,q)
        Jacobian 9 points (from Q-GCM):
        Arakawa and Lamb 1981:
        DOI: http://dx.doi.org/10.1175/1520-0493(1981)109<0018:APEAEC>2.0.CO;2
        """
    
        ### compute PV inversion to streamfunction
        qg.invert_pv()
        
        ### declare local vectors
        local_RHS  = qg.da.createLocalVec()
        #local_dRHS  = qg.da.createLocalVec()
        local_PSI  = qg.da.createLocalVec()
        
        ###
        qg.da.globalToLocal(qg.Q, local_RHS)
        qg.da.globalToLocal(qg.PSI, local_PSI)
        #qg.da.globalToLocal(self._dRHS, local_dRHS)
        #
        q = qg.da.getVecArray(local_RHS)
        psi = qg.da.getVecArray(local_PSI)
        dq = qg.da.getVecArray(self._dRHS)
        #dq = qg.da.getVecArray(local_dRHS)
        #
        mx, my, mz = qg.da.getSizes()
        #
        #D = qg.da.getVecArray(qg.grid.D)
        local_D  = qg.da.createLocalVec()
        qg.da.globalToLocal(qg.grid.D, local_D)
        D = qg.da.getVecArray(local_D)
        kdx = qg.grid._k_dx
        kdy = qg.grid._k_dy
        kf = qg.grid._k_f
        #
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
        #
        # advect PV:
        # RHS= -u x dq/dx - v x dq/dy = -J(psi,q) = - (-dpsi/dy x dq/dx + dpsi/dx x dq/dy)
        # 
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if (i==0    or j==0 or
                        i==mx-1 or j==my-1):
                        # lateral boundaries
                        dq[i, j, k] = 0.
                    else:
                        ### Jacobian
                        #
                        # naive approach leads to noodling (see Arakawa 1966, J_pp)
                        # dq[i, j, k] = - ( -dpsidy * dqdx + dpsidx * dqdy)
                        # 
                        # Arakawa Jacobian
                        #
                        J_pp =   ( q[i+1,j,k] - q[i-1,j,k] ) * ( psi[i,j+1,k] - psi[i,j-1,k] ) \
                               - ( q[i,j+1,k] - q[i,j-1,k] ) * ( psi[i+1,j,k] - psi[i-1,j,k] )
                        J_pp *= 0.25
                        #
                        J_pc =   q[i+1,j,k] * (psi[i+1,j+1,k]-psi[i+1,j-1,k]) \
                               - q[i-1,j,k] * (psi[i-1,j+1,k]-psi[i-1,j-1,k]) \
                               - q[i,j+1,k] * (psi[i+1,j+1,k]-psi[i-1,j+1,k]) \
                               + q[i,j-1,k] * (psi[i+1,j-1,k]-psi[i-1,j-1,k])
                        J_pc *= 0.25
                        #
                        J_cp =   q[i+1,j+1,k] * (psi[i,j+1,k]-psi[i+1,j,k]) \
                               - q[i-1,j-1,k] * (psi[i-1,j,k]-psi[i,j-1,k]) \
                               - q[i-1,j+1,k] * (psi[i,j+1,k]-psi[i-1,j,k]) \
                               + q[i+1,j-1,k] * (psi[i+1,j,k]-psi[i,j-1,k])
                        J_cp *= 0.25
                        #
                        dq[i, j, k] = ( J_pp + J_pc + J_cp )/3. /D[i,j,kdx]/D[i,j,kdy]
                        #
                        ### Add advection of planetary vorticity, shouldn't f-f0 be part of q though !!!
                        Jp_pp =   ( D[i+1,j,kf] - D[i-1,j,kf] ) * ( psi[i,j+1,k] - psi[i,j-1,k] ) \
                               - ( D[i,j+1,kf] - D[i,j-1,kf] ) * ( psi[i+1,j,k] - psi[i-1,j,k] )
                        Jp_pp *= 0.25
                        #
                        Jp_pc =  D[i+1,j,kf] * (psi[i+1,j+1,k]-psi[i+1,j-1,k]) \
                               - D[i-1,j,kf] * (psi[i-1,j+1,k]-psi[i-1,j-1,k]) \
                               - D[i,j+1,kf] * (psi[i+1,j+1,k]-psi[i-1,j+1,k]) \
                               + D[i,j-1,kf] * (psi[i+1,j-1,k]-psi[i-1,j-1,k])
                        Jp_pc *= 0.25
                        #
                        Jp_cp =  D[i+1,j+1,kf] * (psi[i,j+1,k]-psi[i+1,j,k]) \
                               - D[i-1,j-1,kf] * (psi[i-1,j,k]-psi[i,j-1,k]) \
                               - D[i-1,j+1,kf] * (psi[i,j+1,k]-psi[i-1,j,k]) \
                               + D[i+1,j-1,kf] * (psi[i+1,j,k]-psi[i,j-1,k])
                        Jp_cp *= 0.25
                        #
                        dq[i, j, k] += ( Jp_pp + Jp_pc + Jp_cp )/3. /D[i,j,kdx]/D[i,j,kdy]
                        ### Dissipation
                        dq[i, j, k] +=   self.K*(q[i+1,j,k]-2.*q[i,j,k]+q[i-1,j,k])/D[i,j,kdx]/D[i,j,kdx] \
                                       + self.K*(q[i,j+1,k]-2.*q[i,j,k]+q[i,j-1,k])/D[i,j,kdy]/D[i,j,kdy]


    def set_rhs_bdy(self, qg):
        """ Reset rhs at boundaries such that drhs/dn=0 """
        #
        rhs = qg.da.getVecArray(qg.Q)

        #
        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        # south bdy
        if (ys == 0):
            j = 0
            for k in range(zs, ze):
                for i in range(xs, xe):
                    rhs[i, j, k] = rhs[i, j + 1, k]
        # north bdy
        if (ye == my):
            j = my - 1
            for k in range(zs, ze):
                for i in range(xs, xe):
                    rhs[i, j, k] = rhs[i, j - 1, k]
        # west bdy
        if (xs == 0):
            i = 0
            for k in range(zs, ze):
                for j in range(ys, ye):
                    rhs[i, j, k] = rhs[i + 1, j, k]
        # east bdy
        if (xe == mx):
            i = mx - 1
            for k in range(zs, ze):
                for j in range(ys, ye):
                    rhs[i, j, k] = rhs[i - 1, j, k]
        return