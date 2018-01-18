! ---------------------------------------------------------------------
!
!    PV advection
!
! --------------------------------------------------------------------

#include <petsc/finclude/petscdm.h>

module timestepper_adv

  use petsc
  implicit none

  type gridinfo
     PetscInt mx,xs,xe,xm,gxs,gxe,gxm
     PetscInt my,ys,ye,ym,gys,gye,gym
     PetscInt mz,zs,ze,zm,gzs,gze,gzm
  end type gridinfo

contains

  subroutine GetGridInfo(da, grd, ierr)
    implicit none
    DM            da
    type(gridinfo) grd
    PetscErrorCode ierr
    !
    call DMDAGetInfo(da, PETSC_NULL_INTEGER, &
         &           grd%mx, grd%my, grd%mz, &
         &           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         &           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         &           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
         &           PETSC_NULL_INTEGER,ierr); CHKERRQ(ierr)
    call DMDAGetCorners(da, &
         &              grd%xs,grd%ys,grd%zs, &
         &              grd%xm,grd%ym,grd%zm,ierr); CHKERRQ(ierr)
    call DMDAGetGhostCorners(da, &
         &                   grd%gxs,grd%gys,grd%gzs, &
         &                   grd%gxm,grd%gym,grd%gzm,ierr); CHKERRQ(ierr)

    grd%xs  = grd%xs+1
    grd%ys  = grd%ys+1
    grd%zs  = grd%zs+1
    grd%gxs = grd%gxs+1
    grd%gys = grd%gys+1
    grd%gzs = grd%gzs+1

    grd%xe  = grd%xs+grd%xm-1
    grd%ye  = grd%ys+grd%ym-1
    grd%ze  = grd%zs+grd%zm-1
    grd%gxe = grd%gxs+grd%gxm-1
    grd%gye = grd%gys+grd%gym-1
    grd%gze = grd%gzs+grd%gzm-1

  end subroutine GetGridInfo

  subroutine uniformLocal(grd, rhs, q, psi, dx, dy, bdytype, ierr)
    implicit none
    type(gridinfo) grd
    PetscScalar    rhs(grd%xs:grd%xe,grd%ys:grd%ye)
    PetscScalar    q(grd%xs:grd%xe,grd%ys:grd%ye)
    PetscScalar    psi(grd%xs:grd%xe,grd%ys:grd%ye)
    PetscReal      dx, dy
    PetscErrorCode ierr
    !
    PetscInt       i, j, k
    PetscReal      sc

    sc = 1.0/dx/dy*0.25/3.0

    do k=grd%ks,grd%ke
       do j=grd%ys,grd%ye
          do i=grd%xs,grd%xe
             if ((i==1 .or. j==1 .or. i==grd%mx .or. j==grd%my) .and. bdytype==0 ) then
                ! boundary points
                rhs(i,j,k) = 0.0
             else
                ! interior grid points
                ! Arakawa Jacobian
                rhs(i,j,k) = rhs(i,j,k) + &
                             (  (q[i+1,j,k] - q[i-1,j,k] ) * ( psi[i,j+1,k] - psi[i,j-1,k] ) &
                                -(q[i,j+1,k] - q[i,j-1,k] ) * ( psi[i+1,j,k] - psi[i-1,j,k] ) &
                              +( q[i+1,j,k] * (psi[i+1,j+1,k]-psi[i+1,j-1,k]) &
                                -q[i-1,j,k] * (psi[i-1,j+1,k]-psi[i-1,j-1,k]) &
                                -q[i,j+1,k] * (psi[i+1,j+1,k]-psi[i-1,j+1,k]) &
                                +q[i,j-1,k] * (psi[i+1,j-1,k]-psi[i-1,j-1,k]) ) &
                              +( q[i+1,j+1,k] * (psi[i,j+1,k]-psi[i+1,j,k]) &
                                -q[i-1,j-1,k] * (psi[i-1,j,k]-psi[i,j-1,k]) &
                                -q[i-1,j+1,k] * (psi[i,j+1,k]-psi[i-1,j,k]) &
                                +q[i+1,j-1,k] * (psi[i+1,j,k]-psi[i,j-1,k]) ) ) *sc
             end if
          end do
       end do
    end do
    ierr = 0

  end subroutine uniformLocal

end module timestepper_adv

! --------------------------------------------------------------------

subroutine uniform(da, RHS, Q, PSI, dx, dy, bdytype, ierr)
  use timestepper_adv
  implicit none
  DM da
  Vec Q
  Vec PSI
  PetscReal dx, dy, dz
  PetscErrorCode ierr
  ! should try to condensate previous declarations
  !
  type(gridinfo)      :: grd
  PetscScalar,pointer :: rhsrhs(:)
  PetscScalar,pointer :: qq(:)
  PetscScalar,pointer :: psipsi(:)

  call VecGetArrayF90(RHS,rhsrhs,ierr); CHKERRQ(ierr)
  call VecGetArrayF90(Q,qq,ierr); CHKERRQ(ierr)
  call VecGetArrayF90(PSI,psipsi,ierr); CHKERRQ(ierr)

  call GetGridInfo(da,grd,ierr); CHKERRQ(ierr)
  call computeADV_uniformLocal(grd,rhsrhs,qq,psipsi,dx,dy,dz,bdytype,ierr); CHKERRQ(ierr)

  call VecRestoreArrayF90(RHS,rhsrhs,ierr); CHKERRQ(ierr)

end subroutine uniform

! --------------------------------------------------------------------


! Local Variables:
! mode: f90
! End: