! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To calculate heat and moisture fluxes in the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
MODULE mym_ex_flux_tq_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MYM_EX_FLUX_TQ_MOD'
CONTAINS

SUBROUTINE mym_ex_flux_tq(                                                     &
      bl_levels, nSCMDpkgs, L_SCMDiags,                                        &
      tl, qw, rhokh, rhogamt, rhogamq, rdz,                                    &
      ftl, fqw)

USE atm_fields_bounds_mod, ONLY: tdims, pdims
USE model_domain_mod,      ONLY: model_type, mt_single_column
USE planet_constants_mod,  ONLY: cp, grcp
USE s_scmop_mod,           ONLY: default_streams,                              &
                                 t_avg, d_bl, scmdiag_bl
USE scmoutput_mod,         ONLY: scmoutput

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! INTENT IN Variables
INTEGER, INTENT(IN) ::                                                         &
   bl_levels
                 ! Max. no. of "boundary" levels

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                         &
   nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                         &
   L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

REAL(KIND=real_umphys), INTENT(IN) ::                                          &
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),         &
                   ! Liquid/frozen water temperture (K)
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),         &
                   ! Total water content (kg/kg)
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels),                                                           &
                   ! Exchange coeffs for scalars
                   ! between K and K-1 on theta levels.
                   ! i.e. the coeffs are defined on rho levels
   rhogamt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels),                                                       &
                   ! Counter gradient term for FTL on rho levels
   rhogamq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels),                                                       &
                   ! Counter gradient Term for FQW on
   rdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)
                   ! RDZ(,1) is the reciprocal
                   ! height of level 1, i.e. of the
                   ! middle of layer 1.  For K > 1,
                   ! RDZ(,K) is the reciprocal of the
                   ! vertical distance from level
                   ! K-1 to level K.

! INTENT OUT Variables
REAL(KIND=real_umphys), INTENT(IN OUT) ::                                      &
   ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),         &
                   ! FTL(,K) contains net turb
                   ! sensible heat flux into layer K
                   ! from below; so FTL(,1) is the
                   ! surface sensible heat, H. (W/m2)
                   ! defined on rho levels
   fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, bl_levels)
                   ! Moisture flux between layers
                   ! (kg per square metre per sec).
                   ! FQW(,1) is total water flux
                   ! from surface, 'E'.
                   ! defined on rho levels

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'MYM_EX_FLUX_TQ'

! LOCAL VARIABLES.

INTEGER ::                                                                     &
   i, j, k

REAL(KIND=real_umphys) ::                                                      &
   grad_ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                   ! Gradient part of FTL
                   ! K*dth/dz
   grad_fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                   ! Gradient part of FQW
                   ! K*dq/dz
   count_grad_ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  bl_levels),                                                  &
                   ! Counter gradient part of FTL
   count_grad_fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  bl_levels)
                   ! Counter gradient part of FQW

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO k = 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      grad_ftl(i,j,k)=0.0
      grad_fqw(i,j,k)=0.0
      count_grad_ftl(i,j,k)=0.0
      count_grad_fqw(i,j,k)=0.0
    END DO
  END DO
END DO

DO k = 2, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      grad_ftl(i,j,k)= - rhokh(i,j,k) *                                        &
                ( ( ( tl(i,j,k) - tl(i,j,k-1) ) * rdz(i,j,k) )                 &
                                                        + grcp )
      grad_fqw(i,j,k)= - rhokh(i,j,k) *                                        &
                    ( qw(i,j,k) - qw(i,j,k-1) ) * rdz(i,j,k)
      count_grad_ftl(i,j,k) = -rhogamt(i,j,k)
      count_grad_fqw(i,j,k) = -rhogamq(i,j,k)
      ftl(i,j,k) = grad_ftl(i,j,k) + count_grad_ftl(i,j,k)
      fqw(i,j,k) = grad_fqw(i,j,k) + count_grad_fqw(i,j,k)
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
IF ( l_scmdiags(scmdiag_bl) .AND.                                              &
     model_type == mt_single_column ) THEN

  DO k=1, bl_levels
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        grad_ftl(i,j,k)= cp * grad_ftl(i,j,k)
        count_grad_ftl(i,j,k) = cp * count_grad_ftl(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(grad_ftl,'Grad_ftl',                                          &
       'Down gradient flux of TL','W/m2',                                      &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(count_grad_ftl,'CG_ftl',                                      &
       'Counter gradient part of flux of TL','W/m2',                           &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(grad_fqw,'Grad_fqw',                                          &
       'Down-gradient flux of QW','kg/m2/s',                                   &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(count_grad_fqw,'CG_fqw',                                      &
       'Counter gradient part of  flux of QW','kg/m2/s',                       &
       t_avg,d_bl,default_streams,'',routinename)

END IF ! scmdiag_bl / model_type

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mym_ex_flux_tq
END MODULE mym_ex_flux_tq_mod
