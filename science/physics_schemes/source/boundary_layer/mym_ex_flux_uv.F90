! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To calculate momentum fluxes in the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
MODULE mym_ex_flux_uv_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MYM_EX_FLUX_UV_MOD'
CONTAINS

SUBROUTINE mym_ex_flux_uv(                                                     &
        dimsi, dimsi_s, dimso, bl_levels,                                      &
        rdz_u_v, rhokm_u_v, rhogamuv_uv, u_v, tau_xy_fd_uv,                    &
        tau_x_y, tau_grad, tau_count_grad)

USE atm_fields_bounds_mod, ONLY: array_dims
USE jules_surface_mod, ONLY: formdrag, explicit_stress
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Intent IN Variables
TYPE(array_dims), INTENT(IN) ::                                                &
   dimsi,      & ! Array dimensions for the inputs
   dimsi_s,    & ! Array dimensions for input u or v (has haloes).
   dimso         ! Array dimensions for the outputs and work variables

INTEGER, INTENT(IN) :: bl_levels
                 ! Max. no. of "boundary" levels

REAL(KIND=real_umphys), INTENT(IN) ::                                          &
   rdz_u_v (dimsi%i_start:dimsi%i_end,                                         &
            dimsi%j_start:dimsi%j_end, 2:bl_levels),                           &
                 ! Reciprocal of the vertical
                 ! distance from level K-1 to
                 ! level K. (K > 1) on wind levels
   rhokm_u_v (dimsi%i_start:dimsi%i_end,                                       &
              dimsi%j_start:dimsi%j_end, bl_levels),                           &
                 ! Exchange coefficients for
                 ! momentum, on UV-grid with
                 ! first and last j_end ignored.
                 ! for K>=2, between rho level K and K-1.
                 ! i.e. assigned at theta level K-1
   rhogamuv_uv(dimsi%i_start:dimsi%i_end,                                      &
               dimsi%j_start:dimsi%j_end, 2:bl_levels),                        &
                 ! Counter Gradient Term for U or V
                 ! defined on UV-grid
   u_v(dimsi_s%i_start:dimsi_s%i_end,                                          &
       dimsi_s%j_start:dimsi_s%j_end,bl_levels),                               &
                 ! Westerly_Southerly component of wind.
   tau_xy_fd_uv(dimsi%i_start:dimsi%i_end,                                     &
                dimsi%j_start:dimsi%j_end, bl_levels)
                 ! X/Y-component of form-drag stress
                 !    at a UV point

! Intent INOUT Variables
REAL(KIND=real_umphys), INTENT(IN OUT) ::                                      &
   tau_x_y (dimso%i_start:dimso%i_end,                                         &
            dimso%j_start:dimso%j_end, bl_levels)
                 ! explicit x_y-component of
                 ! turbulent stress at levels
                 ! k-1/2; eg. TAUX(,1) is surface
                 ! stress. UV-grid, 1st and last j_end
                 ! set to "missing data". (N/sq m)

! Intent OUT Variables
REAL(KIND=real_umphys), INTENT(OUT) ::                                         &
   tau_grad(dimso%i_start:dimso%i_end,                                         &
            dimso%j_start:dimso%j_end,bl_levels),                              &
                 ! k*du/dz grad stress (kg/m/s2)
   tau_count_grad(dimso%i_start:dimso%i_end,                                   &
                  dimso%j_start:dimso%j_end,bl_levels)
                 ! Counter gradient stress (kg/m/s2)

! LOCAL VARIABLES.
INTEGER ::                                                                     &
   i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MYM_EX_FLUX_UV'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

k=1
DO j = dimso%j_start, dimso%j_end
  DO i = dimso%i_start, dimso%i_end
    tau_grad(i,j,k) = 0.0
    tau_count_grad(i,j,k) = 0.0
  END DO
END DO

DO k = 2, bl_levels
  DO j = dimso%j_start, dimso%j_end
    DO i = dimso%i_start, dimso%i_end

      tau_grad(i,j,k) = rhokm_u_v(i,j,k) *                                     &
                     ( u_v(i,j,k) - u_v(i,j,k-1) ) *rdz_u_v(i,j,k)
      tau_count_grad(i,j,k) = rhogamuv_uv(i, j, k)
      tau_x_y(i,j,k) = tau_grad(i,j,k) + tau_count_grad(i,j,k)

    END DO
  END DO
END DO

! Add explicit orographic stress, noting that the surface stress
! is to be added later

IF (formdrag  ==  explicit_stress) THEN
  DO k = 2, bl_levels
    DO j = dimso%j_start, dimso%j_end
      DO i = dimso%i_start, dimso%i_end
        tau_x_y(i,j,k) = tau_x_y(i,j,k) + tau_xy_fd_uv(i,j,k)
      END DO
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mym_ex_flux_uv
END MODULE mym_ex_flux_uv_mod
