! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To calculate tri-diagonal matrix elements due to
!           diffusion for the prognostic variables in the MY model

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
MODULE mym_diff_matcoef_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MYM_DIFF_MATCOEF_MOD'
CONTAINS

SUBROUTINE mym_diff_matcoef(bl_levels, coef, z_uv, z_tq, dfm, aa, bb, cc)

USE atm_fields_bounds_mod, ONLY: pdims, tdims_s, tdims
USE mym_option_mod, ONLY:                                                      &
      l_my_extra_level, my_z_extra_fact, tke_levels
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Intent IN Variables
INTEGER, INTENT(IN) ::                                                         &
   bl_levels
                 ! Max. no. of "boundary" level

REAL(KIND=real_umphys), INTENT(IN) ::                                          &
   coef
                 ! factor for the diffusion coefficients to those for
                 ! momentum

REAL(KIND=real_umphys), INTENT(IN) ::                                          &
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
       bl_levels+1),                                                           &
                 ! Z_UV(*,K) is height of u level k
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
       bl_levels),                                                             &
                 ! Z_TQ(*,K) is height of theta level k
   dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,            &
       bl_levels)
                 ! diffusion coefficients for momentum

! Intent OUT variables
REAL(KIND=real_umphys), INTENT(OUT) ::                                         &
                ! coefficients of tri-diagonal equations
                ! due to diffusion
   aa(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels),         &
                ! coefs of fields on level K-1
   bb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels),         &
                ! coefs of fields on level K
   cc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels)
                ! coefs of fields on level K+1

! Local variables
INTEGER ::                                                                     &
   i, j, k, k_start
                ! Loop indexes
REAL(KIND=real_umphys) ::                                                      &
   km_m1,                                                                      &
                ! diffusion coefficient on lower level by one
   km_p1
                ! diffusion coefficient on upper level by one

REAL(KIND=real_umphys) ::                                                      &
   r_dr_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,               &
            tke_levels),                                                       &
                ! reciprocal of grid spaces of rho levels
   r_dr_theta(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
              tke_levels),                                                     &
                ! reciprocal of grid spaces of theta levels
   weight1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
   weight2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels)
                ! weight to interporate variables on theta levels
                ! onto rho levels

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MYM_DIFF_MATCOEF'

! Calculate and save r_dr and weight

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

k = 1
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    r_dr_theta(i, j, k) = 1.0 / (z_uv(i, j, k + 1) - z_uv(i, j, k))
    r_dr_rho(i, j, k) = 1.0 / z_tq(i, j, k)
    weight1(i, j, k) = z_uv(i, j, k) * r_dr_rho(i, j, k)
    weight2(i, j, k) = (z_tq(i, j, k) - z_uv(i, j, k)) * r_dr_rho(i, j, k)
  end do
end do

DO k = 2, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      r_dr_theta(i, j, k) = 1.0 / (z_uv(i, j, k + 1) - z_uv(i, j, k))
      r_dr_rho(i, j, k) = 1.0 / (z_tq(i, j, k) - z_tq(i, j, k - 1))

      weight1(i, j, k) = (z_uv(i, j, k) - z_tq(i, j, k - 1)) * r_dr_rho(i, j, k)
      weight2(i, j, k) = (z_tq(i, j, k) - z_uv(i, j, k)) * r_dr_rho(i, j, k)
    END DO
  END DO
END DO

! Calculate aa, bb, cc
k = 2
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    km_m1 = coef * dfm(i, j, k)
    km_p1 = coef *                                                             &
                (weight1(i, j, k) * dfm(i, j, k + 1)                           &
              +  weight2(i, j, k) * dfm(i, j, k))

    cc(i, j, k) = km_p1 * r_dr_rho(i, j, k)                                    &
                                    * r_dr_theta(i, j, k - 1)
    aa(i, j, k) = km_m1 * r_dr_rho(i, j, k - 1)                                &
                                    * r_dr_theta(i, j, k - 1)
    bb(i, j, k) = -aa(i, j, k) - cc(i, j, k)
  END DO
END DO

DO k = 3, tke_levels - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      km_m1 = coef *                                                           &
                (weight1(i, j, k - 1) * dfm(i, j, k)                           &
              +  weight2(i, j, k - 1) * dfm(i, j, k - 1))
      km_p1 = coef *                                                           &
                (weight1(i, j, k) * dfm(i, j, k + 1)                           &
              +  weight2(i, j, k) * dfm(i, j, k))

      cc(i, j, k) = km_p1 * r_dr_rho(i, j, k)                                  &
                                       * r_dr_theta(i, j, k - 1)
      aa(i, j, k) = km_m1 * r_dr_rho(i, j, k - 1)                              &
                                       * r_dr_theta(i, j, k - 1)
      bb(i, j, k) = -aa(i, j, k) - cc(i, j, k)

    END DO
  END DO
END DO

k = tke_levels
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    km_m1 = coef *                                                             &
            (weight1(i, j, k - 1) * dfm(i, j, k)                               &
          +  weight2(i, j, k - 1) * dfm(i, j, k - 1))

    km_p1 = coef * weight2(i, j, k) * dfm(i, j, k)

    cc(i, j, k) = km_p1 * r_dr_rho(i, j, k)                                    &
                                       * r_dr_theta(i, j, k - 1)
    aa(i, j, k) = km_m1 * r_dr_rho(i, j, k - 1)                                &
                                       * r_dr_theta(i, j, k - 1)
    bb(i, j, k) = -aa(i, j, k) - cc(i, j, k)

  END DO
END DO

IF (l_my_extra_level) THEN
  k_start = 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      aa(i, j, 1) = 0.0
      cc(i, j, 1) = coef * dfm(i, j, 2)                                        &
              / (z_tq(i, j, 1)                                                &
                  * my_z_extra_fact) ** 2

      bb(i, j, 1) = - aa(i, j, 1) - cc(i, j, 1)
    END DO
  END DO
ELSE
  k_start = 2
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      aa(i, j, 1) = 0.0
      bb(i, j, 1) = 0.0
      cc(i, j, 1) = 0.0
    END DO
  END DO
END IF

DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    aa(i, j, k_start) = 0.0
    cc(i, j, tke_levels) = 0.0
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mym_diff_matcoef
END MODULE mym_diff_matcoef_mod
