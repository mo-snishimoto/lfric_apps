! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose: To calculate gradient functions at the surface used in
!           evaluating the production terms of the prognostic
!           variables at the lowest layer.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
MODULE mym_calcphi_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MYM_CALCPHI_MOD'
CONTAINS

SUBROUTINE mym_calcphi(bl_levels, z_tq, r_mosurf, pmz, phh)

USE atm_fields_bounds_mod, ONLY: tdims
USE mym_const_mod, ONLY: two_thirds, pr
USE mym_option_mod, ONLY:                                                      &
      businger, bh1991, my_lowest_pd_surf,                                     &
      l_my_extra_level, my_z_extra_fact
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Intent IN Variables
INTEGER, INTENT(IN) ::                                                         &
   bl_levels
                 ! number of boundary layer levels

REAL(KIND=real_umphys), INTENT(IN) ::                                          &
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                 ! Z_TQ(*,K) is height of theta
                 !    level k.
   r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                 ! reciprocal of Monin-Obkhov length

! Intent OUT Variables
REAL(KIND=real_umphys), INTENT(OUT) ::                                         &
   pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                 ! gradient function for momentum
                 ! at surface minus non-dimensional height
   phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                 ! gradient function for scalars
                 ! at surface

! Local variables
INTEGER ::                                                                     &
   i, j
                 ! Loop indexes

REAL(KIND=real_umphys) ::                                                      &
   zeta,                                                                       &
                 ! non-dimensional height
   tmp
                 ! work variable

REAL(KIND=real_umphys) ::                                                      &
   z_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                 ! height of the lowest layer

REAL(KIND=real_umphys), PARAMETER ::                                           &
                 ! coefficients appeared
                 !               in Beljaars and Holtslag(1991)
   bel_a = 1.0,                                                                &
   bel_b = 2.0 / 3.0,                                                          &
   bel_c = 5.0,                                                                &
   bel_d = 0.35

REAL(KIND=real_umphys), PARAMETER ::                                           &
   my_zeta_max = 2.0
                ! upper limit for zeta

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MYM_CALCPHI'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_my_extra_level) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      z_1(i, j) = z_tq(i, j, 1) * my_z_extra_fact
    END DO
  END DO
ELSE
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      z_1(i, j) = z_tq(i, j, 1)
    END DO
  END DO
END IF

IF (my_lowest_pd_surf == businger) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      zeta = MIN(z_1(i, j) * r_mosurf(i, j), my_zeta_max)
      IF (zeta >= 0.0) THEN
        pmz(i, j) = 1.0 + 4.7 * zeta
        phh(i, j) = pr + 4.7 * zeta
      ELSE
        pmz(i, j) = 1.0 / SQRT(SQRT(1.0 - 15.0 * zeta))
        phh(i, j) = pr / SQRT(1.0 - 9.0 * zeta)
      END IF
      pmz(i, j) = pmz(i, j) - zeta
    END DO
  END DO
ELSE IF (my_lowest_pd_surf == bh1991) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      zeta = MIN(z_1(i, j) * r_mosurf(i, j), my_zeta_max)
      IF (zeta >= 0) THEN
        tmp = bel_b * EXP(-bel_d * zeta)                                       &
             * (bel_d * zeta - bel_c - 1.0)
        pmz(i, j) = 1.0 - zeta * (tmp - bel_a)
        phh(i, j) = 1.0 - zeta * (tmp -                                        &
             SQRT(1.0 + two_thirds * bel_a * zeta))
      ELSE
        tmp = SQRT(1.0 - 16.0 * zeta)
        pmz(i, j) = 1.0 / SQRT(tmp)
        phh(i, j) = 1.0 / tmp
      END IF
      pmz(i, j) = pmz(i, j) - zeta
    END DO
  END DO
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE mym_calcphi
END MODULE mym_calcphi_mod
