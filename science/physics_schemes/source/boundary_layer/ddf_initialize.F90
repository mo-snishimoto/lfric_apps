! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To set the initial TKE in the first order closure model

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
MODULE ddf_initialize_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DDF_INITIALIZE_MOD'
CONTAINS

SUBROUTINE ddf_initialize(                                                     &
      bl_levels,                                                               &
      z_uv, z_tq, dbdz, dvdzm, r_mosurf, fb_surf, u_s, h_pbl,                  &
      e_trb)

USE atm_fields_bounds_mod, ONLY: tdims, pdims, tdims_s
USE mym_const_mod, ONLY: e_trb_max
USE mym_option_mod, ONLY: tke_levels, l_my_extra_level,                        &
                          my_z_extra_fact, my_lowest_pd_surf,                  &
                          tke_cm_mx, tke_cm_fa
USE parkind1, ONLY: jprb, jpim
USE planet_constants_mod, ONLY: vkman
USE yomhook, ONLY: lhook, dr_hook
USE ddf_mix_length_mod, ONLY: ddf_mix_length
USE mym_calcphi_mod, ONLY: mym_calcphi
USE mym_diff_matcoef_mod, ONLY: mym_diff_matcoef
USE mym_implic_mod, ONLY: mym_implic
IMPLICIT NONE

! Intent In Variables
INTEGER, INTENT(IN) ::                                                         &
   bl_levels
                  ! Max. no. of "boundary" levels

REAL(KIND=real_umphys), INTENT(IN) ::                                          &
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels+1),                                                          &
                  ! Z_UV(*,K) is height of u level k
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                  ! IN Z_TQ(*,K) is height of theta level k.
                  ! Cloud ice (kg per kg air)
   dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                  ! Buoyancy gradient across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                  ! Modulus of wind shear at theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                  ! reciprocal of Monin-Obkhov length
   fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
                  ! Surface flux buoyancy over density (m^2/s^3)
   u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                  ! Surface friction velocity
   h_pbl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                  ! height of PBL determined by vertical profile
                  ! of SL

REAL(KIND=real_umphys), INTENT(OUT) ::                                         &
   e_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
                                                  bl_levels)
                  ! TKE defined on theta levels K-1

! Local variables
INTEGER ::                                                                     &
   i, j, k, ll,                                                                &
   itr_ini

REAL(KIND=real_umphys) ::                                                      &
   r_pr,                                                                       &
   elq,                                                                        &
   sm,                                                                         &
   sh,                                                                         &
   gm,                                                                         &
   gh

REAL(KIND=real_umphys) ::                                                      &
   ekw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
   elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
   coef_cm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
   coef_ce(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
   pdk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
   dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,            &
       tke_levels),                                                            &
   aa(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
   bb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
   cc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
   pdk0(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
   pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
   phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

REAL(KIND=real_umphys), PARAMETER ::                                           &
   pr = 0.7,                                                                   &
                  ! Prandtl number
                  ! only in the initialization,
                  ! constant prandtl number is assumed.
   diff_fact = 2.0
                  ! factor of a diffusion coef of E_TRB to that of
                  ! momentum

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DDF_INITIALIZE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

r_pr = 1.0 / pr

IF (my_lowest_pd_surf == 0) THEN
  l_my_extra_level = .FALSE.
  my_z_extra_fact = 1.0
END IF

! initial guess for e_trb, assuming neutral layer
! and set some parameters
DO k = 2, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF (z_tq(i, j, k - 1) < h_pbl(i, j)) THEN
        coef_cm(i, j, k) = tke_cm_mx
      ELSE
        coef_cm(i, j, k) = tke_cm_fa
      END IF
      sm = coef_cm(i, j, k)
      sh = coef_cm(i, j, k) * r_pr
      gm = dvdzm(i, j, k) ** 2
      gh = -dbdz(i, j, k)
      pdk(i, j, k) = sm * gm + sh * gh
      IF (pdk(i, j, k) <= 0.0) THEN
        pdk(i, j, k) = 0.0
        e_trb(i, j, k) = 0.0
      ELSE
        e_trb(i, j, k) = 1.0e-5
      END IF
    END DO
  END DO
END DO

IF (my_lowest_pd_surf > 0) THEN
  CALL mym_calcphi(                                                            &
        bl_levels, z_tq, r_mosurf, pmz, phh)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      pdk0(i, j) = 1.0 * u_s(i, j) ** 3 * pmz(i, j)                            &
         / (vkman * z_tq(i, j, 1))
    END DO
  END DO
END IF  ! IF MY_lowest_pd_surf

itr_ini = tke_levels + 1

DO ll = 1, itr_ini
  CALL ddf_mix_length(                                                         &
    tdims%i_end, tdims%j_end, 0, 0, bl_levels,                                 &
    z_uv, z_tq, dbdz, r_mosurf, fb_surf, h_pbl, e_trb,                         &
    elm, coef_ce, ekw)

  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (e_trb(i, j, k) <= 0.0) THEN
          ekw(i, j, k) = 0.0
        END IF
        dfm(i, j, k) = coef_cm(i, j, k) * ekw(i, j, k) * elm(i, j, k)
      END DO
    END DO
  END DO

  CALL mym_diff_matcoef(                                                       &
        bl_levels, diff_fact, dfm, aa, bb, cc)

  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (bb(i, j, k) == 0.0) THEN
          aa(i, j, k) = 0.0
          bb(i, j, k) = 1.0
          cc(i, j, k) = 0.0
          e_trb(i, j, k) = 0.0
        ELSE
          elq = ekw(i, j, k) * elm(i, j, k)
          aa(i, j, k) = - aa(i, j, k)
          bb(i, j, k) = - bb(i, j, k)                                          &
                  + ekw(i, j, k) * coef_ce(i, j, k)                            &
                                     / MAX(elm(i, j, k), 1.0e-20)
          bb(i, j, k) = SIGN(MAX(ABS(bb(i, j, k)), 1.0e-20),                   &
                                    bb(i, j, k))

          cc(i, j, k) = - cc(i, j, k)
          e_trb(i, j, k) = elq * pdk(i, j, k)
        END IF
      END DO
    END DO
  END DO

  IF (my_lowest_pd_surf > 0) THEN
    k = 2
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (bb(i, j, k) /= 0.0 .AND. pdk(i, j, k) > 0.0) THEN
          e_trb(i, j, k) = pdk0(i, j)
        END IF
      END DO
    END DO
  END IF   ! IF MY_lowest_pd_surf > 0

  CALL mym_implic(                                                             &
                  tke_levels, 2, tke_levels, aa, bb, cc, e_trb)
END DO  ! DO ll = 1, itr_ini

DO k = 2, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      e_trb(i, j, k) = MIN(                                                    &
                            MAX(e_trb(i, j, k), 1.0e-20),                      &
                                e_trb_max)
    END DO
  END DO
END DO

DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    e_trb(i, j, 1) = 0.0
  END DO
END DO

DO k = tke_levels + 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      e_trb(i, j, k) = 0.0
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ddf_initialize
END MODULE ddf_initialize_mod
