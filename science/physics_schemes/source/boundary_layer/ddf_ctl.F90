! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: The main subroutine for the first order closure model
!           based on Deardorff (1980).

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
MODULE ddf_ctl_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DDF_CTL_MOD'
CONTAINS

SUBROUTINE ddf_ctl(                                                            &
! IN levels/switches
      bl_levels, BL_diag,                                                      &
! IN fields
      z_uv,z_tq, u_p, v_p, qw, tl, t, q, qcl, qcf,                             &
      p_theta_levels, p_half, bq_gb, bt_gb, rho_mix, rho_wet_tq,               &
      dtldzm, dqwdzm, dudz, dvdz, dbdz, dvdzm, delta_smag, u_s, fb_surf, pstar,&
! INOUT fields
      e_trb, rhokm, rhokh, zhpar_shcu)

USE atm_fields_bounds_mod, ONLY: tdims_l, tdims, pdims, tdims_s
USE bl_diags_mod, ONLY: strnewbldiag
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE model_domain_mod, ONLY: model_type, mt_single_column
USE mym_const_mod, ONLY: e_trb_max
USE mym_option_mod, ONLY: my_ini_dbdz_min, tke_cm_mx, l_shcu_buoy,             &
      l_my_condense, tke_cm_fa, my_lowest_pd_surf, tke_levels,                 &
      l_my_ini_zero, l_my_initialize

USE parkind1, ONLY: jprb, jpim
USE planet_constants_mod, ONLY: vkman, kappa, pref, c_virtual, grcp, g
USE yomhook, ONLY: lhook, dr_hook

USE ddf_initialize_mod, ONLY: ddf_initialize
USE ddf_mix_length_mod, ONLY: ddf_mix_length
USE mym_calcphi_mod, ONLY: mym_calcphi
USE mym_condensation_mod, ONLY: mym_condensation
USE mym_const_set_mod, ONLY: mym_const_set
USE mym_shcu_buoy_mod, ONLY: mym_shcu_buoy
USE mym_update_fields_mod, ONLY: mym_update_fields
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
   u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),         &
                  ! U on P-grid.
   v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),         &
                  ! V on P-grid.
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),         &
                  ! Total water content
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),         &
                  ! Ice/liquid water temperature
   t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),          &
                  ! Temperature
   q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                  ! specific humidity
   qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
       tdims_l%k_start:bl_levels),                                             &
                  ! Cloud liquid water
   qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
       tdims_l%k_start:bl_levels),                                             &
                  ! Cloud ice (kg per kg air)
   p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  0:bl_levels+1),                                              &
                  ! Pressure at theta level
   p_half(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                 &
          bl_levels),                                                          &
                  ! Pressure on rho levels (Pa)
   bq_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                  ! A grid-box mean buoyancy param
                  ! on T,q-levels (full levels).
   bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                  ! A grid-box mean buoyancy param
                  ! on T,q-levels (full levels).
   rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           bl_levels+1),                                                       &
                  ! density on UV (ie. rho) levels;
                  ! used in RHOKH so dry density if
                  ! L_mr_physics is true
   rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
              bl_levels),                                                      &
                  ! density on TQ (ie. theta) levels;
                  ! used in RHOKM so wet density
   dtldzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
           2:bl_levels),                                                       &
                  ! gradient of TL across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   dqwdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          2:bl_levels),                                                        &
                  ! gradient of QW across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   dudz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:bl_levels),                                                          &
                  ! Gradient of u at theta levels.
                  !(:,:,K) repserents the value on theta level K-1
   dvdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:bl_levels),                                                          &
                  ! Gradient of v at theta levels.
                  !(:,:,K) repserents the value on theta level K-1
   dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                  ! Buoyancy gradient across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                  ! Modulus of wind shear at theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   delta_smag(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                  ! IN delta_x used by Smagorinsky
   u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                  ! Surface friction velocity
   fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
                  ! Surface flux buoyancy over density (m^2/s^3)
   pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                  ! surface pressure

REAL(KIND=real_umphys), INTENT(IN OUT) ::                                      &
   e_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
                                                  bl_levels),                  &
                  ! TKE defined on theta levels K-1
   rhokm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,          &
                                                  bl_levels),                  &
                  ! Exchange coeffs for momentum
                  ! between K and K-1 on rho levels.
                  ! i.e. the coeffs are defined on theta level K-1.
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels),                                                           &
                  ! Exchange coeffs for scalars
                  ! between K and K-1 on theta levels.
                  ! i.e. the coeffs are defined on rho levels
   zhpar_shcu(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                  ! Height of mixed layer used to evaluate
                  ! the non-gradient buoyancy flux

!  Declaration of BL diagnostics.
TYPE (strnewbldiag), INTENT(IN OUT) :: BL_diag

! Local Variables
INTEGER ::                                                                     &
   i, j, k
                  ! Loop indexes

REAL(KIND=real_umphys) ::                                                      &
   r_weight1,                                                                  &
                  ! weight factor to interpolate variables on rho
                  ! levels onto theta levels
   weight2,                                                                    &
                  ! weight factor to interpolate variables on rho
                  ! levels onto theta levels
   weight3,                                                                    &
                  ! weight factor to interpolate variables on rho
                  ! levels onto theta levels
   taux,                                                                       &
                  ! stress of x-direction
   tauy,                                                                       &
                  ! stress of y-direction
   r_pr,                                                                       &
                  ! reciprocal of the Prandtl number
   coef_cm
                  ! coefficient appeared in determining a diffusion
                  ! coefficients

INTEGER ::                                                                     &
   flag_calc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                  ! flag to indicate whether the column should be
                  ! calculated

REAL(KIND=real_umphys) ::                                                      &
   r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                  ! reciprocal of Monin-Obkhov length
   rhokh_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                  ! density on TQ (ie. theta) levels;
                  ! used in RHOKM so wet density
   prod(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                  ! production term
   disp_coef(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             tke_levels),                                                      &
                  ! coefficients of E_TRB in a dissipation term
   pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                  ! gradient function for momentum at the surface
                  ! minus non-dimensional height (height / MO length)
   phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                  ! gradient function for scalars at the surface
   elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! mixing length
   ekw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! sqrt(e_trb) on theta level K-1
   coef_ce(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                  ! coefficient appeared in a dissipation term
   sl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! static energy
   h_pbl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
                  ! height of PBL determined by vertical profile
                  ! of SL
   tsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
       bl_levels),                                                             &
                  ! Self covariance of liquid potential temperature
                  ! (thetal'**2) defined on theta levels K-1
   qsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
       bl_levels),                                                             &
                  ! Self covariance of total water
                  ! (qw'**2) defined on theta levels K-1
   cov(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
       bl_levels),                                                             &
                  ! Correlation between thetal and qw
                  ! (thetal'qw') defined on theta levels K-1
   vt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! Buoyancy parameter for FTL (excluding g/thetav)
                  ! on theta level K-1
   vq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! Buoyancy parameter for FQW (excluding g/thetav)
                  ! on theta level K-1
   tv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! Virtual temperature on theta level K-1
   dbdz_l(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
                  ! Buoyancy gradient across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   exner(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
                  ! exner function on theta level K-1
   gtr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! G/thetav on theta level K-1
   q1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! normalized excessive water from the saturation
                  ! on theta level K-1
   cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! cloud fraction derived by the bi-normal
                  ! distribution on theta level K-1
   ql(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! condensed liquid water derived by the bi-normal
                  ! distribution on theta level K-1
   prod_m(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
                  ! production term by wind shear
   prod_h(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
                  ! production term by buoyancy
   wb_ng(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
                  ! buoyancy flux related to the skewness
                  ! on theta level K-1
   frac_shcu(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             tke_levels)
                  ! cloud fraction corrected by shallow cumulus
                  ! process on theta level K-1

LOGICAL, SAVE ::                                                               &
   l_first = .TRUE.
                  ! flag to indicate if it is the first execution

INTEGER, PARAMETER ::                                                          &
  levflag = 2
                  ! For using subroutines for the MY model.

REAL(KIND=real_umphys), PARAMETER ::                                           &
  c_corr = 2.0
                  ! coefficient appeared in parameterizing the width
                  ! of the bi-normal distribution function

REAL(KIND=real_umphys), PARAMETER ::                                           &
   diff_fact = 2.0
                  ! factor of a diffusion coef of E_TRB to that of
                  ! momentum

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'DDF_CTL'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Calculate Monin-Obukov Length
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    r_mosurf(i,j)= -vkman*fb_surf(i,j)                                         &
                     / MAX(u_s(i,j)*u_s(i,j)*u_s(i,j), TINY(1.0))
  END DO
END DO

! Calculate gradient functions
IF (my_lowest_pd_surf > 0) THEN
  CALL mym_calcphi(                                                            &
        bl_levels, z_uv, r_mosurf, pmz, phh)
END IF

! Calculate static energy to determine the top of mixed layer
DO k = 1, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      sl(i, j, k) = tl(i, j, k) + grcp * z_tq(i, j, k)
      sl(i, j, k) = sl(i, j, k) * (1.0 + c_virtual * q(i, j, k)                &
                             - qcl(i, j, k) - qcf(i, j, k))
    END DO
  END DO
END DO

! Determine the height of the top of mixed layer
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    flag_calc(i, j) = 1
    h_pbl(i, j) = z_tq(i, j, 1)
  END DO
END DO
DO k = 2, tke_levels - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF (flag_calc(i, j) == 1) THEN
        IF (sl(i, j, k) > sl(i, j, 1)) THEN
          h_pbl(i, j) = z_tq(i, j, k - 1)                                      &
               + (z_tq(i, j, k) - z_tq(i, j, k - 1))                           &
                 * (sl(i, j, 1) - sl(i, j, k - 1))                             &
               / (sl(i, j, k) - sl(i, j, k - 1))
          flag_calc(i, j) = 0
        END IF
      END IF
    END DO
  END DO
END DO

! Initialization. Executed only once.
! In the initialization, balance between production and dissipation
! is assumed. Diffusion coeffients required to determine production
! terms are calculated with stability functions.
IF (l_first) THEN
  CALL mym_const_set

  ! IF the first value of e_trb has been set to be missing by the
  ! reconfiguration, the initialization for the whole domain
  ! is essential.
  IF (l_my_initialize) THEN
    IF (l_my_ini_zero) THEN
      DO k = 1, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            e_trb(i, j, k) = 0.0
          END DO
        END DO
      END DO
    ELSE  ! not l_my_ini_zero
      ! Initialize the prognostic variables by assuming the balance
      ! between production and dissipation terms

      ! In the initialization, DBDZ by the LS cloud scheme is used.
      ! To avoid to diagnose huge TKE, the lower limit for DBDZ
      ! is imposed.
      DO k = 2, tke_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            dbdz_l(i, j, k) = MAX(my_ini_dbdz_min, dbdz(i, j, k))
          END DO
        END DO
      END DO
      CALL ddf_initialize(                                                     &
         bl_levels,                                                            &
         z_uv, z_tq, dbdz_l, dvdzm, delta_smag, r_mosurf, fb_surf, u_s, h_pbl, &
         e_trb)
      ! Above tke_levels, the prognostic variables should be zeros.
      DO k = tke_levels + 1, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            e_trb(i, j, k) = 0.0
          END DO
        END DO
      END DO
    END IF  ! test if l_my_ini_zero
  END IF  ! test if l_my_initialize

  IF (l_shcu_buoy .and. l_my_initialize) THEN
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        zhpar_shcu(i, j) = z_tq(i, j, tke_levels-1)
      END DO
    END DO
  END IF
  l_first = .FALSE.
END IF

CALL ddf_mix_length(                                                           &
    tdims%i_end,tdims%j_end,tdims_l%halo_i,tdims_l%halo_j, bl_levels,          &
    z_uv, z_tq, dbdz, delta_smag, r_mosurf, fb_surf, h_pbl, e_trb,             &
    elm, coef_ce, ekw)

  ! Calculate diffusion coefficients
DO k = 2, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF (z_tq(i, j, k - 1) < h_pbl(i, j)) THEN
        coef_cm = tke_cm_mx
      ELSE
        coef_cm = tke_cm_fa
      END IF

      r_pr = 1.0 + 2.0 * elm(i, j, k)                                          &
              / (z_uv(i, j, k) - z_uv(i, j, k - 1))
      rhokm(i, j, k) = coef_cm * elm(i, j, k) * ekw(i, j, k)
      rhokh_tq(i, j, k) = rhokm(i, j, k) * r_pr
    END DO
  END DO
END DO

! Set virtual temperature, exner function and g/thetav
DO k = 2, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      tv(i, j, k) = t(i, j, k - 1)                                             &
                   * (1.0 + c_virtual * q(i, j, k - 1)                         &
                         - qcl(i, j, k - 1) - qcf(i, j, k - 1))
      exner(i, j, k) =                                                         &
                (p_theta_levels(i, j, k - 1) / pref) ** kappa
      gtr(i, j, k) = g / tv(i, j, k) * exner(i, j, k)
    END DO
  END DO
END DO

! The covariances to be required by mym_condensation
! are diagnosed assuming balance between
! production and dissipation.
IF (l_my_condense .OR. l_shcu_buoy) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        tsq(i, j, k) = c_corr * elm(i, j, k) ** 2                              &
                               * dtldzm(i, j, k) ** 2
        qsq(i, j, k) = c_corr * elm(i, j, k) ** 2                              &
                               * dqwdzm(i, j, k) ** 2
        cov(i, j, k) = c_corr * elm(i, j, k) ** 2                              &
                               * dtldzm(i, j, k) * dqwdzm(i, j, k)
      END DO
    END DO
  END DO

  CALL mym_condensation(                                                       &
  ! IN levels/switches
            bl_levels, levflag,                                                &
            BL_diag,                                                           &
  ! IN fields
            qw, tl, t, p_theta_levels, tsq, qsq, cov,                          &
  ! OUT fields
            vt, vq, q1, cld, ql)
END IF

IF (.NOT. l_my_condense) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ! convert buoy params from the UM notation to the MY notaation
        vt(i, j, k) = bt_gb(i, j, k - 1) * tv(i, j, k)
        vq(i, j, k) = bq_gb(i, j, k - 1) * tv(i, j, k)                         &
                                             / exner(i, j, k)
      END DO
    END DO
  END DO
END IF

IF (l_shcu_buoy) THEN
  CALL mym_shcu_buoy(                                                          &
  ! IN levels/switches
             bl_levels, BL_diag,                                               &
  ! IN fields
             fb_surf, u_s, pstar, z_tq, z_uv, p_theta_levels, p_half,          &
             u_p, v_p, t, q, qcl, qcf, q1, cld,                                &
  ! INOUT / OUT fields
             zhpar_shcu, frac_shcu, wb_ng)
ELSE
  DO k = 1, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        wb_ng(i,j,k) = 0.0
        frac_shcu(i,j,k) = cld(i,j,k)
      END DO
    END DO
  END DO
END IF

! Calculate production terms and coefficient of dissipation term.
DO k = 2, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dbdz_l(i, j, k) = gtr(i, j, k)                                           &
                          * (vt(i, j, k) * dtldzm(i, j, k)                     &
                                 + vq(i, j, k) * dqwdzm(i, j, k))
      prod_h(i, j, k) = - rhokh_tq(i, j, k) * dbdz_l(i, j, k)                  &
                                 + wb_ng(i, j, k)

      taux = rhokm(i, j, k) * dudz(i, j, k)
      tauy = rhokm(i, j, k) * dvdz(i, j, k)

      prod_m(i, j, k) = taux * dudz(i, j, k) + tauy * dvdz(i, j, k)


      prod(i, j, k) = prod_m(i, j, k) + prod_h(i, j, k)
      disp_coef(i, j, k) = coef_ce(i, j, k) * ekw(i, j, k)                     &
                         / MAX(elm(i, j, k), 1.0e-20)

    END DO
  END DO
END DO

! Overwrite the production term at the lowest level by
! the one evaluated with surface fluxes.
IF (my_lowest_pd_surf > 0) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      prod(i, j, 2) = u_s(i, j) ** 3 * pmz(i, j)                               &
                             / (vkman * z_tq(i, j, 1))
    END DO
  END DO
END IF

CALL mym_update_fields(                                                        &
        bl_levels, diff_fact, z_uv, z_tq, rhokm, prod, disp_coef, e_trb)

DO k = tke_levels + 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      e_trb(i, j, k) = 0.0
      rhokm(i, j, k) = 0.0
      rhokh_tq(i, j, k) = 0.0
      rhokh(i, j, k) = 0.0
    END DO
  END DO
END DO

DO k = 2, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      e_trb(i, j, k) = MIN(MAX(e_trb(i, j, k), 1.0e-20), e_trb_max)
      rhokm(i, j, k) = rho_wet_tq(i, j, k - 1) * rhokm(i, j, k)
    END DO
  END DO
END DO

! Note "RHO" here is always wet density (RHO_WET_TQ) so
! save multiplication of RHOKH to after interpolation
IF (.NOT. l_mr_physics) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        rhokh_tq(i, j, k) = rho_wet_tq(i, j, k - 1) * rhokh_tq(i, j, k)
      END DO
    END DO
  END DO
END IF

! Interpolate RHOKH_TQ on theta levels to rho levels
DO k = 2, tke_levels - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      r_weight1 = 1.0 / (z_tq(i,j,k) -                                         &
                                   z_tq(i,j, k-1))
      weight2 = (z_tq(i,j,k) -                                                 &
                          z_uv(i,j,k)) * r_weight1
      weight3 = (z_uv(i,j,k) -                                                 &
                          z_tq(i,j,k-1)) * r_weight1
      rhokh(i,j,k) =                                                           &
                          weight3 * rhokh_tq(i,j,k+1)                          &
                         +weight2 * rhokh_tq(i,j,k)
    END DO
  END DO
END DO

k = tke_levels
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    r_weight1 = 1.0 / (z_tq(i,j,k) -                                           &
                                   z_tq(i,j, k-1))
    weight2 = (z_tq(i,j,k) -                                                   &
                          z_uv(i,j,k)) * r_weight1

    rhokh(i, j, k) = weight2 * rhokh_tq(i, j, k)
  END DO
END DO

! Finally multiply RHOKH by dry density
IF (l_mr_physics) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        rhokh(i, j, k) = rho_mix(i, j, k) * rhokh(i, j, k)
      END DO
    END DO
  END DO
END IF

IF (BL_diag%l_dbdz) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        BL_diag%dbdz(i, j, k) = dbdz_l(i, j, k)
      END DO
    END DO
  END DO
END IF

IF (BL_diag%l_dvdzm) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        BL_diag%dvdzm(i,j,k) = dvdzm(i, j, k)
      END DO
    END DO
  END DO
END IF

IF (BL_diag%l_tke_shr_prod) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        BL_diag%tke_shr_prod(i, j, k) = prod_m(i, j, k)
      END DO
    END DO
  END DO
END IF

IF (BL_diag%l_tke_boy_prod) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        BL_diag%tke_boy_prod(i, j, k) = prod_h(i, j, k)
      END DO
    END DO
  END DO
END IF

IF (BL_diag%l_tke_boy_prod) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        BL_diag%tke_dissp(i, j, k) =                                           &
                 coef_ce(i, j, k) * (ekw(i, j, k)) ** 3                        &
                           / MAX(elm(i, j, k), 1.0e-20)
      END DO
    END DO
  END DO
END IF

IF (BL_diag%l_elm) THEN
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        BL_diag%elm(i, j, k) = elm(i, j, k)
      END DO
    END DO
  END DO
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ddf_ctl
END MODULE ddf_ctl_mod
