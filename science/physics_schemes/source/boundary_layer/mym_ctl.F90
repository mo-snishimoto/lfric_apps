! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: The main subroutine for the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  This code is based on the code provided by the authors who wrote
!  the following papers.
!   * Nakanishi, M. and H. Niino, 2009: Development of an improved
!      turbulence closure model for the atmospheric boundary layer.
!      J. Meteor. Soc. Japan, 87, 895-912.
!   * Nakanishi, M. and H. Niino, 2006: An improved Mellor-Yamada
!      Level-3 model: Its numerical stability and application to
!      a regional prediction of advection fog.
!      Boundary-Layer Meteor., 119, 397-407.
!   * Nakanishi, M. and H. Niino, 2004: An improved Mellor-Yamada
!      Level-3 model with condensation physics: Its design and
!      verification.
!      Boundary-Layer Meteor., 112, 1-31.
!   * Nakanishi, M., 2001: Improvement of the Mellor-Yamada
!      turbulence closure model based on large-eddy simulation data.
!      Boundary-Layer Meteor., 99, 349-378.
!   The web site publicising their code:
!    http://www.nda.ac.jp/~naka/MYNN/index.html

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
MODULE mym_ctl_mod

USE um_types, ONLY: real_umphys, real_eps

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MYM_CTL_MOD'
CONTAINS

SUBROUTINE mym_ctl(                                                            &
! IN levels/switches
      bl_levels, levflag, nSCMDpkgs,L_SCMDiags,                                &
      BL_diag, cycleno,                                                        &
! IN fields
      z_uv,z_tq, u_p, v_p, qw, tl, t, q, qcl, qcf, bq_gb, bt_gb,               &
      rho_mix, rho_wet_tq, fqw, ftl,                                           &
      dtldzm, dqwdzm, dudz, dvdz, dbdz, dvdzm,                                 &
      p_theta_levels, p_half, u_s, fb_surf, pstar,                             &
! INOUT fields
      e_trb, tsq_trb, qsq_trb, cov_trb, rhokm, rhokh, zhpar_shcu,              &
! OUT fields
      rhogamu, rhogamv, rhogamt, rhogamq)

USE atm_fields_bounds_mod, ONLY: tdims, pdims, tdims_l, tdims_s
USE bl_diags_mod, ONLY: strnewbldiag
USE dynamics_input_mod, ONLY: numcycles
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE level_heights_mod, ONLY:                                                   &
  r_theta_levels, r_rho_levels
USE mym_option_mod, ONLY: l_my_condense, l_shcu_buoy,                          &
      my_lowest_pd_surf, tke_levels, l_my_initialize, l_my_ini_zero,           &
      my_ini_dbdz_min, l_3dtke
USE missing_data_mod, ONLY: rmdi
USE parkind1, ONLY: jprb, jpim
USE planet_constants_mod, ONLY: vkman, kappa, pref, c_virtual, g
USE turb_diff_ctl_mod, ONLY: visc_m, visc_h
USE yomhook, ONLY: lhook, dr_hook
USE mym_calcphi_mod, ONLY: mym_calcphi
USE mym_condensation_mod, ONLY: mym_condensation
USE mym_const_set_mod, ONLY: mym_const_set
USE mym_initialize_mod, ONLY: mym_initialize
USE mym_shcu_buoy_mod, ONLY: mym_shcu_buoy
USE mym_turbulence_mod, ONLY: mym_turbulence
IMPLICIT NONE

! Intent In Variables

INTEGER, INTENT(IN) ::                                                         &
   bl_levels,                                                                  &
                  ! Max. no. of "boundary" levels
   levflag,                                                                    &
                  ! to indicate the level of the MY model
                  ! 2: level 2.5
                  ! 3: level 3
   cycleno        ! Iteration number (EG outer loop)

REAL(KIND=real_umphys), INTENT(IN) ::                                          &
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels+1),                                                          &
                  ! Z_UV(*,K) is height of u level k
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                  ! Z_TQ(*,K) is height of theta level k.
   u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),         &
                  ! U on P-grid.
   v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),         &
                  ! V on P-grid.
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),         &
                  ! Total water content
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),         &
                  ! Ice/liquid water temperature
   t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                  ! temperature
   q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                  ! specific humidity
   qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
       tdims_l%k_start:bl_levels),                                             &
                  ! Cloud liquid water
   qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
       tdims_l%k_start:bl_levels),                                             &
                  ! Cloud ice (kg per kg air)
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
   fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                  ! Moisture flux between layers
                  ! (kg per square metre per sec).
                  ! FQW(,1) is total water flux
                  ! from surface, 'E'.
                  ! defined on rho levels
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                  ! FTL(,K) contains net turbulent
                  ! sensible heat flux into layer K
                  ! from below; so FTL(,1) is the
                  ! surface sensible heat, H. (W/m2)
                  ! defined on rho levels
   dtldzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          2:bl_levels),                                                        &
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
        2:bl_levels),                                                          &
                  ! Buoyancy gradient across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                  ! Modulus of wind shear at theta levels.
                  ! (:,:,K) represents the value on theta level K-1
   p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  0:bl_levels+1),                                              &
                  ! Pressure on theta levels (Pa)
   p_half(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                 &
          bl_levels),                                                          &
                  ! Pressure on rho levels (Pa)
   u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                  ! Surface friction velocity
   fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
                  ! Surface flux buoyancy over
                  ! density (m^2/s^3)
   pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                  ! surface pressure

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                         &
   nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                         &
   L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

! Intent INOUT variables
REAL(KIND=real_umphys), INTENT(IN OUT) ::                                      &
   e_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
          bl_levels),                                                          &
                  ! TKE defined on theta levels K-1
   tsq_trb(tdims%i_start:tdims%i_end,                                          &
           tdims%j_start:tdims%j_end,bl_levels),                               &
                  ! Self covariance of liquid potential temperature
                  ! (thetal'**2) defined on theta levels K-1
   qsq_trb(tdims%i_start:tdims%i_end,                                          &
           tdims%j_start:tdims%j_end,bl_levels),                               &
                  ! Self covariance of total water
                  ! (qw'**2) defined on theta levels K-1
   cov_trb(tdims%i_start:tdims%i_end,                                          &
           tdims%j_start:tdims%j_end,bl_levels),                               &
                  ! Correlation between thetal and qw
                  ! (thetal'qw') defined on theta levels K-1
   rhokm(tdims_s%i_start:tdims_s%i_end,                                        &
         tdims_s%j_start:tdims_s%j_end,bl_levels),                             &
                  ! Exchange coeffs for momentum
                  ! between K and K-1 on rho levels.
                  ! i.e. the coeffs are defined on theta level K-1.
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels),                                                           &
                  ! Exchange coeffs for scalars
                  ! between K and K-1 on theta levels.
                  ! i.e. the coeffs are defined on rho levels
   zhpar_shcu(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                  ! Height of mixed layer used to evaluate
                  ! the non-gradient buoyancy flux

!  Declaration of BL diagnostics.
TYPE (strnewbldiag), INTENT(IN OUT) :: BL_diag

! Intent Out Variables
REAL(KIND=real_umphys), INTENT(OUT) ::                                         &
   rhogamu(tdims_s%i_start:tdims_s%i_end,                                      &
           tdims_s%j_start:tdims_s%j_end,2:bl_levels),                         &
                  ! Counter gradient terms for TAUX
                  ! defined at theta level K-1
   rhogamv(tdims_s%i_start:tdims_s%i_end,                                      &
           tdims_s%j_start:tdims_s%j_end,2:bl_levels),                         &
                  ! Counter gradient terms for TAUY
                  ! defined at theta level K-1
   rhogamt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels),                                                       &
                  ! Counter gradient terms for FTL
                  ! defined at rho levels
   rhogamq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels)
                  ! Counter gradient terms for FQW
                  ! defined at rho levels

! Local Variables
INTEGER ::                                                                     &
   i, j, k
                  ! loop indexes

LOGICAL, SAVE ::                                                               &
   l_first = .TRUE.
                  ! flag to indicate if it is the first execution
REAL(KIND=real_umphys) ::                                                      &
   r_weight1,                                                                  &
                  ! weight factor to interpolate variables on rho
                  ! levels onto theta levels
   weight2,                                                                    &
                  ! weight factor to interpolate variables on rho
                  ! levels onto theta levels
   weight3
                  ! weight factor to interpolate variables on rho
                  ! levels onto theta levels
REAL(KIND=real_umphys) ::                                                      &
   pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                  ! gradient function for momentum at the surface
                  ! minus non-dimensional height (height / MO length)
   phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                  ! gradient function for scalars at the surface
   r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                  ! reciprocal of Monin-Obukhov length
   qke(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
                            bl_levels),                                        &
                  ! twice of TKE (denoted to q**2) on theta level K-1
   dbdz_l(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          2:tke_levels),                                                       &
                  ! Buoyancy gradient across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   vt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! Buoyancy parameter for FTL (excluding g/thetav)
                  ! on theta level K-1
   vq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! Buoyancy parameter for FQW (excluding g/thetav)
                  ! on theta level K-1
   tv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! Virtual temperature on theta level K-1
   exner(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
                  ! exner function on theta level K-1
   gtr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! G/thetav on theta level K-1
   rhokh_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                  ! exchange coeffs for scalars on theta level K-1
   rhogamt_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
              2:bl_levels),                                                    &
                  ! counter gradient term for FTL on theta level K-1
   rhogamq_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
              2:bl_levels),                                                    &
                  ! counter gradient term for FQW on theta level K-1
   q1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! normalized excessive water from the saturation
                  ! on theta level K-1
   cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! cloud fraction derived by the bi-normal
                  ! distribution on theta level K-1
   ql(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                     &
      tke_levels),                                                             &
                  ! condensed liquid water derived by the bi-normal
                  ! distribution on theta level K-1
   wb_ng(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
                  ! buoyancy flux related to the skewness
                  ! on theta level K-1
   frac_shcu(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             tke_levels)
                  ! cloud fraction corrected by shallow cumulus
                  ! process on theta level K-1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MYM_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Monin-Obkhov length
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    r_mosurf(i,j)= -vkman*fb_surf(i,j)                                         &
                     / MAX(u_s(i,j)*u_s(i,j)*u_s(i,j), TINY(1.0))
  END DO
END DO

IF (l_first) THEN
  CALL mym_const_set

  ! IF the first value of e_trb has been set to be missing by the
  ! reconfiguration, the initialization for the whole domain
  ! is essential.
  IF (l_my_initialize .OR. ABS(e_trb(1, 1, 1) - rmdi) < real_eps) THEN
    IF (l_my_ini_zero) THEN
      DO k = 1, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            e_trb(i, j, k) = 0.0
            tsq_trb(i, j, k) = 0.0
            qsq_trb(i, j, k) = 0.0
            cov_trb(i, j, k) = 0.0
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
      ! Initialize the prognostic variables
      CALL mym_initialize(                                                     &
      ! IN levels
                  bl_levels,                                                   &
      ! IN fields
                  z_uv, z_tq, dbdz_l, dvdzm, dtldzm, dqwdzm,                   &
                  fqw, ftl, u_s, r_mosurf, fb_surf,                            &
      ! INOUT fields
                  e_trb, tsq_trb, qsq_trb, cov_trb)

              ! Above tke_levels, the prognostic variables should be zeros.
      DO k = tke_levels + 1, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            e_trb(i, j, k) = 0.0
            tsq_trb(i, j, k) = 0.0
            qsq_trb(i, j, k) = 0.0
            cov_trb(i, j, k) = 0.0
          END DO
        END DO
      END DO
    END IF  ! test if l_my_ini_zero
  END IF  ! test if l_my_initialize .OR. e_trb == rmdi

  IF (l_shcu_buoy) THEN
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (ABS(zhpar_shcu(i, j) - rmdi) < real_eps) THEN
          ! if missing has been set by the reconfiguration,
          ! it is replaced with z_tq(tke_levels-1).
          zhpar_shcu(i,j) = z_tq(i,j,tke_levels-1)
        END IF
      END DO
    END DO
  END IF
  ! need to initialise variables on every cycle as they will have been
  ! reset to mdi
  IF (cycleno == numcycles) l_first = .FALSE.
END IF  ! IF L_FIRST

! copy e_trb to qke (qke = 2 e_trb)
DO k = 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      qke(i, j, k) = 2.0 * e_trb(i, j, k)
    END DO
  END DO
END DO

! Set virtual temperature, exner function, and g / thetav
DO k = 2, tke_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      tv(i, j, k) = t(i, j, k - 1)                                             &
                   * (1.0 + c_virtual * q(i, j, k - 1)                         &
                         - qcl(i, j, k - 1) - qcf(i, j, k - 1))
      exner(i, j, k) =                                                         &
                (p_theta_levels(i, j, k - 1) / pref) ** kappa
      gtr(i, j, k) = g * exner(i, j, k) / tv(i, j, k)
    END DO
  END DO
END DO

IF (l_my_condense .OR. l_shcu_buoy) THEN
  CALL mym_condensation(                                                       &
  ! IN levels/switches
            bl_levels, levflag, nSCMDpkgs,L_SCMDiags,                          &
            BL_diag,                                                           &
  ! IN fields
            qw, tl, t, p_theta_levels, tsq_trb, qsq_trb, cov_trb,              &
  ! OUT fields
            vt, vq, q1, cld, ql)
END IF

IF (l_my_condense) THEN
  ! Re-evaluate DBDZ with the buoyancy parameters diagnosed by
  ! mym_condensation
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dbdz_l(i,j,k) = gtr(i, j, k)                                           &
                             * ( vt(i, j, k) * dtldzm(i, j, k) +               &
                                   vq(i, j, k) * dqwdzm(i, j, k) )

      END DO
    END DO
  END DO
ELSE
  ! Use the buoyancy parameters and DBDZ by the LS cloud scheme
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ! convert buoy params from the UM notation to the MY notaation
        vt(i, j, k) = bt_gb(i, j, k - 1) * tv(i, j, k)
        vq(i, j, k) = bq_gb(i, j, k - 1) * tv(i, j, k)                         &
                                             / exner(i, j, k)
        dbdz_l(i,j,k) = dbdz(i,j,k)
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

IF (l_shcu_buoy) THEN
  ! Evaluate the non-gradient buoyancy flux

  CALL mym_shcu_buoy(                                                          &
  ! IN levels/switches
             bl_levels,nSCMDpkgs,L_SCMDiags,                                   &
             BL_diag,                                                          &
  ! IN fields
             fb_surf, u_s, pstar,                                              &
             z_tq, z_uv, p_theta_levels, p_half,                               &
             u_p, v_p, t, q, qcl, qcf, q1, cld,                                &
  ! INOUT / OUT fields
             zhpar_shcu,frac_shcu, wb_ng)
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

IF (my_lowest_pd_surf > 0) THEN
  ! Calculate the gradient functions at the surface

  CALL mym_calcphi(                                                            &
        bl_levels, z_tq, r_mosurf, pmz, phh)
END IF

  ! Calculate diffusion coefficients and counter gradient terms,
  ! and integrate the prognostic variables.

CALL mym_turbulence(                                                           &
! IN levels/switches
        bl_levels, levflag, nSCMDpkgs,L_SCMDiags, BL_diag,                     &
! IN fields
        z_uv, z_tq,                                                            &
        vq, vt, gtr, fqw, ftl, wb_ng,                                          &
        dbdz_l, dtldzm, dqwdzm, dvdzm, dudz, dvdz,                             &
        r_mosurf, u_s, fb_surf, pmz, phh,                                      &
! INOUT fields
        qke, tsq_trb, qsq_trb, cov_trb, rhokm, rhokh_tq,                       &
! OUT fields
        rhogamu, rhogamv, rhogamt_tq, rhogamq_tq)

IF (l_3dtke) THEN
  DO k = 1, bl_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        visc_m(i,j,k) = rhokm(i,j,k+1)
        visc_h(i,j,k) = rhokh_tq(i,j,k+1)
      END DO
    END DO
  END DO
END IF

  ! multiply the density
DO k = 2, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      rhokm(i, j, k) = rho_wet_tq(i, j, k - 1) * rhokm(i, j, k)
      rhogamu(i, j, k) = rho_wet_tq(i, j, k - 1) * rhogamu(i, j, k)
      rhogamv(i, j, k) = rho_wet_tq(i, j, k - 1) * rhogamv(i, j, k)
    END DO
  END DO
END DO

! Note "RHO" here is always wet density (RHO_WET_TQ) so
! save multiplication of RHOKH to after interpolation
IF (.NOT. l_mr_physics) THEN
  DO k = 2, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        rhokh_tq(i, j, k) = rho_wet_tq(i, j, k - 1)                            &
                                             * rhokh_tq(i, j, k)
        rhogamt_tq(i, j, k) = rho_wet_tq(i, j, k - 1)                          &
                                             * rhogamt_tq(i, j, k)
        rhogamq_tq(i, j, k) = rho_wet_tq(i, j, k - 1)                          &
                                             * rhogamq_tq(i, j, k)
      END DO
    END DO
  END DO
END IF

! convert qke to e_trb
DO k = 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      e_trb(i, j, k) = 0.5 * qke(i, j, k)
    END DO
  END DO
END DO

! Interpolate RHOKH_TQ, RHOGAMT_TQ and RHOGAMQ_TQ on theta levels
! to rho levels.
DO k = 2, tke_levels - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      r_weight1 = 1.0 / (r_theta_levels(i,j,k) -                               &
                                   r_theta_levels(i,j, k-1))
      weight2 = (r_theta_levels(i,j,k) -                                       &
                          r_rho_levels(i,j,k)) * r_weight1
      weight3 = (r_rho_levels(i,j,k) -                                         &
                          r_theta_levels(i,j,k-1)) * r_weight1
      rhokh(i,j,k) =                                                           &
                          weight3 * rhokh_tq(i,j,k+1)                          &
                         +weight2 * rhokh_tq(i,j,k)
      rhogamt(i,j,k) =                                                         &
                          weight3 * rhogamt_tq(i,j,k+1)                        &
                         +weight2 * rhogamt_tq(i,j,k)
      rhogamq(i,j,k) =                                                         &
                          weight3 * rhogamq_tq(i,j,k+1)                        &
                         +weight2 * rhogamq_tq(i,j,k)
    END DO
  END DO
END DO

k = tke_levels
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    r_weight1 = 1.0 / (r_theta_levels(i,j,k) -                                 &
                                 r_theta_levels(i,j, k-1))
    weight2 = (r_theta_levels(i,j,k) -                                         &
                        r_rho_levels(i,j,k)) * r_weight1
    weight3 = (r_rho_levels(i,j,k) -                                           &
                        r_theta_levels(i,j,k-1)) * r_weight1
    rhokh(i,j,k) = weight2 * rhokh_tq(i,j,k)
    rhogamt(i,j,k) = weight2 * rhogamt_tq(i,j,k)
    rhogamq(i,j,k) = weight2 * rhogamq_tq(i,j,k)
  END DO
END DO

! Finally multiply RHOKH by dry density
IF (l_mr_physics) THEN
  DO k = 2, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        rhokh(i, j, k) = rho_mix(i, j, k) * rhokh(i, j, k)
        rhogamt(i, j, k) = rho_mix(i, j, k) * rhogamt(i, j, k)
        rhogamq(i, j, k) = rho_mix(i, j, k) * rhogamq(i, j, k)
      END DO
    END DO
  END DO
END IF

! Above tke_levels, fluxes should be zero.
DO k = tke_levels + 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      rhokh(i, j, k) = 0.0
      rhogamt(i, j, k) = 0.0
      rhogamq(i, j, k) = 0.0
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE mym_ctl
END MODULE mym_ctl_mod
