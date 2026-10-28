!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Do conversion between physical wind on faces and prognostic u in LFRic
!> @details
module convert_w2_wind_strict_kernel_mod

  use argument_mod,                  only: arg_type, func_type,                &
                                           GH_FIELD, GH_SCALAR,                &
                                           GH_INTEGER, GH_REAL, GH_LOGICAL,    &
                                           GH_READ, GH_WRITE,                  &
                                           CELL_COLUMN,                        &
                                           ANY_SPACE_9,                        &
                                           ANY_DISCONTINUOUS_SPACE_1,          &
                                           ANY_DISCONTINUOUS_SPACE_3,          &
                                           GH_BASIS, GH_DIFF_BASIS,            &
                                           GH_EVALUATOR
  use constants_mod,                 only: r_def, i_def, l_def
  use fs_continuity_mod,             only: W2
  use kernel_mod,                    only: kernel_type
  use reference_element_mod,         only: W, S, E, N, B, T
  use sci_native_jacobian_mod,       only: native_jacobian
  use sci_chi_transform_mod,         only: chi2xyz
  use base_mesh_config_mod,          only: geometry_spherical

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: convert_w2_wind_strict_kernel_type
    private
    type(arg_type) :: meta_args(14) = (/                                       &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2),                        & ! u_out
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! u_in
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! dA
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_SPACE_9),               & ! chi1
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_SPACE_9),               & ! chi2
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_SPACE_9),               & ! chi3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! panel_id
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ew
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ns
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             & ! coord_system
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             & ! geometry
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             & ! topology
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             & ! scaled_radius
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              & ! prog_to_phys
         /)
    type(func_type) :: meta_funcs(2) = (/                                      &
         func_type(W2, GH_BASIS),                                              &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                       &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: convert_w2_wind_strict_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: convert_w2_wind_strict_code

contains

!> @brief Conversion between physical wind on faces and prognostic u in LFRic
!! @param[in]     nlayers             Number of layers
!! @param[in,out] u_out               Output W2 wind field
!! @param[in]     u_in                Input W2 wind field
!! @param[in]     dA                  Area element in W2
!! @param[in]     chi1                First coordinate field in Wchi
!! @param[in]     chi2                Second coordinate field in Wchi
!! @param[in]     chi3                Third coordinate field in Wchi
!! @param[in]     panel_id            The ID number of the current panel
!! @param[in]     face_selector_ew    Field indicating which W/E faces to loop
!! @param[in]     face_selector_ns    Field indicating which N/S faces to loop
!! @param[in]     coord_system        Coordinate system enumeration
!! @param[in]     geometry            Mesh geometry enumeration
!! @param[in]     topology            Mesh topology enumeration
!! @param[in]     scaled_radius       Scaled planetary radius
!! @param[in]     prog_to_phys        Switch for conversion direction
!! @param[in]     ndf_w2              Number of DOFs for W2 space
!! @param[in]     undf_w2             Number of unique DOFs for W2 space
!! @param[in]     map_w2              DoF map for W2
!! @param[in]     basis_w2            W2 basis functions evaluated at W2 DoFs
!! @param[in]     ndf_wx              Number of DOFs for Wchi space
!! @param[in]     undf_wx             Number of unique DOFs for Wchi space
!! @param[in]     map_wx              DoF map for Wchi
!! @param[in]     basis_wx            Wchi basis functions evaluated at W2 DoFs
!! @param[in]     diff_basis_wx       Differential basis functions for Wchi
!!                                    evaluated at W2 DoFs
!! @param[in]     ndf_pid             Number of DOFs for pid space
!! @param[in]     undf_pid            Number of unique DOFs for pid space
!! @param[in]     map_pid             Dofmap for pid space
!! @param[in]     ndf_2d              Number of DOFs for 2D W3 space
!! @param[in]     undf_2d             Number of unique DOFs for 2D W3 space
!! @param[in]     map_2d              DoF map for 2D W3
subroutine convert_w2_wind_strict_code(nlayers,                                &
                                       u_out,                                  &
                                       u_in,                                   &
                                       dA,                                     &
                                       chi1,                                   &
                                       chi2,                                   &
                                       chi3,                                   &
                                       panel_id,                               &
                                       face_selector_ew,                       &
                                       face_selector_ns,                       &
                                       coord_system,                           &
                                       geometry,                               &
                                       topology,                               &
                                       scaled_radius,                          &
                                       prog_to_phys,                           &
                                       ndf_w2, undf_w2, map_w2,                &
                                       basis_w2,                               &
                                       ndf_wx, undf_wx, map_wx,                &
                                       basis_wx,                               &
                                       diff_basis_wx,                          &
                                       ndf_pid, undf_pid, map_pid,             &
                                       ndf_2d, undf_2d, map_2d                 &
                                      )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in) :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
  integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

  real(kind=r_def),    dimension(undf_w2),  intent(inout) :: u_out
  real(kind=r_def),    dimension(undf_w2),  intent(in)    :: u_in, dA
  real(kind=r_def),    dimension(undf_wx),  intent(in)    :: chi1, chi2, chi3
  real(kind=r_def),    dimension(undf_pid), intent(in)    :: panel_id
  integer(kind=i_def), dimension(undf_2d),  intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_2d),  intent(in)    :: face_selector_ns
  integer(kind=i_def),                      intent(in)    :: coord_system
  integer(kind=i_def),                      intent(in)    :: geometry
  integer(kind=i_def),                      intent(in)    :: topology
  real(kind=r_def),                         intent(in)    :: scaled_radius
  logical(kind=l_def),                      intent(in)    :: prog_to_phys
  real(kind=r_def), intent(in) :: basis_w2(3,ndf_w2,ndf_w2)
  real(kind=r_def), intent(in) :: basis_wx(1,ndf_wx,ndf_w2)
  real(kind=r_def), intent(in) :: diff_basis_wx(3,ndf_wx,ndf_w2)

  ! Internal variables
  integer(kind=i_def) :: i, j, k, kk, df, bt, ipanel
  real(kind=r_def)    :: chi1_e(ndf_wx), chi2_e(ndf_wx)
  real(kind=r_def)    :: chi3_e(nlayers,ndf_wx)
  real(kind=r_def)    :: jacobian(nlayers,3,3,2)
  real(kind=r_def)    :: dj(nlayers,2)
  real(kind=r_def)    :: coords(3)
  real(kind=r_def)    :: x, y, z, one_ov_r
  real(kind=r_def)    :: ez(3), v_comp(3), v_phys(3), v_base(3)
  real(kind=r_def)    :: vz_base
  real(kind=r_def)    :: vz_base_wt(nlayers+1)
  real(kind=r_def)    :: rmult_wt(nlayers+1)

  integer(kind=i_def), parameter :: horizontal_dofs(4) = (/W, S, E, N/)
  integer(kind=i_def), parameter :: vertical_dofs(2) = (/B, T/)

  ipanel = int(panel_id(map_pid(1)), i_def)

  rmult_wt(1) = 1.0_r_def
  rmult_wt(2:nlayers) = 0.5_r_def
  rmult_wt(nlayers+1) = 1.0_r_def

  ! Calculate metrics at the cell bottom
  do df = 1, ndf_wx
    chi1_e(df) = chi1(map_wx(df))
    chi2_e(df) = chi2(map_wx(df))
    do k = 1, nlayers
      chi3_e(k,df) = chi3(map_wx(df)+k-1)
    end do
  end do
  ! loop for cell bottom & top
  do kk = 1, 2
    bt = vertical_dofs(kk)
    call native_jacobian(coord_system, geometry, topology, scaled_radius,      &
                         ndf_wx, nlayers, chi1_e, chi2_e, chi3_e, ipanel,      &
                         basis_wx(1,1,bt),  diff_basis_wx(1,1,bt),             &
                         jacobian(1,1,1,kk), dj(1,kk))
  end do

  ! Calculate unit vector parallel to gravity
  if (geometry == geometry_spherical) then
    coords(:) = 0.0_r_def
    do df = 1, ndf_wx
      coords(1) = coords(1) + chi1(map_wx(df))*basis_wx(1,df,B)
      coords(2) = coords(2) + chi2(map_wx(df))*basis_wx(1,df,B)
      coords(3) = coords(3) + chi3(map_wx(df))*basis_wx(1,df,B)
    end do
    call chi2xyz( coords(1), coords(2), coords(3), ipanel, x, y, z )
    ! normalize position vector
    one_ov_r = 1.0_r_def / sqrt(x**2 + y**2 + z**2)
    ez(1) = x * one_ov_r
    ez(2) = y * one_ov_r
    ez(3) = z * one_ov_r
  else
    ez(1) = 0.0_r_def
    ez(2) = 0.0_r_def
    ez(3) = 1.0_r_def
  end if

  if (prog_to_phys) then
    ! Calculate horizontal component
    do j = 1, face_selector_ew(map_2d(1)) + face_selector_ns(map_2d(1))
      df = j
      if (j == 3 .and. face_selector_ns(map_2d(1)) == 2 &
                 .and. face_selector_ew(map_2d(1)) == 1) df = N
      do k = 1, nlayers
        u_out(map_w2(df)+k-1) = u_in(map_w2(df)+k-1) / dA(map_w2(df)+k-1)
      end do
    end do

    ! Calculate vertical component

    u_out(map_w2(B):map_w2(B)+nlayers) = 0.0_r_def
    ! loop for cell bottom & top
    do kk = 1, 2
      bt = vertical_dofs(kk)
      do k = 1, nlayers
        ! 3D vector field in computation space
        v_comp(:) = 0.0_r_def
        do df = 1, ndf_w2
          v_comp(:) = v_comp(:) + u_in(map_w2(df)+k-1)*basis_w2(:,df,bt)
        end do
        ! 3D vector field in physical space
        v_phys(:) = 0.0_r_def
        do j = 1, 3
          v_phys(:) = v_phys(:) + jacobian(k,:,j,kk)*v_comp(j)
        end do
        v_phys(:) = v_phys(:) / dj(k,kk)
        ! Calculate vertical component, using inner product
        u_out(map_w2(bt)+k-1) = u_out(map_w2(bt)+k-1) + &
             rmult_wt(kk+k-1) * ( v_phys(1)*ez(1) + &
                                  v_phys(2)*ez(2) + &
                                  v_phys(3)*ez(3) )
      end do
    end do

  else ! from physical to prognostic u

    ! Calculate horizontal component
    do j = 1, face_selector_ew(map_2d(1)) + face_selector_ns(map_2d(1))
      df = j
      if (j == 3 .and. face_selector_ns(map_2d(1)) == 2 &
                 .and. face_selector_ew(map_2d(1)) == 1) df = N
      do k = 1, nlayers
        u_out(map_w2(df)+k-1) = u_in(map_w2(df)+k-1) * dA(map_w2(df)+k-1)
      end do
    end do

    ! Calculate vertical component

    u_out(map_w2(B):map_w2(B)+nlayers) = u_in(map_w2(B):map_w2(B)+nlayers)
    vz_base_wt(:) = 0.0_r_def
    ! loop for cell bottom & top
    do kk = 1, 2
      bt = vertical_dofs(kk)
      do k = 1, nlayers
        do i = 1, 4
          df = horizontal_dofs(i)
          ! Basis function evaluated at cell bottom or top
          v_base(:) = 0.0_r_def
          do j = 1, 3
            v_base(:) = v_base(:) + jacobian(k,:,j,kk)*basis_w2(j,df,bt)
          end do
          v_base(:) = v_base(:) / dj(k,kk)
          ! Calculate vertical component, using inner product
          vz_base = v_base(1) * ez(1) + v_base(2) * ez(2) + v_base(3) * ez(3)

          u_out(map_w2(bt)+k-1) = u_out(map_w2(bt)+k-1) - &
              rmult_wt(kk+k-1) * u_in(map_w2(df)+k-1) * &
              dA(map_w2(df)+k-1) * vz_base
        end do

        ! Calculate basis function evaluated at cell bottom or top
        v_base(:) = 0.0_r_def
        do j = 1, 3
          v_base(:) = v_base(:) + jacobian(k,:,j,kk) * basis_w2(j,bt,bt)
        end do
        v_base(:) = v_base(:) / dj(k,kk)
        vz_base = v_base(1) * ez(1) + v_base(2) * ez(2) + v_base(3) * ez(3)
        vz_base_wt(kk+k-1) = vz_base_wt(kk+k-1) + &
            rmult_wt(kk+k-1) * vz_base
      end do
    end do

    u_out(map_w2(B):map_w2(B)+nlayers) = &
         u_out(map_w2(B):map_w2(B)+nlayers) / vz_base_wt(1:nlayers+1)

  end if

end subroutine convert_w2_wind_strict_code

end module convert_w2_wind_strict_kernel_mod
