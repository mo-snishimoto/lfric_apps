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
                                           GH_DIFF_BASIS,                      &
                                           GH_EVALUATOR
  use constants_mod,                 only: r_def, i_def, l_def
  use fs_continuity_mod,             only: W2
  use kernel_mod,                    only: kernel_type
  use reference_element_mod,         only: W, S, E, N, B, T

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: convert_w2_wind_strict_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2),                        & ! u_out
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! u_in
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! detj_w2
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! dA
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_SPACE_9),               & ! chi3
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! face_selector_ew
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! face_selector_ns
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              & ! prog_to_phys
         /)
    type(func_type) :: meta_funcs(1) = (/                                      &
         func_type(ANY_SPACE_9, GH_DIFF_BASIS)                                 &
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
!! @param[in]     detj_w2             Jacobian at W2 point
!! @param[in]     dA                  Area element in W2
!! @param[in]     chi3                Third coordinate field in Wchi
!! @param[in]     face_selector_ew    Field indicating which W/E faces to loop
!! @param[in]     face_selector_ns    Field indicating which N/S faces to loop
!! @param[in]     prog_to_phys        Switch for conversion direction
!! @param[in]     ndf_w2              Number of DOFs for W2 space
!! @param[in]     undf_w2             Number of unique DOFs for W2 space
!! @param[in]     map_w2              DoF map for W2
!! @param[in]     ndf_wx              Number of DOFs for Wchi space
!! @param[in]     undf_wx             Number of unique DOFs for Wchi space
!! @param[in]     map_wx              DoF map for Wchi
!! @param[in]     diff_basis_wx       Differential basis functions for Wchi
!!                                    evaluated at W2 DoFs
!! @param[in]     ndf_2d              Number of DOFs for 2D W3 space
!! @param[in]     undf_2d             Number of unique DOFs for 2D W3 space
!! @param[in]     map_2d              DoF map for 2D W3
subroutine convert_w2_wind_strict_code(nlayers,                                &
                                       u_out,                                  &
                                       u_in,                                   &
                                       detj_w2,                                &
                                       dA,                                     &
                                       chi3,                                   &
                                       face_selector_ew,                       &
                                       face_selector_ns,                       &
                                       prog_to_phys,                           &
                                       ndf_w2, undf_w2, map_w2,                &
                                       ndf_wx, undf_wx, map_wx,                &
                                       diff_basis_wx,                          &
                                       ndf_2d, undf_2d, map_2d                 &
                                      )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in) :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
  integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

  real(kind=r_def),    dimension(undf_w2),  intent(inout) :: u_out
  real(kind=r_def),    dimension(undf_w2),  intent(in)    :: u_in, detj_w2, dA
  real(kind=r_def),    dimension(undf_wx),  intent(in)    :: chi3
  integer(kind=i_def), dimension(undf_2d),  intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_2d),  intent(in)    :: face_selector_ns
  logical(kind=l_def),                      intent(in)    :: prog_to_phys
  real(kind=r_def), intent(in) :: diff_basis_wx(3,ndf_wx,ndf_w2)

  ! Internal variables
  integer(kind=i_def) :: j, k, kk, df, bt, idxs, idxe
  real(kind=r_def), dimension(nlayers,2) :: dzdx1, dzdx2, dzdx3, detj
  real(kind=r_def), dimension(nlayers+1) :: rmult_wt, vz_base_wt

  integer(kind=i_def), parameter :: vertical_dofs(2) = (/B, T/)

  rmult_wt(1) = 1.0_r_def
  rmult_wt(2:nlayers) = 0.5_r_def
  rmult_wt(nlayers+1) = 1.0_r_def

  ! Calculate dz/dx[1-3] at the cell bottom & top
  dzdx1(:,:) = 0.0_r_def
  dzdx2(:,:) = 0.0_r_def
  dzdx3(:,:) = 0.0_r_def
  do kk = 1, 2
    bt = vertical_dofs(kk)
    do df = 1, ndf_wx
      idxs = map_wx(df)
      idxe = map_wx(df) + nlayers - 1
      dzdx1(:, kk) = dzdx1(:, kk) + &
        chi3(idxs:idxe)*diff_basis_wx(1, df, bt)
      dzdx2(:, kk) = dzdx2(:, kk) + &
        chi3(idxs:idxe)*diff_basis_wx(2, df, bt)
      dzdx3(:, kk) = dzdx3(:, kk) + &
        chi3(idxs:idxe)*diff_basis_wx(3, df, bt)
    end do
  end do

  ! Calculate jacobian at the bottom
  detj(1,1) = detj_w2(map_w2(B))
  do k = 2, nlayers
    detj(k,1) = detj_w2(map_w2(B)+k-1) * &
                2.0_r_def * dzdx3(k,1) / (dzdx3(k,1) + dzdx3(k-1,2))
  end do

  ! Calculate jacobian at the top
  do k = 1, nlayers - 1
    detj(k,2) = detj_w2(map_w2(T)+k-1) * &
                2.0_r_def * dzdx3(k,2) / (dzdx3(k+1,1) + dzdx3(k,2))
  end do
  detj(nlayers,2) = detj_w2(map_w2(T)+nlayers-1)

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
        u_out(map_w2(bt)+k-1) = u_out(map_w2(bt)+k-1) + ( &
          0.5_r_def*(u_in(map_w2(W)+k-1) + u_in(map_w2(E)+k-1))*dzdx1(k,kk) - &
          0.5_r_def*(u_in(map_w2(S)+k-1) + u_in(map_w2(N)+k-1))*dzdx2(k,kk) + &
                                           u_in(map_w2(bt)+k-1)*dzdx3(k,kk) ) &
          / detj(k,kk)
      end do
    end do
    u_out(map_w2(B):map_w2(B)+nlayers) = &
      rmult_wt(:) * u_out(map_w2(B):map_w2(B)+nlayers)

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
        u_out(map_w2(bt)+k-1) = u_out(map_w2(bt)+k-1) -                &
          rmult_wt(kk+k-1) * 0.5_r_def * (                             &
            ( u_in(map_w2(W)+k-1)*dA(map_w2(W)+k-1) +                  &
              u_in(map_w2(E)+k-1)*dA(map_w2(E)+k-1)                    &
            ) * dzdx1(k,kk) -                                          &
            ( u_in(map_w2(S)+k-1)*dA(map_w2(S)+k-1) +                  &
              u_in(map_w2(N)+k-1)*dA(map_w2(N)+k-1)                    &
            ) * dzdx2(k,kk)                                            &
          ) / detj(k,kk)
        vz_base_wt(k+kk-1) = vz_base_wt(k+kk-1) +                      &
          rmult_wt(kk+k-1) * dzdx3(k,kk) / detj(k,kk)
      end do
    end do

    u_out(map_w2(B):map_w2(B)+nlayers) = &
      u_out(map_w2(B):map_w2(B)+nlayers) / vz_base_wt(1:nlayers+1)

  end if

end subroutine convert_w2_wind_strict_code

end module convert_w2_wind_strict_kernel_mod
