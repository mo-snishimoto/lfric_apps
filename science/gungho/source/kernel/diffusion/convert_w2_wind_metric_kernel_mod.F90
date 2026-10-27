!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Do conversion between physical wind on faces and prognostic u in LFRic
!> @details
module convert_w2_wind_metric_kernel_mod

  use argument_mod,                  only: arg_type, func_type,                &
                                           GH_FIELD, GH_SCALAR,                &
                                           GH_INTEGER, GH_REAL, GH_LOGICAL,    &
                                           GH_READ, GH_WRITE,                  &
                                           CELL_COLUMN,                        &
                                           ANY_SPACE_9,                        &
                                           ANY_DISCONTINUOUS_SPACE_1,          &
                                           GH_DIFF_BASIS, GH_EVALUATOR
  use constants_mod,                 only: r_def, i_def, l_def
  use fs_continuity_mod,             only: W2
  use kernel_mod,                    only: kernel_type
  use reference_element_mod,         only: W, S, E, N, B
  use sci_face_selector_support_mod, only: face_from_face_selector

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: convert_w2_wind_metric_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2),                        & ! u_out
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! u_in
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! dA
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_SPACE_9),               & ! chi3
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! face selector ew
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! face selector ns
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              & ! prog_to_phys
         /)
    type(func_type) :: meta_funcs(1) = (/                                      &
         func_type(ANY_SPACE_9, GH_DIFF_BASIS)                                 &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: convert_w2_wind_metric_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: convert_w2_wind_metric_code

contains

!> @brief Conversion between physical wind on faces and prognostic u in LFRic
!! @param[in]     nlayers             Number of layers
!! @param[in,out] u_out               Output W2 wind field
!! @param[in]     u_in                Input W2 wind field
!! @param[in]     dA                  Area element in W2
!! @param[in]     chi3                Third coordinate field in Wchi
!! @param[in]     face_selector_ew    Field indicating which W/E faces to loop
!! @param[in]     face_selector_ns    Field indicating which N/S faces to loop
!> @param[in]     prog_to_phys        Switch for conversion direction
!! @param[in]     ndf_w2              Number of DOFs for W2 space
!! @param[in]     undf_w2             Number of unique DOFs for W2 space
!! @param[in]     map_w2              DoF map for W2
!! @param[in]     ndf_wx              Number of DOFs for Wchi space
!! @param[in]     undf_wx             Number of unique DOFs for Wchi space
!! @param[in]     map_wx              DoF map for Wchi
!! @param[in]     diff_basis_wx_on_w2 Differential basis functions for Wchi
!!                                    evaluated at W2 points
!! @param[in]     ndf_2d              Number of DOFs for 2D W3 space
!! @param[in]     undf_2d             Number of unique DOFs for 2D W3 space
!! @param[in]     map_2d              DoF map for 2D W3
subroutine convert_w2_wind_metric_code(nlayers,                                &
                                       u_out,                                  &
                                       u_in,                                   &
                                       dA,                                     &
                                       chi3,                                   &
                                       face_selector_ew,                       &
                                       face_selector_ns,                       &
                                       prog_to_phys,                           &
                                       ndf_w2, undf_w2, map_w2,                &
                                       ndf_wx, undf_wx, map_wx,                &
                                       diff_basis_wx_on_w2,                    &
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

  real(kind=r_def),    dimension(undf_w2), intent(inout) :: u_out
  real(kind=r_def),    dimension(undf_w2), intent(in)    :: u_in, dA
  real(kind=r_def),    dimension(undf_wx), intent(in)    :: chi3
  integer(kind=i_def), dimension(undf_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_2d), intent(in)    :: face_selector_ns
  logical(kind=l_def),                     intent(in)    :: prog_to_phys
  real(kind=r_def), intent(in) :: diff_basis_wx_on_w2(3,ndf_wx,ndf_w2)

  ! Internal variables
  integer(kind=i_def) :: j, k, df
  real(kind=r_def)    :: u_at_w, v_at_w
  real(kind=r_def), dimension(nlayers) :: dzdx1, dzdx2, dzdx3

  ! Differential of physical altitude in each direction evaluated at cell bottom
  dzdx1 = 0.0_r_def
  dzdx2 = 0.0_r_def
  dzdx3 = 0.0_r_def
  do df = 1, ndf_wx
    do k = 1, nlayers
      dzdx1(k) = dzdx1(k) + chi3(map_wx(df)+k-1)*diff_basis_wx_on_w2(1, df, B)
      dzdx2(k) = dzdx2(k) + chi3(map_wx(df)+k-1)*diff_basis_wx_on_w2(2, df, B)
      dzdx3(k) = dzdx3(k) + chi3(map_wx(df)+k-1)*diff_basis_wx_on_w2(3, df, B)
    end do
  end do

  if (prog_to_phys) then
    do j = 1, ABS(face_selector_ew(map_2d(1))) + ABS(face_selector_ns(map_2d(1)))
      df = face_from_face_selector(j, face_selector_ew(map_2d(1)), &
                                      face_selector_ns(map_2d(1)))
      do k = 0, nlayers-1
        u_out(map_w2(df)+k) = u_in(map_w2(df)+k) / dA(map_w2(df)+k)
      end do
    end do

    ! Assume no-slip at the ground
    u_out(map_w2(B)) = 0.0_r_def

    do k = 1, nlayers-1
      u_at_w = 0.25_r_def * ( u_in(map_w2(W)+k) + u_in(map_w2(W)+k-1) + &
                              u_in(map_w2(E)+k) + u_in(map_w2(E)+k-1) )
      v_at_w = 0.25_r_def * ( u_in(map_w2(S)+k) + u_in(map_w2(S)+k-1) + &
                              u_in(map_w2(N)+k) + u_in(map_w2(N)+k-1) )
      u_out(map_w2(B)+k) = &
        ( u_in(map_w2(B)+k) + &
          ( u_at_w*dzdx1(k+1) - v_at_w*dzdx2(k+1) ) / dzdx3(k+1) &
        ) / dA(map_w2(B)+k)
    end do
    u_out(map_w2(B)+nlayers) = 0.0_r_def
  else
    do j = 1, ABS(face_selector_ew(map_2d(1))) + ABS(face_selector_ns(map_2d(1)))
      df = face_from_face_selector(j, face_selector_ew(map_2d(1)), &
                                      face_selector_ns(map_2d(1)))
      do k = 0, nlayers-1
        u_out(map_w2(df)+k) = u_in(map_w2(df)+k) * dA(map_w2(df)+k)
      end do
    end do

    u_out(map_w2(B)) = 0.0_r_def
    do k = 1, nlayers-1
      u_at_w = 0.25_r_def * ( u_in(map_w2(W)+k)   * dA(map_w2(W)+k)   + &
                              u_in(map_w2(W)+k-1) * dA(map_w2(W)+k-1) + &
                              u_in(map_w2(E)+k)   * dA(map_w2(E)+k)   + &
                              u_in(map_w2(E)+k-1) * dA(map_w2(E)+k-1) )
      v_at_w = 0.25_r_def * ( u_in(map_w2(S)+k)   * dA(map_w2(S)+k)     + &
                              u_in(map_w2(S)+k-1) * dA(map_w2(S)+k-1) + &
                              u_in(map_w2(N)+k)   * dA(map_w2(N)+k)   + &
                              u_in(map_w2(N)+k-1) * dA(map_w2(N)+k-1) )
      u_out(map_w2(B)+k) = u_in(map_w2(B)+k) * dA(map_w2(B)+k) + &
                        ( -u_at_w*dzdx1(k+1) + v_at_w*dzdx2(k+1) ) / dzdx3(k+1)
    end do
    u_out(map_w2(B)+nlayers) = 0.0_r_def
  end if

end subroutine convert_w2_wind_metric_code

end module convert_w2_wind_metric_kernel_mod
