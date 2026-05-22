!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Check dofmap of W0, W1 and W2 space.
!>
module check_dof_relation_kernel_mod

  use argument_mod,      only : arg_type,                    &
                                GH_FIELD, GH_REAL,           &
                                GH_READ, GH_WRITE,           &
                                STENCIL, CROSS, CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W0, W1, W2
  use kernel_mod,        only : kernel_type

  implicit none
  private

  type, public, extends(kernel_type) :: check_dof_relation_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W0, STENCIL(CROSS)),         &! data_w0
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W1, STENCIL(CROSS)),         &! data_w1
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W2, STENCIL(CROSS)),         &! data_w2
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE,  W2)                          &
! data_w2_dummy
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: check_dof_relation_code
  end type

  public :: check_dof_relation_code

contains

!> @brief check_dof_relation_kernel_mod
!> @param[in]     nlayers       Number of layers in the mesh
!> @param[in]     data_w0       W0 data
!> @param[in]     smap_w0_size  Size of the stencil map for data_w0
!> @param[in]     smap_w0       Stencil map for data_w0
!> @param[in]     data_w1       W1 data
!> @param[in]     smap_w1_size  Size of the stencil map for data_w1
!> @param[in]     smap_w1       Stencil map for data_w1
!> @param[in]     data_w2       W2 data
!> @param[in]     smap_w2_size  Size of the stencil map for data_w2
!> @param[in]     smap_w2       Stencil map for data_w2
!> @param[in]     data_w2_dummy W2 data
!> @param[in]     ndf_w0        Number of DOFs for W0 space
!> @param[in]     undf_w0       Number of unique DOFs for W0 space
!> @param[in]     map_w0        Dofmap for the cell at the base of the column
!> @param[in]     ndf_w1        Number of DOFs for W1 space
!> @param[in]     undf_w1       Number of unique DOFs for W1 space
!> @param[in]     map_w1        Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2        Number of DOFs for W2 space
!> @param[in]     undf_w2       Number of unique DOFs for W2 space
!> @param[in]     map_w2        Dofmap for the cell at the base of the column
subroutine check_dof_relation_code( nlayers,                                   &
                                    data_w0,                                   &
                                    smap_w0_size,                              &
                                    smap_w0,                                   &
                                    data_w1,                                   &
                                    smap_w1_size,                              &
                                    smap_w1,                                   &
                                    data_w2,                                   &
                                    smap_w2_size,                              &
                                    smap_w2,                                   &
                                    data_w2_dummy,                             &
                                    ndf_w0, undf_w0, map_w0,                   &
                                    ndf_w1, undf_w1, map_w1,                   &
                                    ndf_w2, undf_w2, map_w2                    &
                                   )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w0, undf_w0
  integer(kind=i_def), intent(in) :: map_w0(ndf_w0)
  integer(kind=i_def), intent(in) :: ndf_w1, undf_w1
  integer(kind=i_def), intent(in) :: map_w1(ndf_w1)
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: smap_w0_size
  integer(kind=i_def), intent(in) :: smap_w0(ndf_w0,smap_w0_size)
  integer(kind=i_def), intent(in) :: smap_w1_size
  integer(kind=i_def), intent(in) :: smap_w1(ndf_w1,smap_w1_size)
  integer(kind=i_def), intent(in) :: smap_w2_size
  integer(kind=i_def), intent(in) :: smap_w2(ndf_w2,smap_w2_size)
  real(kind=r_def), dimension(undf_w0), intent(in) :: data_w0
  real(kind=r_def), dimension(undf_w1), intent(in) :: data_w1
  real(kind=r_def), dimension(undf_w2), intent(in) :: data_w2
  real(kind=r_def), dimension(undf_w2), intent(out) :: data_w2_dummy

  ! Internal variables
  integer(kind=i_def) :: df, st, df2
  integer(kind=i_def) :: w0_st(undf_w0), w0_df(undf_w0)
  integer(kind=i_def) :: w1_st(undf_w1), w1_df(undf_w1)
  integer(kind=i_def) :: w2_st(undf_w2), w2_df(undf_w2)

  w0_st(:) = -1
  w1_st(:) = -1
  w2_st(:) = -1
  w0_df(:) = -1
  w1_df(:) = -1
  w2_df(:) = -1

  ! Assumed direction for derivatives in this kernel is:
  !  y
  !  ^
  !  |_> x
  !

  ! Layout of dofs for the stencil map
  ! dimensions of map are (ndf, ncell)
  ! Horizontally:
  !
  !   -- 4 --
  !   |     |
  !   1     3
  !   |     |
  !   -- 2 --
  !
  ! df = 5 is in the centre on the bottom face
  ! df = 6 is in the centre on the top face

  ! The layout of the cells in the stencil is:
  !
  !          -----
  !          |   |
  !          | 5 |
  !     ---------------
  !     |    |   |    |
  !     |  2 | 1 |  4 |
  !     ---------------
  !          |   |
  !          | 3 |
  !          -----

  ! If the full stencil isn't available, we must be at the domain edge.
  if (smap_w2_size < 5_i_def) then
    return
  end if

  ! Check W0
  do df = 1, ndf_w0
    search_w0: do st = 2, 5
      do df2 = 1, ndf_w0
        if (map_w0(df) == smap_w0(df2, st)) then
          w0_st(df) = st
          w0_df(df) = df2
          exit search_w0
        end if
      end do
    end do search_w0
  end do

  ! Check W1
  do df = 1, ndf_w1
    search_w1: do st = 2, 5
      do df2 = 1, ndf_w1
        if (map_w1(df) == smap_w1(df2, st)) then
          w1_st(df) = st
          w1_df(df) = df2
          exit search_w1
        end if
      end do
    end do search_w1
  end do

  ! Check W2
  do df = 1, ndf_w2
    search_w2: do st = 2, 5
      do df2 = 1, ndf_w2
        if (map_w2(df) == smap_w2(df2, st)) then
          w2_st(df) = st
          w2_df(df) = df2
          exit search_w2
        end if
      end do
    end do search_w2
  end do

  do df = 1, ndf_w0
    write(6,'(a,i2,a,i2,a,i2,a)') 'In W0 space, ', df, 'th dof is same as ', w0_df(df), 'th dof in ', w0_st(df), 'th stencil cell.'
  end do
  do df = 1, ndf_w1
    write(6,'(a,i2,a,i2,a,i2,a)') 'In W1 space, ', df, 'th dof is same as ', w1_df(df), 'th dof in ', w1_st(df), 'th stencil cell.'
  end do
  do df = 1, ndf_w2
    write(6,'(a,i2,a,i2,a,i2,a)') 'In W2 space, ', df, 'th dof is same as ', w2_df(df), 'th dof in ', w2_df(df), 'th stencil cell.'
  end do

  stop 55

end subroutine check_dof_relation_code

end module check_dof_relation_kernel_mod

