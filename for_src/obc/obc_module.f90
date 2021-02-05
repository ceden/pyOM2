


module obc_module
  implicit none
  logical :: enable_obc_north = .false.
  logical :: enable_obc_south = .false.
  logical :: enable_prescribe_psi_obc_south = .false.
  logical :: enable_prescribe_psi_obc_north = .false.
  logical :: enable_restore_TS_obc_south = .false.
  logical :: enable_restore_TS_obc_north = .false.
  real*8, allocatable :: psi_wall_north(:),psi_wall_south(:)
  real*8, allocatable :: temp_wall_north(:,:),temp_wall_south(:,:)
  real*8, allocatable :: salt_wall_north(:,:),salt_wall_south(:,:)
  real*8 :: obc_tscl = 0.0
  real*8 :: obc_K_h  = 0.0
  logical :: enable_obc_north_damping = .false.
  logical :: enable_obc_south_damping = .false.
  real*8 :: obc_north_damping_len = 1.
  real*8 :: obc_south_damping_len = 1.
  real*8 :: obc_north_damping_amp = 0.
  real*8 :: obc_south_damping_amp = 0.
end module obc_module

subroutine allocate_obc_module
  use main_module
  use obc_module
  if (enable_obc_north) then
     allocate(psi_wall_north(is_pe-onx:ie_pe+onx) ); psi_wall_north=0.0
     allocate(temp_wall_north(is_pe-onx:ie_pe+onx,nz) ); temp_wall_north=0.0
     allocate(salt_wall_north(is_pe-onx:ie_pe+onx,nz) ); salt_wall_north=0.0
  endif
  if (enable_obc_south) then
     allocate(psi_wall_south(is_pe-onx:ie_pe+onx) ); psi_wall_south=0.0
     allocate(temp_wall_south(is_pe-onx:ie_pe+onx,nz) ); temp_wall_south=0.0
     allocate(salt_wall_south(is_pe-onx:ie_pe+onx,nz) ); salt_wall_south=0.0
  endif
end subroutine allocate_obc_module
