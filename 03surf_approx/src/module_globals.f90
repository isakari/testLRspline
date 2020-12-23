module module_globals
  implicit none

  integer,parameter :: ndim=2 !< 次元
  real(8),parameter :: pi=acos(-1.d0) !< pi
  real(8),parameter :: tiny=epsilon(1.d0) !< 小さい数

  logical,parameter :: true=.true.
  logical,parameter :: false=.false.
  
end module module_globals
