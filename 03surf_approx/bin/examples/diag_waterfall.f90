!> 対角の滝
program diag_waterfall
  implicit none
  
  integer(8),parameter :: n=50
  real(8),parameter :: xmin=-2.d0
  real(8),parameter :: xmax=+2.d0
  real(8),parameter :: ymin=-2.d0
  real(8),parameter :: ymax=+2.d0
  
  integer :: i, j
  real(8) :: x, y, z
  
  open(1,file="diag_waterfall.dat")
  write(1,*) n ! x方向のデータの個数-1
  write(1,*) n ! y方向のデータの個数-1
  do j=-1,n+1
     y=ymin+(ymax-ymin)*dble(j)/dble(n)
     do i=-1,n+1
        x=xmin+(xmax-xmin)*dble(i)/dble(n)
        z=-1.d0
        if(y.gt.x) then
           z=1.d0
        end if
        write(1,*) x, y, z
     end do
  end do
  close(1)
end program diag_waterfall
