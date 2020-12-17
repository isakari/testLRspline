!> K.A.Johannessen et al./Compt. Methods Appl. Mech. Engrg 269 (2014) 471-514の2.2.4 Exampleを再現
program main
  use module_globals
  use module_LRspline
  implicit none

  !> global knot
  real(8),allocatable :: tn1(:), tn2(:)

  !> LR-spline
  type(LocallyRefinedSpline),pointer :: lrs

  !> Mesh Line
  type(MeshLine) :: ml

  ! local variables
  type(LocallyRefinedSpline),pointer :: iter
  
  ! initialise
  call init_LRspline(lrs)

  ! tensor mesh = Bspline の global knotを設定
  ncp(1)=6
  ncp(2)=6
  ncp(0)=(ncp(1)+1)*(ncp(2)+1)
  allocate(tn1(0:ncp(1)+ndeg(1)+1))
  allocate(tn2(0:ncp(2)+ndeg(2)+1))
  tn1=[0.d0,0.d0,0.d0,1.d0,2.d0,4.d0,5.d0,6.d0,6.d0,6.d0]/6.d0
  tn2=[0.d0,0.d0,0.d0,1.d0,2.d0,4.d0,5.d0,6.d0,6.d0,6.d0]/6.d0

  ! tensor mesh の local knot を設定 = Bspline の knot
  call set_Bspline(lrs,tn1,tn2)
  
  ! LRmeshをdraw
  call draw_LRmesh(lrs,"LRmesh.res")

  ! 適当なsurfaceを設定
  call set_greville(lrs)
  iter=>lrs
  do while(associated(iter))
     iter%cp(3)=cos(2.d0*pi*iter%cp(1))*sin(2.d0*pi*iter%cp(2))
     iter=>iter%next
  end do

  ! surfaceをdraw
  call draw_surface(lrs,10,"surface.res")

  ! control point を draw
  open(1,file="controlpoints.res")
  iter=>lrs
  do while(associated(iter))
     write(1,*) iter%cp(1:3)
     iter=>iter%next
  end do
  close(1)
  
  ! refine (see Fig. 8)
  ml%idir=2
  ml%pos=3.d0/6.d0
  ml%st=1.d0/6.d0
  ml%en=5.d0/6.d0
  call refine_LRspline(lrs,ml)

  ! LRmeshをdraw
  call draw_LRmesh(lrs,"LRmesh2.res")

  ! surfaceをdraw
  call draw_surface(lrs,10,"surface2.res")

  ! control point を draw
  open(1,file="controlpoints2.res")
  iter=>lrs
  do while(associated(iter))
     write(1,*) iter%cp(1:3)
     iter=>iter%next
  end do
  close(1)

  ! finalise
  deallocate(tn1,tn2)
  call uninit_LRspline(lrs)

end program main
