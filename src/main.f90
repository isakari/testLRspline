!> K.A.Johannessen et al./Compt. Methods Appl. Mech. Engrg 269 (2014) 471-514の2.2.4 Exampleを再現
program main
  use module_globals
!!$  use module_lrbspline
  implicit none

  !> global knot
  real(8),allocatable :: tn1(:), tn2(:)

!!$  type(LocallyRefinedSplineSurface),pointer :: lrs, iter
!!$
!!$  ! local variables
!!$  integer :: i, j
!!$  type(MeshLine) :: ml
!!$  
!!$  ncp(1)=6
!!$  ncp(2)=6
!!$  ncp(0)=(ncp(1)+1)*(ncp(2)+1)
!!$
!!$  ! tensor mesh の global knotを設定
!!$  allocate(tn1(0:ncp(1)+ndeg(1)+1))
!!$  allocate(tn2(0:ncp(2)+ndeg(2)+1))
!!$  tn1=[0.d0,0.d0,0.d0,1.d0,2.d0,4.d0,5.d0,6.d0,6.d0,6.d0]/6.d0
!!$  tn2=[0.d0,0.d0,0.d0,1.d0,2.d0,4.d0,5.d0,6.d0,6.d0,6.d0]/6.d0
!!$
!!$  ! tensor mesh の local knot を設定
!!$  call init_lrs(lrs)
!!$  do j=0,ncp(2)
!!$     do i=0,ncp(1)
!!$        call append_knot(lrs,tn1(i:i+ndeg(1)+1),tn2(j:j+ndeg(2)+1))
!!$        !write(1,*) int(tn1(i:i+ndeg(1)+1)*6), "|", int(tn2(j:j+ndeg(2)+1)*6)
!!$     end do
!!$  end do
!!$
!!$  ! LRmeshをdraw
!!$  call draw_LRmesh(lrs,"LRmesh.res")
!!$
!!$  ! 適当なsurfaceを設定
!!$  call set_greville(lrs)
!!$  iter=>lrs
!!$  do while(associated(iter))
!!$     iter%cp(3)=cos(2.d0*pi*iter%cp(1))*sin(2.d0*pi*iter%cp(2))
!!$     iter%gma=1.d0 ! tensor mesh ならば 1
!!$     iter=>iter%next
!!$  end do
!!$
!!$  ! surfaceをdraw
!!$  call draw_surface(lrs,100,"surface.res")
!!$
!!$  ! control point を draw
!!$  call set_greville(lrs)
!!$  iter=>lrs
!!$  do while(associated(iter))
!!$     ! write(11,*) iter%cp(1:3)
!!$     iter=>iter%next
!!$  end do
!!$  
!!$  ! refine (see Fig. 8)
!!$  ml%idir=2
!!$  ml%pos=3.d0/6.d0
!!$  ml%st=1.d0/6.d0
!!$  ml%en=5.d0/6.d0
!!$  call refine_lrs(lrs,ml)
!!$
!!$  iter=>lrs
!!$  do while(associated(iter))
!!$     ! write(2,*) int(iter%tn1(:)*6),"|",int(iter%tn2(:)*6)
!!$     iter=>iter%next
!!$  end do
!!$  
!!$  ! LRmeshをdraw
!!$  call draw_LRmesh(lrs,"LRmesh2.res")
!!$
!!$  ! surfaceをdraw
!!$  call draw_surface(lrs,100,"surface2.res")
!!$
!!$  ! control point を draw
!!$  call set_greville(lrs)
!!$  iter=>lrs
!!$  do while(associated(iter))
!!$     ! write(12,*) iter%cp(1:3)
!!$     iter=>iter%next
!!$  end do
!!$  
!!$  ! finalise
!!$  deallocate(tn1,tn2)
!!$  call uninit_lrs(lrs)

end program main
