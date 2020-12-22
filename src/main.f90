!> K.A.Johannessen et al./Compt. Methods Appl. Mech. Engrg 269 (2014) 471-514の2.2.4 Exampleを再現
program main
  use module_globals
  use module_meshline
  use module_LRspline
  implicit none

  !> global tensor knot
  real(8),allocatable :: tn1(:), tn2(:)

  !> LR-spline
  type(LocallyRefinedSpline),pointer :: lrs_

  !> Mesh Line
  type(MeshLine),pointer :: ml_

  ! local variables
  type(LocallyRefinedSpline),pointer :: lrs
  type(MeshLine),pointer :: ml

  ! initialise
  call init_LRspline(lrs_)
  call init_meshline(ml_)
  
  ! tensor mesh = Bspline の global knotを設定
  ncp(1)=6
  ncp(2)=6
  ncp(0)=(ncp(1)+1)*(ncp(2)+1)
  allocate(tn1(0:ncp(1)+ndeg(1)+1))
  allocate(tn2(0:ncp(2)+ndeg(2)+1))
  tn1=[0.d0,0.d0,0.d0,1.d0,2.d0,4.d0,5.d0,6.d0,6.d0,6.d0]/6.d0
  tn2=[0.d0,0.d0,0.d0,1.d0,2.d0,4.d0,5.d0,6.d0,6.d0,6.d0]/6.d0

  ! tensor mesh の local knot を設定 = Bspline の knot
  call set_Bspline(lrs_,tn1,tn2)

  ! 適当なsurfaceを設定
  call set_Greville(lrs_)
  lrs=>lrs_
  do while(associated(lrs))
     lrs%cp(3)=cos(2.d0*pi*lrs%cp(1))*sin(2.d0*pi*lrs%cp(2))
     lrs=>lrs%next
  end do

  ! draw itinial Bspline
  call draw_LRmesh(lrs_,"LRmesh.res")
  call draw_controlpoints(lrs_,"controlpoints.res")
  call draw_surface(lrs_,50,"surface.res")

  ! insert meshline
  call append_meshline(ml_,2,3.0d0/6.d0,0.d0/6.d0,1.d0/6.d0)
  call append_meshline(ml_,2,3.0d0/6.d0,1.d0/6.d0,5.d0/6.d0)
  call append_meshline(ml_,2,3.0d0/6.d0,5.d0/6.d0,6.d0/6.d0)
  call append_meshline(ml_,1,3.0d0/6.d0,5.d0/6.d0,6.d0/6.d0)
  call append_meshline(ml_,1,3.0d0/6.d0,0.d0/6.d0,1.d0/6.d0)
  call append_meshline(ml_,1,3.0d0/6.d0,1.d0/6.d0,5.d0/6.d0)
  call append_meshline(ml_,1,3.0d0/6.d0,1.d0/6.d0,5.d0/6.d0)
  call append_meshline(ml_,2,2.9d0/6.d0,1.d0/6.d0,5.d0/6.d0)
  call append_meshline(ml_,2,2.8d0/6.d0,1.d0/6.d0,5.d0/6.d0)
  call append_meshline(ml_,2,2.8d0/6.d0,0.d0/6.d0,1.d0/6.d0)
  call append_meshline(ml_,1,1.5d0/6.d0,1.d0/6.d0,2.9/6.d0)
  call append_meshline(ml_,1,1.5d0/6.d0,5.d0/6.d0,6.d0/6.d0)
  call append_meshline(ml_,1,1.5d0/6.d0,2.9/6.d0,5.d0/6.d0)
  
  open(1,file="meshline.res")
  ml=>ml_
  do while(associated(ml))
     if(ml%idir.eq.1)then !horizontal meshline
        write(1,*) ml%st,ml%pos
        write(1,*) ml%en,ml%pos
        write(1,*); write(1,*)
     else !vertical meshline
        write(1,*) ml%pos,ml%st
        write(1,*) ml%pos,ml%en
        write(1,*); write(1,*)
     end if
     ml=>ml%next
  end do
  close(1)
    
  ! finalise
  deallocate(tn1,tn2)
  call uninit_meshline(ml_)
  call uninit_LRspline(lrs_)

end program main
