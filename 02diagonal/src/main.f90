!> K.A.Johannessen et al./Compt. Methods Appl. Mech. Engrg 269 (2014) 471-514 の diagonal structured refinement を再現
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
  integer :: i, j, k
  real(8) :: pos
  logical :: diag
  
  ! initialise
  call init_LRspline(lrs_)
  call init_meshline(ml_)
  
  ! tensor mesh = Bspline の global knotを設定
  ncp(1)=5
  ncp(2)=5
  ncp(0)=(ncp(1)+1)*(ncp(2)+1)
  allocate(tn1(0:ncp(1)+ndeg(1)+1))
  allocate(tn2(0:ncp(2)+ndeg(2)+1))
  tn1=[0.d0,0.d0,0.d0,1.d0,2.d0,3.d0,4.d0,4.d0,4.d0]/4.d0
  tn2=[0.d0,0.d0,0.d0,1.d0,2.d0,3.d0,4.d0,4.d0,4.d0]/4.d0

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

  do j=1,7 !level 8 まで
     ! set meshlines
     lrs=>lrs_
     do while(associated(lrs))
        diag=true
        do i=0,ndeg(1)+1
           diag=diag.and.(abs(lrs%tn1(i)-lrs%tn2(i)).le.tiny)
        end do
        if(diag) then
           do k=0,ndeg(1)
              pos=(lrs%tn1(k)+lrs%tn1(k+1))*0.5d0
              if((pos.gt.tiny).and.(pos.lt.1.d0-tiny))then
                 call append_meshline(ml_,2,pos,lrs%tn2(0),lrs%tn2(ndeg(2)+1))
              end if
           end do
           do k=0,ndeg(2)
              pos=(lrs%tn2(k)+lrs%tn2(k+1))*0.5d0
              if((pos.gt.tiny).and.(pos.lt.1.d0-tiny))then
                 call append_meshline(ml_,1,pos,lrs%tn1(0),lrs%tn1(ndeg(1)+1))
              end if
           end do
        end if
        lrs=>lrs%next
     end do

     ! draw meshline
     call draw_meshline(ml_,"meshline.res")

     ! refine LR-splines
     ml=>ml_
     do while(associated(ml))
        call refine_LRspline(lrs_,ml_,ml)
        ml=>ml%next
     end do
  end do
     
  ! draw refined LR-spline
  call draw_LRmesh(lrs_,"LRmesh2.res")
  call draw_controlpoints(lrs_,"controlpoints2.res")
  call draw_surface(lrs_,50,"surface2.res")

  ! finalise
  deallocate(tn1,tn2)
  call uninit_meshline(ml_)
  call uninit_LRspline(lrs_)

end program main
