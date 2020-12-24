!> LR-splineによる曲面近似
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

  !> surface data to be approximated
  integer :: ngrid(-1:ndim)
  integer,allocatable :: igrid(:)
  real(8),allocatable :: ff(:,:)!< データ
  real(8),allocatable :: f(:,:)!< データ
  real(8),allocatable :: df(:,:,:) !< データの勾配
  
  ! local variables
  type(LocallyRefinedSpline),pointer :: lrs
  type(MeshLine),pointer :: ml
  integer :: i, j, k, icnt, keast, kwest, ksouth, knorth, k0, info
  real(8) :: pos, t1, t2
  logical :: diag
  integer,allocatable :: ipiv(:)
  real(8),allocatable :: bmat(:,:), d1bmat(:,:), d2bmat(:,:), cv(:,:), amat(:,:)

  ! dataを読み込み
  open(1,file="examples/diag_waterfall.dat")
  read(1,*) ngrid(1)
  read(1,*) ngrid(2)
  ngrid(0)=(ngrid(1)+1)*(ngrid(2)+1) !設計領域内のgird点の数
  ngrid(-1)=(ngrid(1)+3)*(ngrid(2)+3) !設計領域+バッファ領域のgrid点の数
  allocate(ff(ndim+1,ngrid(-1)))
  k=1
  do j=-1,ngrid(2)+1
     do i=-1,ngrid(1)+1
        read(1,*) ff(:,k)
        k=k+1
     end do
  end do
  close(1)

  ! dataの勾配を計算
  allocate(igrid(ngrid(0)),df(2,ndim+1,ngrid(0)))
  icnt=1
  do j=-1,ngrid(2)+1
     do i=-1,ngrid(1)+1
        k=i+(j+1)*(ngrid(1)+3)+2
        if(i.ge.0.and.i.le.ngrid(1).and.j.ge.0.and.j.le.ngrid(2)) then !設計領域内部なら
           igrid(icnt)=k
           icnt=icnt+1
        end if
     end do
  end do
  allocate(f(ndim+1,ngrid(0)))
  do j=0,ngrid(2)
     do i=0,ngrid(1)
        k=i+j*(ngrid(1)+1)+1
        f(:,k)=ff(:,igrid(k))
        ! write(11,*) f(:,k)
     end do
  end do

  do j=0,ngrid(2)
     do i=0,ngrid(1)
        k=i+(j+1)*(ngrid(1)+3)+2
        k0=i+j*(ngrid(1)+1)+1
        keast=i+1+(j+1)*(ngrid(1)+3)+2
        kwest=i-1+(j+1)*(ngrid(1)+3)+2
        ksouth=i+(j+1-1)*(ngrid(1)+3)+2
        knorth=i+(j+1+1)*(ngrid(1)+3)+2
        df(1,:,k0)=(ff(:,keast)-ff(:,kwest))/(2.d0/dble(ngrid(1))) !中央
        df(2,:,k0)=(ff(:,knorth)-ff(:,ksouth))/(2.d0/dble(ngrid(2))) !差分
        if(abs(df(1,3,k0)).gt.10.d0) df(1,3,k0)=0.d0 !不連続点の
        if(abs(df(2,3,k0)).gt.10.d0) df(2,3,k0)=0.d0 !勾配はゼロにしておく
        ! write(12,*) f(:,k0),df(1,:,k0),df(2,:,k0)
     end do
  end do
  
  ! initialise
  call init_LRspline(lrs_)
  call init_meshline(ml_)
  
  ! tensor mesh = Bspline の global knotを設定
  ncp(1)=5
  ncp(2)=5
  ncp(0)=(ncp(1)+1)*(ncp(2)+1)
  allocate(tn1(0:ncp(1)+ndeg(1)+1))
  allocate(tn2(0:ncp(2)+ndeg(2)+1))
  call set_open_uniform_global_knot(ncp,tn1,tn2)

  ! tensor mesh の local knot を設定 = Bspline の knot
  call set_Bspline(lrs_,tn1,tn2)

  ! initialise the control points
  lrs=>lrs_
  do while(associated(lrs))
     lrs%cp(:)=0.d0
     lrs=>lrs%next
  end do
  
  do j=1,4 !level 4 まで
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

  ncp(0)=0
  lrs=>lrs_
  do while(associated(lrs))
     ncp(0)=ncp(0)+1
     lrs=>lrs%next
  end do
  
  ! H1近似
  allocate(bmat(ngrid(0),ncp(0)))
  allocate(d1bmat(ngrid(0),ncp(0)))
  allocate(d2bmat(ngrid(0),ncp(0)))
  allocate(cv(ndim+1,ncp(0)))
  allocate(amat(ncp(0),ncp(0)))
  allocate(ipiv(ncp(0)))
  amat(:,:)=0.d0
  bmat(:,:)=0.d0
  cv(:,:)=0.d0
  icnt=1; lrs=>lrs_
  do while(associated(lrs))
     do j=0,ngrid(2)
        do i=0,ngrid(1)
           k=i+j*(ngrid(1)+1)+1
           t1=dble(i)/dble(ngrid(1))
           t2=dble(j)/dble(ngrid(2))
           bmat(k,icnt)=lrs%gma*basis_func(ndeg(1),lrs%tn1,t1)*basis_func(ndeg(2),lrs%tn2,t2)
           d1bmat(k,icnt)=lrs%gma*diff_basis_func(ndeg(1),lrs%tn1,t1)*basis_func(ndeg(2),lrs%tn2,t2)
           d2bmat(k,icnt)=lrs%gma*basis_func(ndeg(1),lrs%tn1,t1)*diff_basis_func(ndeg(2),lrs%tn2,t2)
        end do
     end do
     icnt=icnt+1; lrs=>lrs%next
  end do
  amat=matmul(transpose(bmat),bmat)&
       +matmul(transpose(d1bmat),d1bmat)/dble(ncp(1)-ndeg(1)+1)**2*0.02d0&
       +matmul(transpose(d2bmat),d2bmat)/dble(ncp(2)-ndeg(2)+1)**2*0.02d0
  call DGETRF(ncp(0),ncp(0),amat,ncp(0),ipiv,info)
  cv(1,:)=matmul(transpose(bmat),f(1,:))&
       +matmul(transpose(d1bmat),df(1,1,:))/dble(ncp(1)-ndeg(1)+1)**2*0.02d0&
       +matmul(transpose(d2bmat),df(2,1,:))/dble(ncp(2)-ndeg(2)+1)**2*0.02d0
  call DGETRS("N",ncp(0),1,amat,ncp(0),ipiv,cv(1,:),ncp(0),info)
  cv(2,:)=matmul(transpose(bmat),f(2,:))&
       +matmul(transpose(d1bmat),df(1,2,:))/dble(ncp(1)-ndeg(1)+1)**2*0.02d0&
       +matmul(transpose(d2bmat),df(2,2,:))/dble(ncp(2)-ndeg(2)+1)**2*0.02d0
  call DGETRS("N",ncp(0),1,amat,ncp(0),ipiv,cv(2,:),ncp(0),info)
  cv(3,:)=matmul(transpose(bmat),f(3,:))&
       +matmul(transpose(d1bmat),df(1,3,:))/dble(ncp(1)-ndeg(1)+1)**2*0.02d0&
       +matmul(transpose(d2bmat),df(2,3,:))/dble(ncp(2)-ndeg(2)+1)**2*0.02d0
  call DGETRS("N",ncp(0),1,amat,ncp(0),ipiv,cv(3,:),ncp(0),info)

  ! 制御変数
  icnt=1; lrs=>lrs_
  do while(associated(lrs))
     lrs%cp(:)=cv(:,icnt)
     icnt=icnt+1; lrs=>lrs%next
  end do
  
  ! draw refined LR-spline
  call draw_LRmesh(lrs_,"LRmesh.res")
  call draw_controlpoints(lrs_,"controlpoints.res")
  call draw_surface(lrs_,50,"surface.res")
  
  ! finalise
  deallocate(tn1,tn2)
  call uninit_meshline(ml_)
  call uninit_LRspline(lrs_)
  deallocate(bmat,d1bmat,d2bmat,cv,amat,ipiv)
  deallocate(igrid,ff,f,df)
  
end program main
