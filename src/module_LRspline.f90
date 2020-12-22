module module_LRspline
  use module_globals
  use module_meshline
  implicit none

  private
  public ndeg, ncp, &
       LocallyRefinedSpline, &
       init_LRspline, uninit_LRspline, &
       set_Greville, &
       set_Bspline, &
       draw_LRmesh, draw_controlpoints, draw_surface
  
  integer,parameter :: ndeg(ndim)=[2,2] !< B-splineの次数
  integer :: ncp(0:ndim) !< ncp(i): x_i方向の基底関数の数-1, ncp(0)=全ての基底関数の数

  !> LR spline surface
  type LocallyRefinedSpline
     real(8) :: tn1(0:ndeg(1)+1) !< x方向のローカルノットベクトル
     real(8) :: tn2(0:ndeg(2)+1) !< y方向のローカルノットベクトル
     real(8) :: gma !< Partition of unity を満たすための係数
     real(8) :: cp(ndim+1) !< 制御変数
     type(LocallyRefinedSpline),pointer :: next
  end type LocallyRefinedSpline

contains

  !> initialise LRspline
  subroutine init_LRspline(lrs_)
    type(LocallyRefinedSpline),pointer,intent(out) :: lrs_ !< LRspline to be initialised
    nullify(lrs_)
  end subroutine init_LRspline

  !> uninitialise LRspline
  subroutine uninit_LRspline(lrs_)
    type(LocallyRefinedSpline),pointer,intent(inout) :: lrs_ !< LRspline to be destroyed

    type(LocallyRefinedSpline),pointer :: lrs, tmp

    if(.not.associated(lrs_))then !lrs_がなければ
       return !何もしない
    end if

    lrs=>lrs_ !=lrs_の先頭へのpointer
    do while(associated(lrs))
       tmp=>lrs !=lrsの現在の先頭へのpointer
       lrs=>lrs%next !=lrsの現在の先頭の次
       deallocate(tmp) !lrsの現在の先頭をdeallocate
    end do

    nullify(lrs_)
    
  end subroutine uninit_LRspline

  !> return the last LocallyRefinedSplineSurface
  function tail(a) result(b)
    type(LocallyRefinedSpline),pointer,intent(in) :: a
    type(LocallyRefinedSpline),pointer :: b
    
    b=>a !bでaの先頭を指し
    do while(associated(b%next)) !bの次がnullになるまで
       b=>b%next !bを辿っていく
    end do

  end function tail

  !> new_lrsのknot(とnext)をset
  subroutine set_new_knot(new_lrs,tn1,tn2,next)
    type(LocallyRefinedSpline),pointer,intent(inout) :: new_lrs
    real(8),intent(in) :: tn1(0:ndeg(1)+1), tn2(0:ndeg(2)+1)
    type(LocallyRefinedSpline),pointer,optional,intent(in) :: next

    nullify(new_lrs);
    allocate(new_lrs) !空の領域をallocate

    if(present(next))then !nextが指定されていれば
       new_lrs%next=>next !new_lrsの次にnextを渡す
    else !そうでなければ
       nullify(new_lrs%next) !new_lrsの次は何もない
    end if

    !knotをnew_lrsに渡す
    new_lrs%tn1=tn1
    new_lrs%tn2=tn2

  end subroutine set_new_knot
  
  !> append an lrs function at the last 
  subroutine append_knot(lrs_,tn1,tn2)
    type(LocallyRefinedSpline),pointer,intent(inout) :: lrs_ !< LRspline
    real(8),intent(in) :: tn1(0:ndeg(1)+1) !< x方向のglobal knot
    real(8),intent(in) :: tn2(0:ndeg(2)+1) !< y方向のglobal knot

    type(LocallyRefinedSpline),pointer :: lrs

    if(associated(lrs_))then !lrs_があれば
       lrs=>tail(lrs_) !aのお尻
       call set_new_knot(lrs%next,tn1,tn2) !の次にtn1,tn2をset
    else !aがなければ
       call set_new_knot(lrs_,tn1,tn2) !a(の先頭)にtn1,tn2をset
    end if

  end subroutine append_knot
  
  !> global knot (tn1|tn2) に対応する Bspline surfaceをset
  subroutine set_Bspline(lrs_,tn1,tn2)
    type(LocallyRefinedSpline),pointer,intent(inout) :: lrs_ !< LRspline
    real(8),intent(in) :: tn1(0:ncp(1)+ndeg(1)+1) !< x方向のglobal knot
    real(8),intent(in) :: tn2(0:ncp(2)+ndeg(2)+1) !< y方向のglobal knot

    integer :: i, j
    type(LocallyRefinedSpline),pointer :: lrs

    do j=0,ncp(2)
       do i=0,ncp(1)
          ! local knot を lrs_ の最後にくっつけていく
          call append_knot(lrs_,tn1(i:i+ndeg(1)+1),tn2(j:j+ndeg(2)+1))
       end do
    end do

    lrs=>lrs_
    do while(associated(lrs))
       lrs%gma=1.d0 ! Bスプラインならば必ず1
       lrs=>lrs%next
    end do
    
  end subroutine set_Bspline

  !> cp(1:2)にGreville座標をset
  subroutine set_Greville(lrs_)
    type(LocallyRefinedSpline),pointer,intent(inout) :: lrs_
    type(LocallyRefinedSpline),pointer :: lrs

    lrs=>lrs_
    do while(associated(lrs))
       lrs%cp(1)=sum(lrs%tn1(1:ndeg(1)))/dble(ndeg(1))
       lrs%cp(2)=sum(lrs%tn2(1:ndeg(2)))/dble(ndeg(2))
       lrs=>lrs%next
    end do

  end subroutine set_Greville

  !> Bスプライン基底関数 (local knot版) B^np[tn(0),tn(1),...,tn(np+1)](t)を計算\\
  !> open knot を仮定
  !> 中間のknotは重ならないことを仮定
  recursive function basis_func(np,tn,t) result(ret)
    integer,intent(in) :: np !< Bスプラインの次数
    real(8),intent(in) :: tn(0:np+1) !< local knot vec
    real(8),intent(in) :: t  !< argument of Bspline
    real(8) :: ret !< result (Bspline value)
    
    integer :: i
    logical :: chk
    real(8) :: w1, w2
    
    ! quick return if possible
    if(np.gt.0)then
       chk=true  ! 左端
       do i=0,np-1
           !local knotの左p+1個が重なっているかcheck
          chk=chk.and.(abs(tn(i+1)-tn(i)).le.tiny)
       end do
       !local knotの左p+1個が重なっていてかつ左端の点の基底を計算する
       if(chk.and.(abs(t).le.tiny))then
          ret=1.d0
          return
       end if

       chk=true ! 右端
       do i=1,np
          !local knotの右p+1個が重なっているかcheck
          chk=chk.and.(abs(tn(i+1)-tn(i)).le.tiny)
       end do
       !local knotの右p+1個が重なっていてかつ右端の点の基底を計算する
       if(chk.and.(abs(t-1.d0).le.tiny))then 
          ret=1.d0
          return
       end if
    end if
       
    if(t.lt.tn(0).or.t.gt.tn(np+1))then
       ret=0.d0
       return
    end if
       
    ! Cox-de Boor
    w1=0.d0
    w2=0.d0
    if(np.eq.0)then
       if(t.ge.tn(0).and.t.lt.tn(1))then
          ret=1.d0
          return
       else
          ret=0.d0
          return
       end if
    else
       if(abs(tn(np)-tn(0)).gt.tiny)then
          w1=(t-tn(0))/(tn(np)-tn(0))*basis_func(np-1,tn(0:np),t)
       end if
       if(abs(tn(np+1)-tn(1)).gt.tiny)then
          w2=(tn(np+1)-t)/(tn(np+1)-tn(1))*basis_func(np-1,tn(1:np+1),t)
       end if
       ret=w1+w2
       return
    end if
    
  end function basis_func
  
  !> draw LRmesh
  !> \brief gnuplotでp file w lp
  subroutine draw_LRmesh(lrs_,fn)
    type(LocallyRefinedSpline),pointer,intent(in) :: lrs_
    character(*),intent(in) :: fn

    type(LocallyRefinedSpline),pointer :: lrs
    integer :: i, j
    
    lrs=>lrs_
    open(1,file=fn)
    do while(associated(lrs))
       do j=0,ndeg(2)
          do i=0,ndeg(1)
             write(1,*) lrs%tn1(i  ),lrs%tn2(j  )
             write(1,*) lrs%tn1(i+1),lrs%tn2(j  )
             write(1,*) lrs%tn1(i+1),lrs%tn2(j+1)
             write(1,*) lrs%tn1(i  ),lrs%tn2(j+1)
             write(1,*) lrs%tn1(i  ),lrs%tn2(j  )
             write(1,*)
             write(1,*)
          end do
       end do
       lrs=>lrs%next
    end do
    close(1)

  end subroutine draw_LRmesh
    
  !> draw control points
  subroutine draw_controlpoints(lrs_,fn)
    type(LocallyRefinedSpline),pointer,intent(in) :: lrs_
    character(*),intent(in) :: fn

    type(LocallyRefinedSpline),pointer :: lrs
    
    open(1,file=fn)
    lrs=>lrs_
    do while(associated(lrs))
       write(1,*) lrs%cp
       lrs=>lrs%next
    end do
    close(1)
    
  end subroutine draw_controlpoints
    
  !> draw surface
  !> \brief gnuplotでp file
  !> \param n (n+1)x(n+1)のgridでsurfaceを描画
  subroutine draw_surface(lrs_,n,fn)
    type(LocallyRefinedSpline),pointer,intent(in) :: lrs_
    integer,intent(in) :: n
    character(*),intent(in) :: fn

    type(LocallyRefinedSpline),pointer :: lrs
    integer :: i, j
    real(8) :: t1, t2, b1, b2, c(ndim+1)
    
    open(1,file=fn)
    do j=0,n
       t2=dble(j)/dble(n)
       do i=0,n
          t1=dble(i)/dble(n)
          c(:)=0.d0
          lrs=>lrs_
          do while(associated(lrs))
             b1=basis_func(ndeg(1),lrs%tn1,t1)
             if(abs(b1).le.tiny) goto 1
             b2=basis_func(ndeg(2),lrs%tn2,t2)
             if(abs(b2).le.tiny) goto 1
             c(:)=c(:)+b1*b2*lrs%gma*lrs%cp(:)
1            lrs=>lrs%next
          end do
          write(1,*) c(:)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine draw_surface
  
end module module_LRspline
