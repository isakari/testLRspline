module module_LRspline
  use module_globals
  implicit none

  private
  public ndeg, ncp, &
       LocallyRefinedSpline, &
       MeshLine, &
       init_LRspline, uninit_LRspline, &
       set_Bspline, &
       set_greville, &
       refine_LRspline, &
       draw_LRmesh, draw_surface
  
  integer,parameter :: ndeg(ndim)=[2,2] !< B-splineの次数
  integer :: ncp(0:ndim) !< 基底関数の数-1

  !> LR spline surface
  type LocallyRefinedSpline
     real(8) :: tn1(0:ndeg(1)+1) !< x方向のローカルノットベクトル
     real(8) :: tn2(0:ndeg(2)+1) !< y方向のローカルノットベクトル
     real(8) :: gma !< Partition of unity を満たすための係数
     real(8) :: cp(ndim+1) !< 制御変数
     type(LocallyRefinedSpline),pointer :: next
  end type LocallyRefinedSpline

  !> mesh line
  !> \brief if idir=1, meshline is (st,pos)-->(en,pos)
  !> \brief if idir=2, meshline is (pos,st)-->(pos,en)
  type MeshLine
     integer :: idir !< meshlineの方向 1(=x) or 2(=y)
     real(8) :: pos !< meshlineの位置
     real(8) :: st !< meshlineの始点
     real(8) :: en !< meshlineの終点
  end type MeshLine

contains

  !> initialise LRspline
  subroutine init_LRspline(a)
    type(LocallyRefinedSpline),pointer,intent(out) :: a !< LRspline to be initialised
    nullify(a)
  end subroutine init_LRspline

  !> uninitialise the list 
  subroutine uninit_LRspline(a)
    type(LocallyRefinedSpline),pointer,intent(inout) :: a !< LRspline to be destroyed

    type(LocallyRefinedSpline),pointer :: iter, tmp

    if(.not.associated(a))then
       return
    end if

    iter=>a
    do while(associated(iter))
       tmp=>iter !=aの現在の先頭へのpointer
       iter=>iter%next !=aの現在の先頭の次
       deallocate(tmp) !aの現在の先頭をdeallocate
    end do

    nullify(a)
    
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
  subroutine set_knot(new_lrs,tn1,tn2,next)
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

  end subroutine set_knot
  
  !> append an lrs function at the last 
  subroutine append_knot(a,tn1,tn2)
    type(LocallyRefinedSpline),pointer,intent(inout) :: a !< LRspline
    real(8),intent(in) :: tn1(0:ndeg(1)+1) !< x方向のglobal knot
    real(8),intent(in) :: tn2(0:ndeg(2)+1) !< y方向のglobal knot

    type(LocallyRefinedSpline),pointer :: iter

    if(associated(a))then !aがあれば
       iter=>tail(a) !aのお尻
       call set_knot(iter%next,tn1,tn2) !の次にtn1,tn2をset
    else !aがなければ
       call set_knot(a,tn1,tn2) !a(の先頭)にtn1,tn2をset
    end if

  end subroutine append_knot

  
  !> global knot (tn1|tn2) に対応する Bspline surfaceをset
  subroutine set_Bspline(a,tn1,tn2)
    type(LocallyRefinedSpline),pointer,intent(inout) :: a !< LRspline
    real(8),intent(in) :: tn1(0:ncp(1)+ndeg(1)+1) !< x方向のglobal knot
    real(8),intent(in) :: tn2(0:ncp(2)+ndeg(2)+1) !< y方向のglobal knot

    integer :: i, j
    type(LocallyRefinedSpline),pointer :: iter

    do j=0,ncp(2)
       do i=0,ncp(1)
          ! local knot を a の最後にくっつけていく
          call append_knot(a,tn1(i:i+ndeg(1)+1),tn2(j:j+ndeg(2)+1))
       end do
    end do

    iter=>a
    do while(associated(iter))
       iter%gma=1.d0 ! Bスプラインならば必ず1
       iter=>iter%next
    end do
    
  end subroutine set_Bspline

  !> cp(1:2)にGreville座標をset
  subroutine set_greville(a)
    type(LocallyRefinedSpline),pointer,intent(inout) :: a
    type(LocallyRefinedSpline),pointer :: iter

    iter=>a
    do while(associated(iter))
       iter%cp(1)=sum(iter%tn1(1:ndeg(1)))/dble(ndeg(1))
       iter%cp(2)=sum(iter%tn2(1:ndeg(2)))/dble(ndeg(2))
       iter=>iter%next
    end do

  end subroutine set_greville

  !> quick sort
  pure recursive function qsort(x) result(res)
    real(8),intent(in) :: x(:)
    real(8) :: res(size(x))
    integer :: n
    n=size(x)
    if(n>1)then
       res=[qsort(pack((x(2:)),x(2:)<x(1))),& !x(2:)のなかでx(1)より小さいものをソート
            x(1),& !x(1)を並べる
            qsort(pack(x(2:),x(2:)>=x(1)))]!x(2:)のなかでx(1)以上のものをソート
    else
       res=x
    end if
  end function qsort

  !> insert a LR-spline function at the given index
  subroutine insert_LRspline(a,idx,tn1,tn2,cp,gma)
    type(LocallyRefinedSpline),pointer,intent(inout) :: a !< LRspline
    integer,intent(in) :: idx !< aの(idx)のうしろに lrs function をはさみ、tn1&tn2を渡す
    real(8),intent(in) :: tn1(0:ndeg(1)+1) !< 挿入するローカルノットベクトル(xi)
    real(8),intent(in) :: tn2(0:ndeg(2)+1) !< 挿入するローカルノットベクトル(eta)
    real(8),intent(in) :: cp(ndim+1) !< 挿入するローカルノットベクトルに対応する制御変数
    real(8),intent(in) :: gma !< 挿入するローカルノットベクトルに対応する重み
    
    integer :: i
    type(LocallyRefinedSpline),pointer :: iter, prev, new_lrs

    ! find a node after which a new node is inserted
    nullify(prev)
    iter=>a
    do i=0,idx
       !                   prev            iter
       !(idx-1)---idx-1--->(idx)---idx--->(idx+1)
       prev=>iter
       iter=>iter%next
    end do

    ! create a new LR spline function
    nullify(new_lrs)
    !new_lrs-->iterとつなぐ
    call set_knot(new_lrs,tn1,tn2,next=iter)
    new_lrs%cp=cp
    new_lrs%gma=gma
    
    ! prev-->new_lrs とつなぐ
    prev%next=>new_lrs

  end subroutine insert_LRspline

  !> remove a LR-spline function at the given index
  subroutine remove_LRspline(a,idx)
    type(LocallyRefinedSpline),pointer,intent(inout) :: a
    integer,intent(in) :: idx

    integer :: i
    type(LOcallyRefinedSpline), pointer :: iter, prev

    if(idx<0) then
       stop "aho"
    end if

    ! find a node to be removed
    nullify(prev)
    iter=>a
    do i=0,idx-1
       ! prev            iter
       !(idx-1)--idx-1-->(idx)
       prev=>iter
       iter=>iter%next
    end do

    prev%next=>iter%next
    deallocate(iter)

  end subroutine remove_LRspline
  
  !> Algorithm 1 (local xi split)
  !> \brief LR-spline aの中の LR-spline function b を xi=t により refine する
  subroutine local_split_1(a,b,idx,t)
    type(LocallyRefinedSpline),pointer,intent(inout) :: a !< LR-spline
    type(LocallyRefinedSpline),intent(in) :: b !< LR-spline function to be refined
    integer,intent(inout) :: idx !< b is idx-th function in a
    real(8),intent(in) :: t !< inserted knot value

    real(8) :: alp(2)
    real(8) :: tn(0:ndeg(1)+2) !xi
    real(8) :: tn1a(0:ndeg(1)+1) !xi1
    real(8) :: tn1b(0:ndeg(1)+1) !xi2
    real(8) :: tn2a(0:ndeg(2)+1) !eta1
    real(8) :: tn2b(0:ndeg(2)+1) !eta2
    real(8) :: dis1(0:ndeg(1)+1), dis2(0:ndeg(2)+1)
    logical :: exist
    type(LocallyRefinedSpline),pointer :: iter
    
    ! alpha
    alp(1)=1.d0
    if(t.le.b%tn1(ndeg(1))) alp(1)=(t-b%tn1(0))/(b%tn1(ndeg(1))-b%tn1(0))
    alp(2)=1.d0
    if(b%tn1(1).le.t) alp(2)=(b%tn1(ndeg(1)+1)-t)/(b%tn1(ndeg(1)+1)-b%tn1(1))
    ! write(*,*) alp(:)
    
    ! xi, xi1, xi2
    tn(0)=t
    tn(1:ndeg(1)+2)=b%tn1(:)
    tn=qsort(tn)
    ! write(*,*) tn
    tn1a(:)=tn(0:ndeg(1)+1)
    tn1b(:)=tn(1:ndeg(1)+2)
    ! write(*,*) tn1a
    ! write(*,*) tn1b

    ! eta1, eta2
    tn2a(:)=b%tn2(:)
    tn2b(:)=b%tn2(:)
    ! write(*,*) tn2a
    ! write(*,*) tn2b

    ! 8-15行目
    exist=false
    iter=>a
    do while(associated(iter))
       dis1=abs(iter%tn1-tn1a)
       dis2=abs(iter%tn2-tn2a)
       if((maxval(dis1).le.tiny).and.(maxval(dis2).le.tiny))then
          exist=true
          iter%cp(:)=(iter%cp(:)*iter%gma+b%cp(:)*b%gma*alp(1))/(iter%gma+alp(1)*b%gma)
          iter%gma=iter%gma+alp(1)*b%gma
          exit
       end if
       iter=>iter%next
    end do
    if(exist.eqv.false)then ! idxの前に local knot を追加
       call insert_LRspline(a,idx-1,tn1a,tn2a,b%cp,alp(1)*b%gma)
       idx=idx+1 !追加した分
    end if

    ! 16-23行目
    exist=false
    iter=>a
    do while(associated(iter))
       dis1=abs(iter%tn1-tn1b)
       dis2=abs(iter%tn2-tn2b)
       if((maxval(dis1).le.tiny).and.(maxval(dis2).le.tiny))then
          exist=true
          iter%cp(:)=(iter%cp(:)*iter%gma+b%cp(:)*b%gma*alp(2))/(iter%gma+alp(2)*b%gma)
          iter%gma=iter%gma+alp(2)*b%gma
          exit
       end if
       iter=>iter%next
    end do
    if(exist.eqv.false)then ! idxの前に local knot を追加
       call insert_LRspline(a,idx-1,tn1b,tn2b,b%cp,alp(2)*b%gma)
       idx=idx+1 !追加した分
    end if

    ! 24行目
    call remove_LRspline(a,idx)
    idx=idx-1 !削除した分
    
  end subroutine local_split_1
  
  !> Algorithm 2
  subroutine refine_LRspline(a,ml)
    type(LocallyRefinedSpline),pointer,intent(inout) :: a !< LR-spline
    type(MeshLine),intent(in) :: ml !< mesh line

    type(LocallyRefinedSpline),pointer :: iter
    type(LocallyRefinedSpline) :: b
    integer :: icnt

    if(.not.(ml%idir.eq.1.or.ml%idir.eq.2)) stop "Error in ml%idir @ refine_lrs"

    ! step1
    icnt=0
    iter=>a
    do while(associated(iter))
       if(ml%idir.eq.2)then ! vertical meshline insertion
          if((iter%tn1(0)<ml%pos)& !mlがiterのsupportを分断するなら
               .and.(ml%pos<iter%tn1(ndeg(1)+1))&
               .and.(ml%st.le.iter%tn2(0))&
               .and.(iter%tn2(ndeg(2)+1)).le.ml%en)then
             b=iter
             nullify(b%next)
             call local_split_1(a,b,icnt,ml%pos)
          end if
       else ! horizontal meshline insertion
!!$          if((iter%tn2(0)<ml%pos)& !mlがiterのsupportを分断するなら
!!$               .and.(ml%pos<iter%tn2(ndeg(2)+1)).and.&
!!$               (ml%st.le.iter%tn1(0))&
!!$               .and.(iter%tn1(ndeg(1)+1)).le.ml%en)then 
!!$             call local_split_2(a,iter,ml%pos,icnt) !local eta split
!!$          end if
       end if
       iter=>iter%next
       icnt=icnt+1
    end do

!!$    !test
!!$    icnt=0
!!$    iter=>a
!!$    do while(associated(iter))
!!$       write(*,*) int(iter%tn1*6),"|",int(iter%tn2*6)
!!$       icnt=icnt+1
!!$       iter=>iter%next
!!$    end do
    
  end subroutine refine_LRspline
  

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
       chk=.true.  ! 左端
       do i=0,np-1
           !local knotの左p+1個が重なっているかcheck
          chk=and(chk,abs(tn(i+1)-tn(i)).le.tiny)
       end do
       !local knotの左p+1個が重なっていてかつ左端の点の基底を計算する
       if((chk.eqv..true.).and.(abs(t).le.tiny))then
          ret=1.d0
          return
       end if

       chk=.true. ! 右端
       do i=1,np
          !local knotの右p+1個が重なっているかcheck
          chk=and(chk,abs(tn(i+1)-tn(i)).le.tiny)
       end do
       !local knotの右p+1個が重なっていてかつ右端の点の基底を計算する
       if((chk.eqv..true.).and.(abs(t-1.d0).le.tiny))then 
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
  subroutine draw_LRmesh(a,fn)
    type(LocallyRefinedSpline),pointer,intent(in) :: a
    character(*),intent(in) :: fn

    type(LocallyRefinedSpline),pointer :: iter    
    integer :: i, j
    
    iter=>a
    open(1,file=fn)
    do while(associated(iter))
       do j=0,ndeg(2)
          do i=0,ndeg(1)
             write(1,*) iter%tn1(i  ),iter%tn2(j  )
             write(1,*) iter%tn1(i+1),iter%tn2(j  )
             write(1,*) iter%tn1(i+1),iter%tn2(j+1)
             write(1,*) iter%tn1(i  ),iter%tn2(j+1)
             write(1,*) iter%tn1(i  ),iter%tn2(j  )
             write(1,*)
             write(1,*)
          end do
       end do
       iter=>iter%next
    end do
    close(1)
    
  end subroutine draw_LRmesh

  !> draw surface
  !> \brief gnuplotでp file
  !> \param n (n+1)x(n+1)のgridでsurfaceを描画
  subroutine draw_surface(a,n,fn)
    type(LocallyRefinedSpline),pointer,intent(in) :: a
    integer,intent(in) :: n
    character(*),intent(in) :: fn

    type(LocallyRefinedSpline),pointer :: iter    
    integer :: i, j
    real(8) :: t1, t2, b1, b2, c(ndim+1)
    
    open(1,file=fn)
    do j=0,n
       t2=dble(j)/dble(n)
       do i=0,n
          t1=dble(i)/dble(n)
          c(:)=0.d0
          iter=>a
          do while(associated(iter))
             b1=basis_func(ndeg(1),iter%tn1,t1)
             if(abs(b1).le.tiny) goto 1
             b2=basis_func(ndeg(2),iter%tn2,t2)
             if(abs(b2).le.tiny) goto 1
             c(:)=c(:)+b1*b2*iter%cp(:)
1            iter=>iter%next
          end do
          write(1,*) c(:)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine draw_surface
  
end module module_LRspline
