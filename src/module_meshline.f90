module module_meshline
  use module_globals
  implicit none

  private
  public MeshLine, &
       init_meshline, uninit_meshline, &
       append_meshline

  !> mesh line
  !> \brief if idir=1, meshline is (st,pos)-->(en,pos)
  !> \brief if idir=2, meshline is (pos,st)-->(pos,en)
  type MeshLine
     integer :: idir !< meshlineの方向 1(=x) or 2(=y)
     real(8) :: pos !< meshlineの位置
     real(8) :: st !< meshlineの始点
     real(8) :: en !< meshlineの終点
     type(MeshLine),pointer :: next
  end type MeshLine

contains

  !> initalise meshline
  subroutine init_meshline(ml_)
    type(MeshLine),pointer,intent(out) :: ml_ !< MeshLine to be initalised
    nullify(ml_)
  end subroutine init_meshline

  !> uninitialise meshline
  subroutine uninit_meshline(ml_)
    type(MeshLine),pointer,intent(inout) :: ml_ !< MeshLine to be destroyed

    type(MeshLine),pointer :: ml, tmp

    if(.not.associated(ml_))then !ml_がなければ
       return !なにもしない
    end if

    ml=>ml_ !ml_の先頭へのpointer
    do while(associated(ml))
       tmp=>ml !=mlの現在の先頭へのpointer
       ml=>ml%next !=mlの現在の先頭の次
       deallocate(tmp) !mlの現在の先頭をdeallocate
    end do

    nullify(ml_)
    
  end subroutine uninit_meshline

  !> return the last meshline
  function tail(a) result(b)
    type(MeshLine),pointer,intent(in) :: a
    type(MeshLine),pointer :: b
    
    b=>a !bでaの先頭を指し
    do while(associated(b%next)) !bの次がnullになるまで
       b=>b%next !bを辿っていく
    end do
    
  end function tail
  
  !> new_meshline(とそのnext)をset
  subroutine set_new_meshline(new_ml,idir,pos,st,en,next)
    type(MeshLine),pointer,intent(inout) :: new_ml !< setするmeshline
    integer,intent(in) :: idir !< 1=horizontal, 2=vertical meshline
    real(8),intent(in) :: pos !< eta=pos (idir=1), xi=pos (idir=2)
    real(8),intent(in) :: st !< (st,pos)-->(en,pos) if idir=1
    real(8),intent(in) :: en !< (pos,st)-->(pos,en) if idir=2
    type(MeshLine),pointer,optional,intent(in) :: next

    nullify(new_ml);
    allocate(new_ml) !空の領域をallocate

    if(present(next))then !nextが指定されていれば
       new_ml%next=>next !new_mlの次にnextを渡す
    else !そうでなければ
       nullify(new_ml%next) !new_mlの次は何もない
    end if

    !与えられたmeshlineの情報をnew_mlに渡す
    new_ml%idir=idir
    new_ml%pos=pos
    new_ml%st=st
    new_ml%en=en

  end subroutine set_new_meshline

  !> check if (idir,pos,st,en) is included in ml
  function is_included(ml,idir,pos,st,en) result(res)
    type(MeshLine),pointer,intent(in) :: ml !< existing meshline
    integer,intent(in) :: idir !< direction of the meshline to be inserted (ml2)
    real(8),intent(in) :: pos !< position of the meshline to be inserted
    real(8),intent(in) :: st !< starting point of the meshline to be inserted
    real(8),intent(in) :: en !< end point of the meshline to be inserted
    logical :: res !< T if ml2 included in ml, F otherwise

    res=(idir.eq.ml%idir)& ! 同じ方向を向いていて
         .and.(abs(pos-ml%pos).le.tiny)& ! 位置が同じで
         .and.(ml%st.le.st).and.(en.le.ml%en) ! ml2 \subset ml ならば T
    
  end function is_included

  !> check if (idir,pos,st,en) is on the same line as ml
  function is_on_the_same_line(ml,idir,pos) result(res)
    type(MeshLine),pointer,intent(in) :: ml !< existing meshline
    integer,intent(in) :: idir !< direction of the meshline to be inserted (ml2)
    real(8),intent(in) :: pos !< position of the meshline to be inserted
    logical :: res !< T if ml2 included in ml, F otherwise

    res=(idir.eq.ml%idir)& ! 同じ方向を向いていて
         .and.(abs(pos-ml%pos).le.tiny) ! 位置が同じならば T
    
  end function is_on_the_same_line

  function is_elongating_to_left(ml,st,en) result(res)
    type(MeshLine),pointer,intent(in) :: ml !< existing meshline
    real(8),intent(in) :: st
    real(8),intent(in) :: en
    logical :: res

    res=(st.le.ml%st).and.(ml%st.le.en).and.(en.le.ml%en)

  end function is_elongating_to_left

  function is_elongating_to_right(ml,st,en) result(res)
    type(MeshLine),pointer,intent(in) :: ml !< existing meshline
    real(8),intent(in) :: st
    real(8),intent(in) :: en
    logical :: res

    res=(ml%st.le.st).and.(st.le.ml%en).and.(ml%en.le.en)

  end function is_elongating_to_right

  !> remove a meshline at the given index
  subroutine remove_meshline(a,idx)
    type(MeshLine),pointer,intent(inout) :: a
    integer,intent(in) :: idx

    integer :: i
    type(MeshLine), pointer :: iter, prev

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

    if(.not.associated(prev))then ! first node
       a=>a%next
    else
       prev%next=>iter%next
    end if    
    deallocate(iter)

  end subroutine remove_meshline
  
  !> append a meshline at the last 
  subroutine append_meshline(ml_,idir,pos,st_,en_)
    type(MeshLine),pointer,intent(inout) :: ml_ !< meshline
    integer,intent(in) :: idir !< 1=horizontal, 2=vertical meshline
    real(8),intent(in) :: pos !< eta=pos (idir=1), xi=pos (idir=2)
    real(8),intent(in) :: st_ !< (st,pos)-->(en,pos) if idir=1
    real(8),intent(in) :: en_ !< (pos,st)-->(pos,en) if idir=2

    integer :: idel(2), icnt
    real(8) :: st, en
    type(MeshLine),pointer :: ml
    st=st_
    en=en_

    if(associated(ml_))then !ml_があれば
       ml=>ml_
       do while(associated(ml))
          if(is_included(ml,idir,pos,st,en)) return !(idir,pos,st,en)\subset mlなら何もしない
          ml=>ml%next
       end do
       idel(:)=-1
       icnt=0
       ml=>ml_
       do while(associated(ml))
          if(is_on_the_same_line(ml,idir,pos)) then !追加する予定のメッシュラインが既存のmeshline(ml)と同一直線上にあり、
             if(is_elongating_to_left(ml,st,en)) then !st--ml%st--en--ml%enのとき
                en=ml%en !mlの終点はmlの終点とする
                idel(1)=icnt !で、mlは消す
             elseif(is_elongating_to_right(ml,st,en)) then !ml%st--st--ml%en--enのとき
                st=ml%st !mlの始点はmlの始点とする
                idel(2)=icnt !で、mlは消す
             end if
             if(idel(1).ne.-1.and.idel(2).ne.-1) goto 1
          end if
          icnt=icnt+1
          ml=>ml%next
       end do
1      continue
       ml=>tail(ml_) !ml_のお尻
       call set_new_meshline(ml%next,idir,pos,st,en) !の次にmeshlineをset
       if(maxval(idel).ge.0) call remove_meshline(ml_,maxval(idel)) !消すべきものを
       if(minval(idel).ge.0) call remove_meshline(ml_,minval(idel)) !後ろから順に消す
    else !ml_がなければ
       call set_new_meshline(ml_,idir,pos,st,en) !ml_(の先頭)にmeshlineをset
    end if

  end subroutine append_meshline
  
end module module_meshline
