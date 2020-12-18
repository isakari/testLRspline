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
  subroutine init_meshline(a)
    type(MeshLine),pointer,intent(out) :: a !< MeshLine to be initalised
    nullify(a)
  end subroutine init_meshline

  !> uninitialise meshline
  subroutine uninit_meshline(a)
    type(MeshLine),pointer,intent(inout) :: a !< MeshLine to be destroyed

    type(MeshLine),pointer :: iter, tmp

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
  subroutine append_meshline(a,idir,pos,st_,en_)
    type(MeshLine),pointer,intent(inout) :: a !< meshline
    integer,intent(in) :: idir !< 1=horizontal, 2=vertical meshline
    real(8),intent(in) :: pos !< eta=pos (idir=1), xi=pos (idir=2)
    real(8),intent(in) :: st_ !< (st,pos)-->(en,pos) if idir=1
    real(8),intent(in) :: en_ !< (pos,st)-->(pos,en) if idir=2
    
    type(MeshLine),pointer :: iter
    integer :: icnt, idel(2)
    real(8) :: st, en
    st=st_
    en=en_

    if(associated(a))then !aがあれば
       iter=>a
       do while(associated(iter))
          if((idir.eq.iter%idir).and.(abs(pos-iter%pos).le.tiny)& !mlがiterに含まれるなら
               .and.(iter%st.le.st).and.(en.le.iter%en)) return ! 何もしない
          iter=>iter%next
       end do
       idel(:)=-1
       icnt=0
       iter=>a
       do while(associated(iter))
          if((idir.eq.iter%idir)& !追加する予定のメッシュライン(ml)が
               .and.abs(pos-iter%pos).le.tiny)then !既存のmeshline(iter)と同一直線上にあり
             if((st.le.iter%st)&
                  !.and.(iter%st.le.en)&
                  .and.(en.le.iter%en))then !st--iter%st--en--iter%enのとき
                en=iter%en !mlの終点はiterの終点とする
                idel(1)=icnt !で、iterは消す
             elseif((iter%st.le.st)&
                  !.and.(st.le.iter%en)&
                  .and.(iter%en.le.en))then !iter%st--st--iter%en--enのとき
                st=iter%st !mlの始点はiterの始点とする
                idel(2)=icnt !で、iterは消す
             end if
             if(idel(1).ne.-1.and.idel(2).ne.-1) goto 1
          end if
          icnt=icnt+1
          iter=>iter%next
       end do
1      continue
       iter=>tail(a) !aのお尻
       call set_new_meshline(iter%next,idir,pos,st,en) !の次にmeshlineをset
       if(maxval(idel).ge.0) call remove_meshline(a,maxval(idel)) !消すべきものを
       if(minval(idel).ge.0) call remove_meshline(a,minval(idel)) !後ろから順に消す
    else !aがなければ
       call set_new_meshline(a,idir,pos,st,en) !a(の先頭)にmeshlineをset
    end if

  end subroutine append_meshline
  
end module module_meshline
