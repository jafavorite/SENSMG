      subroutine wrgmom(ib,iconv)
c write binary moments of sum of generalized adjoints.
c the moments are in mesh centers.
      use COMS
      implicit none
      integer ib,iconv
      character hname*8
      integer*4 ia(8)
      integer i,j,k,l,ng
      character afile*20
c
      write(afile,'("a",i2.2,"/a",i2.2,"_mom")')ib,ib
      write(*,'("writing ",a,"...")')trim(afile)
      write(iuo,'("writing ",a,"...")')trim(afile)
c
      open(iuga,file=afile,status='unknown',position='rewind',
#if ( defined INTEL )
     1 form='binary'
#elif ( defined GNU )
     1 access='stream'  ! untested capability
#endif
     2 )
      write(hname,'("a",i2.2)')ib
      if(calc2d.eq.0.or.calc2d.eq.-1)then
        ia(1)=1
      else if(calc2d.eq.1)then
        ia(1)=2
      end if
      ia(2)=neg
      ia(3)=it
      ia(4)=jt
      ia(5)=1
      ia(6)=nmom
      ia(7)=iter
      ia(8)=iconv
      write(iuga)hname
      write(iuga)(ia(i),i=1,8)
      do ng=1,neg
        write(iuga)(((gmom(ng,0,i,j),i=1,ia(3)),j=1,ia(4)),
     1   k=1,ia(5))
        if(ia(6).gt.1)write(iuga)((((gmom(ng,l,i,j),i=1,ia(3)),
     1   j=1,ia(4)),k=1,ia(5)),l=1,ia(6)-1)
      end do ! ng
      close(iuga)
c debug print
c     if(ib.eq.1)then
c       write(iuo,'("debug print, iteration",i3)')iter
c       l=0
c       do i=1,it
c         write(100+iter,'(i6,1p100e14.6)')i,0.5d0*(rfm(i)+rfm(i-1)),
c    1     gmom(1:neg,l,i,1),sum(gmom(1:neg,l,i,1))
c       end do ! i
c     end if
c
      return
      end subroutine
