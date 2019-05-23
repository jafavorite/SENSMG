      subroutine wrfixsrc(ig,ib,nl,gams)
c write binary fixsrc file.
      use COMS
      implicit none
      integer ig,ib,nl
      real*8 gams(neg,0:nmom-1,it,jt)
      character hname*8,huse(2)*8
      integer i,j,k,l,ng,nwd0,nwd
      parameter (nwd=10)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      if(ig.eq.1.or.ig.eq.2)then
        write(afile,'("a",i2.2,"/a",i2.2,"_fixsrc")')ib,ib
        write(hname,'("a",i2.2)')ib
      else if(ig.eq.3)then
        if(ib.eq.1.or.ib.eq.2)then
          write(afile,'("smf/",i2.2,"_fixsrc")')ib
        else if(ib.eq.3)then
          write(afile,'("sma/",i2.2,"_fixsrc")')ib
        end if
        write(hname,'(i2.2,"fixsrc")')ib
      end if
      write(*,'("writing ",a,"...")')trim(afile)
      write(iuo,'("writing ",a,"...")')trim(afile)
c
c 5.97 opens this file with form='binary',access='sequential' but this
c still works.
      open(iuf,file=afile,status='unknown',position='rewind',
#if ( defined INTEL )
     1 form='unformatted',access='stream'
#elif ( defined GNU )
     1 access='stream'  ! untested capability
#endif
     2 )
      if(calc2d.eq.0.or.calc2d.eq.-1)then
        ia(1)=1
      else if(calc2d.eq.1)then
        ia(1)=2
      end if
      ia(2)=neg
      ia(3)=it
      ia(4)=jt
      ia(5)=1
      ia(6)=0 ! unused
      ia(7)=0 ! unused
      ia(9)=nl+1 ! nmq is indexed from 0
      if(partisnver(1:1).eq."5")then
        ia(8)=0 ! unused
        nwd0=9
      else
        ia(8)=lng
        ia(10)=1 ! flag for group-dependent data
        nwd0=10
      end if
      ivers=0
      huse(1:2)=" "
      write(iuf)hname,(huse(i),i=1,2),ivers
      write(iuf)(ia(i),i=1,nwd0)
      do ng=1,neg
        do l=0,nl
          write(iuf)(((gams(ng,l,i,j),i=1,ia(3)),j=1,ia(4)),
     1     k=1,ia(5))
        end do ! l
      end do ! ng
      close(iuf)
c debug print
c     if(ib.eq.1)then
c       write(iuo,'("debug print, iteration",i3)')iter
c       l=0
c       do i=1,it
c         write(100+iter,'(i6,1p100e14.6)')i,0.5d0*(rfm(i)+rfm(i-1)),
c    1     gams(1:neg,l,i,1),sum(gams(1:neg,l,i,1))
c       end do ! i
c     end if
c
      return
      end subroutine
