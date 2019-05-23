      subroutine output2_bin(nm1,sden)
c write binary sensrx file.
      use COMS
      implicit none
      integer nm1
      real*8 sden(neg,nm1)
      character hname*8,huse(2)*8
      integer i,ng,nwd
      parameter (nwd=4)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      write(afile,'("sensrx")')
      write(hname,'("sensrx")')
      write(*,'("writing ",a,"...")')trim(afile)
      write(iuo,'("writing ",a,"...")')trim(afile)
c
c 5.97 opens this file with form='binary',access='sequential' but this
c still works.
      open(iub,file=afile,status='unknown',position='rewind',
#if ( defined INTEL )
     1 form='unformatted',access='stream'
#elif ( defined GNU )
     1 access='stream'  ! untested capability
#endif
     2 )
      ia(1)=neg
      ia(2)=nr
      ia(3)=nz
      ia(4)=nm1
      ivers=0
      huse(1:2)=" "
      write(iub)hname,(huse(i),i=1,2),ivers
      write(iub)(ia(i),i=1,nwd)
      do i=1,nm1
        do ng=1,neg
          write(iub)sden(ng,i)
        end do ! ng
      end do ! i
c     close(iub) ! keep open for surfaces
c
      return
      end subroutine
c
      subroutine rdoutput2_bin(nm1,sden)
c read binary sensrx file.
      use COMS
      implicit none
      integer nm1
      real*8 sden(neg,nm1)
      character hname*8,huse(2)*8
      integer i,ng,nwd
      parameter (nwd=4)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      write(afile,'("sensrx")')
      write(hname,'("sensrx")')
      write(*,'("reading ",a,"...")')trim(afile)
      write(iuo,'("reading ",a,"...")')trim(afile)
c
c 5.97 opens this file with form='binary',access='sequential' but this
c still works.
      open(iub,file=afile,status='old',position='rewind',
#if ( defined INTEL )
     1 form='unformatted',access='stream',
#elif ( defined GNU )
     1 access='stream',  ! untested capability
#endif
     2 err=900)
      read(iub,err=910)hname,(huse(i),i=1,2),ivers
      read(iub,err=910)(ia(i),i=1,nwd)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(10i8)')(ia(i),i=1,nwd)
      endif
      if(ia(1).ne.neg)go to 800
      if(ia(2).ne.nr)go to 805
      if(ia(3).ne.nz)go to 810
      if(ia(4).ne.nm1)go to 815
c
      do i=1,nm1
        do ng=1,neg
          read(iub)sden(ng,i)
        end do ! ng
      end do ! i
c     close(iub) ! keep open for surfaces
      write(*,'("done but still open.")')
      go to 999
c
  800 write(*,'("error. ia(1).ne.neg in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      write(*,'("neg=",i8)')neg
      call stoponerror
  805 write(*,'("error. ia(2).ne.nr in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("nr=",i8)')nr
      call stoponerror
  810 write(*,'("error. ia(3).ne.nz in ",a,"; ia(3)=",i8)')
     1 trim(afile),ia(3)
      write(*,'("nz=",i8)')nz
      call stoponerror
  815 write(*,'("error. ia(4).ne.nm1 in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      write(*,'("nm1=",i8)')nm1
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
