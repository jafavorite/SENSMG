      subroutine wrsmbin(r2,mrsrc,mrfis,ipsrc,ipfis,chivec)
c write binary file for feyny or sm2 data.
      use COMS
      implicit none
      real*8 r2,mrsrc(0:nm),mrfis(0:nm),ipsrc(it),ipfis(it),
     1 chivec(neg,nm)
      character hname*8,huse(2)*8
      integer i,ng,nwd
      parameter (nwd=5)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      write(afile,'("senssm")')
      write(hname,'("senssm")')
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
      if(calc2d.eq.0.or.calc2d.eq.-1)then
        ia(1)=1
      else if(calc2d.eq.1)then
        ia(1)=2
      end if
      ia(2)=neg
      ia(3)=it
      ia(4)=jt
      ia(5)=nm
      ivers=0
      huse(1:2)=" "
      write(iub)hname,(huse(i),i=1,2),ivers
      write(iub)(ia(i),i=1,nwd)
      write(iub)yasym,sm2,r2,lkg
      write(iub)(mrsrc(i),i=0,ia(5)),(mrfis(i),i=0,ia(5))
      write(iub)(ipsrc(i),i=1,ia(3)),(ipfis(i),i=1,ia(3))
      do ng=1,neg
        write(iub)(chivec(ng,i),i=1,ia(5))
      end do ! ng
      close(iub)
      return
      end subroutine
c
      subroutine rdsmbin(r2,mrsrc,mrfis,ipsrc,ipfis,chivec)
c read binary file for feyny or sm2 data.
      use COMS
      implicit none
      real*8 r2,mrsrc(0:nm),mrfis(0:nm),ipsrc(it),ipfis(it),
     1 chivec(neg,nm)
      character hname*8,huse(2)*8
      integer i,ng,nwd
      parameter (nwd=5)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      write(afile,'("senssm")')
      write(hname,'("senssm")')
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
      if((calc2d.eq.0.or.calc2d.eq.-1).and.ia(1).ne.1)go to 795
      if(calc2d.eq.1.and.ia(1).ne.2)go to 797
      if(ia(2).ne.neg)go to 800
      if(ia(3).ne.it)go to 805
      if(ia(4).ne.jt)go to 810
      if(ia(5).ne.nm)go to 815
c
      read(iub)yasym,sm2,r2,lkg
      read(iub)(mrsrc(i),i=0,ia(5)),(mrfis(i),i=0,ia(5))
      read(iub)(ipsrc(i),i=1,ia(3)),(ipfis(i),i=1,ia(3))
      do ng=1,neg
        read(iub)(chivec(ng,i),i=1,ia(5))
      end do ! ng
      close(iub)
      go to 999
c
  795 write(*,'("error. ia(1).ne.1 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  797 write(*,'("error. ia(1).ne.2 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  800 write(*,'("error. ia(2).ne.neg in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      call stoponerror
  805 write(*,'("error. ia(3).ne.it in ",a,"; ia(3)=",i8)')
     1 trim(afile),ia(3)
      call stoponerror
  810 write(*,'("error. ia(4).ne.jt in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      call stoponerror
  815 write(*,'("error. ia(5).ne.nm in ",a,"; ia(5)=",i8)')
     1 trim(afile),ia(5)
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
