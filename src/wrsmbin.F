      subroutine wrsmbin(r2,ipsrc,ipfis,ipsrc_s,ipfis_s)
c write binary file for feyny or sm2 data.
      use COMS
      implicit none
      real*8 r2,ipsrc(it),ipfis(it),ipsrc_s(0:nr),ipfis_s(0:nr)
      character hname*8,huse(2)*8
      integer i,ng,nwd
      parameter (nwd=7)
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
      ia(3)=nr
      ia(4)=nz
      ia(5)=it
      ia(6)=jt
      ia(7)=nm
      ivers=0
      huse(1:2)=" "
      write(iub)hname,(huse(i),i=1,2),ivers
      write(iub)(ia(i),i=1,nwd)
      write(iub)yasym,sm2,r2,lkg
      write(iub)(ipsrc(i),i=1,ia(5))
      write(iub)(ipfis(i),i=1,ia(5))
      write(iub)((chisrc(ng,i),i=0,ia(7)),ng=1,ia(2))
      write(iub)(ipsrc_s(i),i=0,ia(3))
      write(iub)(ipfis_s(i),i=0,ia(3))
      close(iub)
      return
      end subroutine
c
      subroutine rdsmbin(r2,ipsrc,ipfis,ipsrc_s,ipfis_s)
c read binary file for feyny or sm2 data.
      use COMS
      implicit none
      real*8 r2,ipsrc(it),ipfis(it),ipsrc_s(0:nr),ipfis_s(0:nr)
      character hname*8,huse(2)*8
      integer i,ng,nwd
      parameter (nwd=7)
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
      if(ia(3).ne.nr)go to 805
      if(ia(4).ne.nz)go to 810
      if(ia(5).ne.it)go to 815
      if(ia(6).ne.jt)go to 820
      if(ia(7).ne.nm)go to 825
c
      read(iub)yasym,sm2,r2,lkg
      read(iub)(ipsrc(i),i=1,ia(5))
      read(iub)(ipfis(i),i=1,ia(5))
      read(iub)((chisrc(ng,i),i=0,ia(7)),ng=1,ia(2))
      read(iub)(ipsrc_s(i),i=0,ia(3))
      read(iub)(ipfis_s(i),i=0,ia(3))
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
      write(*,'("neg=",i8)')neg
      call stoponerror
  805 write(*,'("error. ia(3).ne.nr in ",a,"; ia(3)=",i8)')
     1 trim(afile),ia(3)
      write(*,'("nr=",i8)')nr
      call stoponerror
  810 write(*,'("error. ia(4).ne.nz in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      write(*,'("nz=",i8)')nz
      call stoponerror
  815 write(*,'("error. ia(5).ne.it in ",a,"; ia(5)=",i8)')
     1 trim(afile),ia(5)
      write(*,'("it=",i8)')it
      call stoponerror
  820 write(*,'("error. ia(6).ne.jt in ",a,"; ia(6)=",i8)')
     1 trim(afile),ia(6)
      write(*,'("jt=",i8)')jt
      call stoponerror
  825 write(*,'("error. ia(7).ne.nm in ",a,"; ia(7)=",i8)')
     1 trim(afile),ia(7)
      write(*,'("nm=",i8)')nm
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
