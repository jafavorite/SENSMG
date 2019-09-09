      subroutine wrawbin
c write binary file for atomic weights, densities, weight function,
c and sigsp1.
      use COMS
      implicit none
      character hname*8,huse(2)*8
      integer i,j,lg,ng,nwd
      parameter (nwd=5)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      write(afile,'("sensaw")')
      write(hname,'("sensaw")')
      write(*,'("writing ",a,"...")')trim(afile)
      write(iuo,'("writing ",a,"...")')trim(afile)
c
c 5.97 opens this file with form='binary',access='sequential' but this
c still works.
      open(iub1,file=afile,status='unknown',position='rewind',
#if ( defined INTEL )
     1 form='unformatted',access='stream'
#elif ( defined GNU )
     1 access='stream'  ! untested capability
#endif
     2 )
      ia(1)=neg
      ia(2)=nm
      ia(3)=nel
      ia(4)=niso
      if(itrcor.eq.1.or.itrcor.eq.2)then
        ia(5)=1
      else
        ia(5)=0
      end if
      ivers=0
      huse(1:2)=" "
      write(iub1)hname,(huse(i),i=1,2),ivers
      write(iub1)(ia(i),i=1,nwd)
      write(iub1)(atwt(i),i=1,niso)
      write(iub1)((blk(j,i),i=1,nel),j=2,3)
      write(iub1)(rho(i),i=0,nm)
      write(iub1)(rhoa(i),i=0,nm)
      write(iub1)(zaidfull(i),i=1,niso)
      write(iub1)((wgtfct(ng,i),i=1,niso),ng=1,neg)
      if(ia(5).eq.1)write(iub1)(((sigsp1(lg,ng,nm+i),i=1,niso),
     1 ng=1,neg),lg=1,neg)
      close(iub1)
      return
      end subroutine
c
      subroutine rdawbin
c read binary file for atomic weights, densities, weight function,
c and sigsp1.
      use COMS
      implicit none
      character hname*8,huse(2)*8
      integer i,j,lg,ng,nwd
      parameter (nwd=5)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      write(afile,'("sensaw")')
      write(hname,'("sensaw")')
      write(*,'("reading ",a,"...")')trim(afile)
      write(iuo,'("reading ",a,"...")')trim(afile)
c
c 5.97 opens this file with form='binary',access='sequential' but this
c still works.
      open(iub1,file=afile,status='old',position='rewind',
#if ( defined INTEL )
     1 form='unformatted',access='stream',
#elif ( defined GNU )
     1 access='stream',  ! untested capability
#endif
     2 err=900)
      read(iub1,err=910)hname,(huse(i),i=1,2),ivers
      read(iub1,err=910)(ia(i),i=1,nwd)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(10i8)')(ia(i),i=1,nwd)
      endif
      if(ia(1).ne.neg)go to 800
      if(ia(2).ne.nm)go to 805
      if(ia(3).ne.nel)go to 810
      if(ia(4).ne.niso)go to 815
c
      read(iub1)(atwt(i),i=1,niso)
      read(iub1)((blk(j,i),i=1,nel),j=2,3)
      read(iub1)(rho(i),i=0,nm)
      read(iub1)(rhoa(i),i=0,nm)
      read(iub1)(zaidfull(i),i=1,niso)
      read(iub1)((wgtfct(ng,i),i=1,niso),ng=1,neg)
      if(ia(5).eq.1)read(iub1)(((sigsp1(lg,ng,nm+i),i=1,niso),
     1 ng=1,neg),lg=1,neg)
      close(iub1)
c
      go to 999
c
  800 write(*,'("error. ia(1).ne.neg in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      write(*,'("neg=",i8)')neg
      call stoponerror
  805 write(*,'("error. ia(2).ne.nm in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("nm=",i8)')nm
      call stoponerror
  810 write(*,'("error. ia(3).ne.nel in ",a,"; ia(3)=",i8)')
     1 trim(afile),ia(3)
      write(*,'("nel=",i8)')nel
      call stoponerror
  815 write(*,'("error. ia(4).ne.niso in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      write(*,'("niso=",i8)')niso
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
