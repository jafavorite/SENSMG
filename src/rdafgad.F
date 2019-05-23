      subroutine rdafgad(ib,iconv)
c read binary angular fluxes of sum of generalized adjoints.
c the angular fluxes are on mesh edges.
      use COMS
      implicit none
      integer ib,iconv
      character hname*8
      integer nwd
      parameter (nwd=8)
      integer*4 ia(nwd)
      integer i,m,ng
      character afile*20
c
      write(afile,'("a",i2.2,"/a",i2.2,"_afgad")')ib,ib
      write(*,'("reading ",a,"...")')trim(afile)
      write(iuo,'("reading ",a,"...")')trim(afile)
c
      open(iuga,file=afile,status='old',position='rewind',
#if ( defined INTEL )
     1 form='binary',
#elif ( defined GNU )
     1 access='stream', ! untested capability
#endif
     2 err=900)
      read(iuga)hname
      read(iuga)(ia(i),i=1,nwd)
      do ng=1,neg
        do m=1,ndir
          read(iuga)(afgad(ng,m,i,1),i=1,ia(3)+1)
        end do ! m
      end do ! ng
      close(iuga)
      iconv=ia(8)
c
      go to 999
c
  900 write(*,'("error opening ",a,".")')trim(afile)
      write(iuo,'("error opening ",a,".")')trim(afile)
      call stoponerror
  999 return
      end subroutine
