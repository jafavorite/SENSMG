      subroutine wrafgad(ib,iconv)
c write binary angular fluxes of sum of generalized adjoints.
c the angular fluxes are on mesh edges.
      use COMS
      implicit none
      integer ib,iconv
      character hname*8
      integer i,m,ng,nwd
      parameter (nwd=8)
      integer*4 ia(nwd)
      character afile*20
c
      write(afile,'("a",i2.2,"/a",i2.2,"_afgad")')ib,ib
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
      ia(5)=1 ! would be kt
      ia(6)=ndir
      ia(7)=iter
      ia(8)=iconv
      write(iuga)hname
      write(iuga)(ia(i),i=1,nwd)
      do ng=1,neg
        do m=1,ndir
          write(iuga)(afgad(ng,m,i,1),i=1,ia(3)+1)
        end do ! m 
      end do ! ng
      close(iuga)
c
      return
      end subroutine
