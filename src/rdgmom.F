      subroutine rdgmom(ib,iconv)
c read binary moments of sum of generalized adjoints.
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
      read(iuga)(ia(i),i=1,8)
      do ng=1,neg
        read(iuga)(((gmom(ng,0,i,j),i=1,ia(3)),j=1,ia(4)),
     1   k=1,ia(5))
        if(ia(6).gt.1)read(iuga)((((gmom(ng,l,i,j),i=1,ia(3)),
     1   j=1,ia(4)),k=1,ia(5)),l=1,ia(6)-1)
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
